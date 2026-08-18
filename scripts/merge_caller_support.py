#!/usr/bin/env python3
import argparse
import glob
import gzip
import os
import re
from collections import Counter, OrderedDict
from itertools import combinations


CALLER_PREFIX = {
    "mutect2": "M2",
    "deepsomatic": "DS",
    "clairsto": "CT",
}


# Which FILTER values let a caller cast a vote under --require-pass.
#
# A per-caller set rather than a bare `== "PASS"` test, because the three callers use
# different FILTER vocabularies (Mutect2: germline/weak_evidence/strand_bias/...,
# DeepSomatic: GERMLINE/RefCall/NoCall/PON, ClairS-TO: NonSomatic/LowQual/...), and a
# future caller may have a second label that legitimately counts as somatic-PASS.
CALLER_PASS = {
    "mutect2": {"PASS"},
    "deepsomatic": {"PASS"},
    "clairsto": {"PASS"},
}


EVIDENCE_SUFFIXES = ("FILTER", "QUAL", "AF", "DP", "AD")


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(info_text):
    info = OrderedDict()
    if not info_text or info_text == ".":
        return info
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def format_info(info):
    parts = []
    for key, value in info.items():
        if value is True:
            parts.append(key)
        elif value is None or value == "":
            continue
        else:
            parts.append(f"{key}={sanitize_info_value(str(value))}")
    return ";".join(parts) if parts else "."


def sanitize_info_value(value):
    return (
        value.replace(" ", "_")
        .replace(";", ",")
        .replace("\t", "_")
        .replace("\n", "_")
    )


def header_id(line):
    match = re.match(r"##(INFO|FORMAT|FILTER|contig)=<ID=([^,>]+)", line)
    if match:
        return match.group(1), match.group(2)
    return None


class Record:
    def __init__(self, caller, fields, line_no):
        self.caller = caller
        self.fields = fields
        self.line_no = line_no
        self.chrom = fields[0]
        self.pos = fields[1]
        self.ref = fields[3]
        self.alt = fields[4]

    @property
    def key(self):
        return self.chrom, int(self.pos), self.ref, self.alt

    @property
    def qual(self):
        return self.fields[5] if len(self.fields) > 5 else "."

    @property
    def filt(self):
        return self.fields[6] if len(self.fields) > 6 else "."

    @property
    def info(self):
        return parse_info(self.fields[7] if len(self.fields) > 7 else ".")

    def sample_values(self):
        if len(self.fields) < 10:
            return {}, {}
        fmt_keys = self.fields[8].split(":")
        fmt_vals = self.fields[9].split(":")
        values = OrderedDict(zip(fmt_keys, fmt_vals))
        return values, self.info

    def evidence(self):
        fmt, info = self.sample_values()
        af = first_value(fmt, ("AF", "VAF", "T_AF", "TUMOR_AF"))
        dp = first_value(fmt, ("DP", "T_DP"))
        ad = first_value(fmt, ("AD", "T_AD"))

        if not af:
            af = first_value(info, ("AF", "VAF", "TLOD_AF", "ALLELE_FRACTION"))
        if not dp:
            dp = first_value(info, ("DP", "TDP"))
        if not ad:
            ad = first_value(info, ("AD",))

        if (not af or not dp) and ad:
            calc_dp, calc_af = calc_from_ad(ad)
            if not dp:
                dp = calc_dp
            if not af:
                af = calc_af

        return {
            "FILTER": self.filt,
            "QUAL": self.qual,
            "AF": af or ".",
            "DP": dp or ".",
            "AD": ad or ".",
        }


def first_value(container, keys):
    for key in keys:
        if key in container and container[key] not in ("", "."):
            return container[key]
    return ""


def calc_from_ad(ad_value):
    try:
        parts = [float(x) for x in ad_value.split(",") if x not in ("", ".")]
    except ValueError:
        return "", ""
    if len(parts) < 2:
        return "", ""
    total = sum(parts)
    if total <= 0:
        return str(int(total)), ""
    return str(int(total)), f"{parts[1] / total:.6g}"


def read_vcf(caller, path):
    meta = []
    chrom_header = None
    records = []
    contigs = []
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if line.startswith("##"):
                meta.append(line)
                if line.startswith("##contig=<ID="):
                    hid = header_id(line)
                    if hid:
                        contigs.append(hid[1])
                continue
            if line.startswith("#CHROM"):
                chrom_header = line
                continue
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_no} has fewer than 8 VCF columns")
            for alt in fields[4].split(","):
                split_fields = fields[:]
                split_fields[4] = alt
                records.append(Record(caller, split_fields, line_no))
    if not chrom_header:
        raise ValueError(f"{path} is missing a #CHROM header")
    return meta, chrom_header, contigs, records


def merge_headers(vcf_parts, sample):
    seen_id = set()
    seen_exact = set()
    merged = []
    contig_order = OrderedDict()
    for meta, _, contigs, _ in vcf_parts:
        for contig in contigs:
            contig_order.setdefault(contig, len(contig_order))
        for line in meta:
            key = header_id(line)
            if key:
                if key in seen_id:
                    continue
                seen_id.add(key)
            elif line in seen_exact:
                continue
            seen_exact.add(line)
            merged.append(line)

    merged.append("##source=merge_caller_support.py")
    merged.append('##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of tumor-only callers supporting this normalized allele">')
    merged.append('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Tumor-only callers supporting this normalized allele">')
    merged.append('##INFO=<ID=SELECTED_CALLER,Number=1,Type=String,Description="Caller record selected as the representative VCF row">')
    for prefix in CALLER_PREFIX.values():
        for suffix in EVIDENCE_SUFFIXES:
            merged.append(f'##INFO=<ID={prefix}_{suffix},Number=.,Type=String,Description="{prefix} caller {suffix} evidence before consensus merge">')

    chrom_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if sample:
        chrom_header += f"\tFORMAT\t{sample}"
    return merged, chrom_header, contig_order


def build_index(records):
    by_key = OrderedDict()
    for record in records:
        by_key.setdefault(record.key, OrderedDict())[record.caller] = record
    return by_key


def somatic_index(by_key):
    """Drop non-PASS caller records so only somatic-PASS callers can vote.

    Why this is needed: harmonize_filter_vcf.py gates on AF / ALT-read-count / DP and
    never reads the FILTER column, so GERMLINE, NonSomatic, RefCall, NoCall and PON
    records arrive here untouched. Germline SNPs sit at AF~0.5 with high depth, clear
    every threshold, and DeepSomatic and ClairS-TO each emit them independently. The
    result was a concordance table reporting support_2 (48,598) ABOVE support_1
    (39,468) for S1-DA-01 -- impossible for a real 2-caller intersection -- because
    that one caller pair contributed 43,193 germline sites. 42,187 of the 57,436
    "concordant" records for that sample were FILTER=GERMLINE against 2,369 PASS.

    Returns an index of the same shape as build_index(), so every helper below works
    on it unchanged. Sites left with no supporting caller are dropped entirely: they
    are not "support_0", they are simply not somatic candidates.
    """
    gated = OrderedDict()
    for key, callers in by_key.items():
        keep = OrderedDict()
        for caller, record in callers.items():
            if record.filt in CALLER_PASS.get(caller, {"PASS"}):
                keep[caller] = record
        if keep:
            gated[key] = keep
    return gated


def load_candidate_keys(candidate_dir):
    if not candidate_dir:
        return None

    candidate_paths = sorted(glob.glob(os.path.join(candidate_dir, "000*.vcf*")))
    if not candidate_paths:
        return set()

    keys = set()
    for path in candidate_paths:
        _, _, _, records = read_vcf("candidate", path)
        for record in records:
            keys.add(record.key)
    return keys


def support_records(by_key, min_callers, caller_order):
    selected = []
    for key, callers in by_key.items():
        support = [caller for caller in caller_order if caller in callers]
        if len(support) >= min_callers:
            chosen_caller = next(caller for caller in caller_order if caller in callers)
            selected.append((key, chosen_caller, support, callers))
    return selected


def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == len(alt):
        return "mnv"
    if len(ref) < len(alt):
        return "insertion" if alt.startswith(ref) else "complex"
    if len(ref) > len(alt):
        return "deletion" if ref.startswith(alt) else "complex"
    return "complex"


def parse_number(value):
    if value in (None, "", "."):
        return None
    value = str(value).split(",")[0]
    try:
        return float(value)
    except ValueError:
        return None


def fmt_number(value):
    if value is None:
        return "."
    return f"{value:.6g}"


def mean(values):
    return sum(values) / len(values) if values else None


def median(values):
    if not values:
        return None
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def pct(count, denominator):
    if denominator == 0:
        return "."
    return f"{100 * count / denominator:.3f}"


def caller_key_sets(by_key, caller_order):
    key_sets = {caller: set() for caller in caller_order}
    for key, callers in by_key.items():
        for caller in caller_order:
            if caller in callers:
                key_sets[caller].add(key)
    return key_sets


def caller_combination_counts(by_key, caller_order):
    counts = Counter()
    support_counts = Counter()
    for callers in by_key.values():
        support = tuple(caller for caller in caller_order if caller in callers)
        counts[support] += 1
        support_counts[len(support)] += 1
    return counts, support_counts


def write_caller_characteristics(path, sample, by_key, caller_order):
    header = [
        "sample", "caller", "total_variants", "snv", "mnv", "insertion",
        "deletion", "complex", "mean_af", "median_af", "mean_dp",
        "median_dp", "mean_qual", "median_qual", "singleton_variants",
        "singleton_pct", "concordant_2plus_variants", "concordant_2plus_pct",
        "all_three_variants", "all_three_pct", "filter_counts",
    ]
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for caller in caller_order:
            records = [callers[caller] for callers in by_key.values() if caller in callers]
            total = len(records)
            type_counts = Counter(variant_type(record.ref, record.alt) for record in records)
            filter_counts = Counter()
            af_values = []
            dp_values = []
            qual_values = []
            singleton = 0
            concordant_2plus = 0
            all_three = 0

            for callers in by_key.values():
                if caller not in callers:
                    continue
                support_size = len(callers)
                if support_size == 1:
                    singleton += 1
                if support_size >= 2:
                    concordant_2plus += 1
                if support_size == len(caller_order):
                    all_three += 1

            for record in records:
                evidence = record.evidence()
                filter_counts[evidence["FILTER"] or "."] += 1
                af = parse_number(evidence["AF"])
                dp = parse_number(evidence["DP"])
                qual = parse_number(record.qual)
                if af is not None:
                    af_values.append(af)
                if dp is not None:
                    dp_values.append(dp)
                if qual is not None:
                    qual_values.append(qual)

            filter_summary = ",".join(f"{key}={value}" for key, value in sorted(filter_counts.items())) or "."
            row = [
                sample, caller, str(total),
                str(type_counts["snv"]), str(type_counts["mnv"]),
                str(type_counts["insertion"]), str(type_counts["deletion"]),
                str(type_counts["complex"]),
                fmt_number(mean(af_values)), fmt_number(median(af_values)),
                fmt_number(mean(dp_values)), fmt_number(median(dp_values)),
                fmt_number(mean(qual_values)), fmt_number(median(qual_values)),
                str(singleton), pct(singleton, total),
                str(concordant_2plus), pct(concordant_2plus, total),
                str(all_three), pct(all_three, total),
                filter_summary,
            ]
            out.write("\t".join(row) + "\n")


def write_concordance_summary(path, sample, by_key, caller_order):
    combo_counts, support_counts = caller_combination_counts(by_key, caller_order)
    union_total = len(by_key)
    header = ["sample", "category", "callers", "support_count", "variant_count", "pct_union"]
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for support_size in range(1, len(caller_order) + 1):
            count = support_counts[support_size]
            out.write("\t".join([
                sample, "support_size", f"support_{support_size}",
                str(support_size), str(count), pct(count, union_total),
            ]) + "\n")

        for support_size in range(1, len(caller_order) + 1):
            for combo in combinations(caller_order, support_size):
                count = combo_counts[combo]
                out.write("\t".join([
                    sample, "caller_combination", "+".join(combo),
                    str(support_size), str(count), pct(count, union_total),
                ]) + "\n")


def write_pairwise_concordance(path, sample, by_key, caller_order):
    key_sets = caller_key_sets(by_key, caller_order)
    header = [
        "sample", "caller_a", "caller_b", "caller_a_total",
        "caller_b_total", "shared", "caller_a_only", "caller_b_only",
        "union", "jaccard", "shared_pct_of_a", "shared_pct_of_b",
    ]
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for caller_a, caller_b in combinations(caller_order, 2):
            a_keys = key_sets[caller_a]
            b_keys = key_sets[caller_b]
            shared = len(a_keys & b_keys)
            union = len(a_keys | b_keys)
            row = [
                sample, caller_a, caller_b, str(len(a_keys)), str(len(b_keys)),
                str(shared), str(len(a_keys - b_keys)), str(len(b_keys - a_keys)),
                str(union), fmt_number(shared / union if union else None),
                pct(shared, len(a_keys)), pct(shared, len(b_keys)),
            ]
            out.write("\t".join(row) + "\n")


def sort_key(item, contig_order):
    key = item[0]
    chrom, pos, ref, alt = key
    return contig_order.get(chrom, len(contig_order) + 1), chrom, pos, ref, alt


def write_outputs(output_vcf, support_tsv, selected, header_meta, chrom_header, contig_order):
    selected = sorted(selected, key=lambda item: sort_key(item, contig_order))
    with open(output_vcf, "w", encoding="utf-8") as vcf, open(support_tsv, "w", encoding="utf-8") as tsv:
        for line in header_meta:
            vcf.write(line + "\n")
        vcf.write(chrom_header + "\n")

        tsv.write("\t".join([
            "sample", "chrom", "pos", "ref", "alt", "support", "callers",
            "selected_caller", "selected_qual", "selected_filter",
            "mutect2_af", "mutect2_dp", "mutect2_ad",
            "deepsomatic_af", "deepsomatic_dp", "deepsomatic_ad",
            "clairsto_af", "clairsto_dp", "clairsto_ad",
        ]) + "\n")

        sample = chrom_header.split("\t")[-1] if chrom_header.count("\t") >= 9 else ""
        for key, chosen_caller, support, callers in selected:
            chosen = callers[chosen_caller]
            fields = chosen.fields[:]
            if len(fields) >= 10 and sample:
                fields[9] = chosen.fields[9]

            info = chosen.info
            for custom_key in ["SUPPORT", "CALLERS", "SELECTED_CALLER"]:
                info.pop(custom_key, None)
            for prefix in CALLER_PREFIX.values():
                for suffix in EVIDENCE_SUFFIXES:
                    info.pop(f"{prefix}_{suffix}", None)

            info["SUPPORT"] = str(len(support))
            info["CALLERS"] = ",".join(support)
            info["SELECTED_CALLER"] = chosen_caller

            evidence_by_caller = {}
            for caller in support:
                ev = callers[caller].evidence()
                evidence_by_caller[caller] = ev
                prefix = CALLER_PREFIX[caller]
                for suffix in EVIDENCE_SUFFIXES:
                    info[f"{prefix}_{suffix}"] = ev[suffix]

            fields[7] = format_info(info)
            vcf.write("\t".join(fields) + "\n")

            def ev(caller, field):
                return evidence_by_caller.get(caller, {}).get(field, ".")

            tsv.write("\t".join([
                sample,
                key[0], str(key[1]), key[2], key[3],
                str(len(support)), ",".join(support),
                chosen_caller, chosen.qual, chosen.filt,
                ev("mutect2", "AF"), ev("mutect2", "DP"), ev("mutect2", "AD"),
                ev("deepsomatic", "AF"), ev("deepsomatic", "DP"), ev("deepsomatic", "AD"),
                ev("clairsto", "AF"), ev("clairsto", "DP"), ev("clairsto", "AD"),
            ]) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Merge tumor-only caller VCFs with support/evidence INFO fields.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--min-callers", type=int, default=2)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--support-tsv", required=True)
    parser.add_argument("--all-output-vcf")
    parser.add_argument("--all-support-tsv")
    parser.add_argument("--caller-characteristics-tsv")
    parser.add_argument("--concordance-summary-tsv")
    parser.add_argument("--pairwise-concordance-tsv")
    parser.add_argument("--candidate-dir", help="Optional bcftools isec output directory used to gate concordant calls.")
    # Defaults to OFF at the script level on purpose: a run already in flight generated
    # its CONCORDANCE_MERGE command line from an older main.nf that does not pass this
    # flag, and must keep behaving exactly as before. main.nf sets the pipeline-level
    # default to true, so new runs are gated. See somatic_index() for why.
    parser.add_argument("--require-pass", action="store_true",
                        help="Only let a caller vote when its own FILTER is PASS. "
                             "Without this, germline sites dominate the concordance "
                             "counts and support_2 can exceed support_1.")
    parser.add_argument("--caller", action="append", required=True, help="caller_id:path.vcf.gz")
    args = parser.parse_args()

    caller_paths = OrderedDict()
    for value in args.caller:
        if ":" not in value:
            raise ValueError(f"--caller must be caller_id:path, got {value}")
        caller, path = value.split(":", 1)
        if caller not in CALLER_PREFIX:
            raise ValueError(f"Unsupported caller_id {caller}; expected one of {','.join(CALLER_PREFIX)}")
        caller_paths[caller] = path

    caller_order = [caller for caller in ("mutect2", "deepsomatic", "clairsto") if caller in caller_paths]
    vcf_parts = []
    all_records = []
    for caller in caller_order:
        part = read_vcf(caller, caller_paths[caller])
        vcf_parts.append(part)
        all_records.extend(part[3])

    header_meta, chrom_header, contig_order = merge_headers(vcf_parts, args.sample)
    by_key = build_index(all_records)

    # Everything downstream -- the vote, the union file, the characteristics table, the
    # concordance summary and the pairwise table -- must run off the same index, or the
    # tables stop describing the VCFs beside them. That mismatch is precisely what made
    # the published summaries unusable.
    if args.require_pass:
        by_key = somatic_index(by_key)
        header_meta = list(header_meta) + [
            "##somatic_gate=caller_native_FILTER==PASS",
        ]

    concordant = support_records(by_key, args.min_callers, caller_order)
    candidate_keys = load_candidate_keys(args.candidate_dir)
    if candidate_keys is not None:
        concordant = [item for item in concordant if item[0] in candidate_keys]
    write_outputs(args.output_vcf, args.support_tsv, concordant, header_meta, chrom_header, contig_order)

    if args.all_output_vcf and args.all_support_tsv:
        union = support_records(by_key, 1, caller_order)
        write_outputs(args.all_output_vcf, args.all_support_tsv, union, header_meta, chrom_header, contig_order)

    if args.caller_characteristics_tsv:
        write_caller_characteristics(args.caller_characteristics_tsv, args.sample, by_key, caller_order)
    if args.concordance_summary_tsv:
        write_concordance_summary(args.concordance_summary_tsv, args.sample, by_key, caller_order)
    if args.pairwise_concordance_tsv:
        write_pairwise_concordance(args.pairwise_concordance_tsv, args.sample, by_key, caller_order)


if __name__ == "__main__":
    main()
