#!/usr/bin/env python3
import argparse
import gzip
from collections import Counter


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def parse_info(info_text):
    info = {}
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


def first_value(container, keys):
    for key in keys:
        if key in container and container[key] not in ("", "."):
            return container[key]
    return ""


def parse_sample_fields(fields):
    if len(fields) < 10:
        return {}, {}
    fmt_keys = fields[8].split(":")
    fmt_vals = fields[9].split(":")
    return dict(zip(fmt_keys, fmt_vals)), parse_info(fields[7] if len(fields) > 7 else ".")


def parse_float(value):
    if value in (None, "", "."):
        return None
    try:
        return float(str(value).split(",")[0])
    except ValueError:
        return None


def parse_int(value):
    if value in (None, "", "."):
        return None
    try:
        return int(round(float(str(value).split(",")[0])))
    except ValueError:
        return None


def parse_number_list(value):
    try:
        return [float(x) for x in str(value).split(",") if x not in ("", ".")]
    except ValueError:
        return []


def calc_from_ad(ad_value):
    parts = parse_number_list(ad_value)
    if len(parts) < 2:
        return None, None, None
    total = sum(parts)
    if total <= 0:
        return int(total), None, None
    alt_reads = int(round(max(parts[1:])))
    af = max(parts[1:]) / total
    return int(round(total)), af, alt_reads


def extract_metrics(fields):
    fmt, info = parse_sample_fields(fields)
    af = parse_float(first_value(fmt, ("AF", "VAF", "T_AF", "TUMOR_AF")))
    if af is None:
        af = parse_float(first_value(info, ("AF", "VAF", "TLOD_AF", "ALLELE_FRACTION")))

    dp = parse_int(first_value(fmt, ("DP", "T_DP")))
    if dp is None:
        dp = parse_int(first_value(info, ("DP", "TDP")))

    ad_value = first_value(fmt, ("AD", "T_AD"))
    alt_reads = None
    if ad_value:
        ad_dp, ad_af, ad_alt = calc_from_ad(ad_value)
        if dp is None:
            dp = ad_dp
        if af is None:
            af = ad_af
        alt_reads = ad_alt

    if alt_reads is None:
        alt_reads = parse_int(first_value(fmt, ("AO", "NV", "ALT_COUNT")))
    if alt_reads is None:
        alt_reads = parse_int(first_value(info, ("AO", "NV", "ALT_COUNT")))

    return af, dp, alt_reads


def variant_bucket(ref, alt):
    return "snv" if len(ref) == 1 and len(alt) == 1 else "indel"


def thresholds_for(bucket, args):
    if bucket == "snv":
        return args.snv_min_af, args.snv_min_alt_reads, args.snv_min_dp
    return args.indel_min_af, args.indel_min_alt_reads, args.indel_min_dp


def write_summary(path, sample, caller, total_records, passed_records, reasons):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join([
            "sample",
            "caller",
            "total_records",
            "passed_records",
            "failed_records",
            "missing_metrics",
            "failed_af",
            "failed_alt_reads",
            "failed_dp",
        ]) + "\n")
        handle.write("\t".join([
            sample,
            caller,
            str(total_records),
            str(passed_records),
            str(total_records - passed_records),
            str(reasons["missing_metrics"]),
            str(reasons["failed_af"]),
            str(reasons["failed_alt_reads"]),
            str(reasons["failed_dp"]),
        ]) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Apply a shared tumor-only post-filter to a normalized VCF.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--caller", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--snv-min-af", required=True, type=float)
    parser.add_argument("--indel-min-af", required=True, type=float)
    parser.add_argument("--snv-min-alt-reads", required=True, type=int)
    parser.add_argument("--indel-min-alt-reads", required=True, type=int)
    parser.add_argument("--snv-min-dp", required=True, type=int)
    parser.add_argument("--indel-min-dp", required=True, type=int)
    args = parser.parse_args()

    total_records = 0
    passed_records = 0
    reasons = Counter()

    with open_text(args.vcf) as in_handle, open(args.out, "w", encoding="utf-8") as out_handle:
        for raw in in_handle:
            if raw.startswith("#"):
                out_handle.write(raw)
                continue

            total_records += 1
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 8:
                reasons["missing_metrics"] += 1
                continue

            bucket = variant_bucket(fields[3], fields[4])
            min_af, min_alt_reads, min_dp = thresholds_for(bucket, args)
            af, dp, alt_reads = extract_metrics(fields)

            if af is None or dp is None or alt_reads is None:
                reasons["missing_metrics"] += 1
                continue

            failed = False
            if af < min_af:
                reasons["failed_af"] += 1
                failed = True
            if alt_reads < min_alt_reads:
                reasons["failed_alt_reads"] += 1
                failed = True
            if dp < min_dp:
                reasons["failed_dp"] += 1
                failed = True

            if failed:
                continue

            passed_records += 1
            out_handle.write(raw)

    write_summary(args.summary, args.sample, args.caller, total_records, passed_records, reasons)


if __name__ == "__main__":
    main()
