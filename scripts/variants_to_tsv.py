#!/usr/bin/env python3
import argparse
import gzip
from collections import OrderedDict


IMPACT_RANK = {
    "HIGH": 0,
    "MODERATE": 1,
    "LOW": 2,
    "MODIFIER": 3,
}


EXTRA_INFO_FIELDS = [
    "CLNSIG", "CLNDN", "CLNREVSTAT", "CLNVC", "CLNHGVS",
    "COSMIC", "COSMIC_ID", "HOTSPOT", "CIVIC_ID", "ONCOKB",
    "AF", "AF_POPMAX", "gnomAD_AF", "gnomAD_AF_POPMAX",
]


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
            info[item] = "TRUE"
    return info


def choose_ann(info):
    anns = info.get("ANN", "")
    if not anns:
        return {}
    best = None
    best_rank = 99
    for ann in anns.split(","):
        fields = ann.split("|")
        impact = fields[2] if len(fields) > 2 else ""
        rank = IMPACT_RANK.get(impact, 98)
        if best is None or rank < best_rank:
            best = fields
            best_rank = rank
    if best is None:
        return {}
    return {
        "effect": best[1] if len(best) > 1 else "",
        "impact": best[2] if len(best) > 2 else "",
        "gene": best[3] if len(best) > 3 else "",
        "gene_id": best[4] if len(best) > 4 else "",
        "feature": best[6] if len(best) > 6 else "",
        "biotype": best[7] if len(best) > 7 else "",
        "hgvs_c": best[9] if len(best) > 9 else "",
        "hgvs_p": best[10] if len(best) > 10 else "",
    }


def main():
    parser = argparse.ArgumentParser(description="Create a tumor-only variant review TSV from an annotated VCF.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    columns = [
        "sample", "chrom", "pos", "ref", "alt", "qual", "filter",
        "support", "callers", "selected_caller",
        "gene", "effect", "impact", "hgvs_c", "hgvs_p", "feature", "biotype",
        "mutect2_af", "mutect2_dp", "mutect2_ad",
        "deepsomatic_af", "deepsomatic_dp", "deepsomatic_ad",
        "clairsto_af", "clairsto_dp", "clairsto_ad",
    ] + EXTRA_INFO_FIELDS

    with open_text(args.vcf) as handle, open(args.out, "w", encoding="utf-8") as out:
        out.write("\t".join(columns) + "\n")
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            ann = choose_ann(info)

            row = {
                "sample": args.sample,
                "chrom": fields[0],
                "pos": fields[1],
                "ref": fields[3],
                "alt": fields[4],
                "qual": fields[5],
                "filter": fields[6],
                "support": info.get("SUPPORT", "."),
                "callers": info.get("CALLERS", "."),
                "selected_caller": info.get("SELECTED_CALLER", "."),
                "gene": ann.get("gene", "."),
                "effect": ann.get("effect", "."),
                "impact": ann.get("impact", "."),
                "hgvs_c": ann.get("hgvs_c", "."),
                "hgvs_p": ann.get("hgvs_p", "."),
                "feature": ann.get("feature", "."),
                "biotype": ann.get("biotype", "."),
                "mutect2_af": info.get("M2_AF", "."),
                "mutect2_dp": info.get("M2_DP", "."),
                "mutect2_ad": info.get("M2_AD", "."),
                "deepsomatic_af": info.get("DS_AF", "."),
                "deepsomatic_dp": info.get("DS_DP", "."),
                "deepsomatic_ad": info.get("DS_AD", "."),
                "clairsto_af": info.get("CT_AF", "."),
                "clairsto_dp": info.get("CT_DP", "."),
                "clairsto_ad": info.get("CT_AD", "."),
            }
            for key in EXTRA_INFO_FIELDS:
                row[key] = info.get(key, ".")

            out.write("\t".join(str(row.get(col, ".")) for col in columns) + "\n")


if __name__ == "__main__":
    main()
