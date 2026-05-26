#!/usr/bin/env python3

import csv
import glob
import math
import os
from collections import defaultdict

import matplotlib.pyplot as plt

BASE = "/media/daniel/Espacio/RefinedIBD"
FILE_LC = os.path.join(BASE, "LC/chr*_ibd_merged.ibd")
FILE_HC = os.path.join(BASE, "HC/chr*_ibd_merged.ibd")

OUT_GLOBAL = os.path.join(BASE, "IBD_metrics_global.csv")
OUT_PER_PAIR = os.path.join(BASE, "IBD_metrics_per_pair.csv")
OUT_BY_LENGTH = os.path.join(BASE, "IBD_metrics_by_length.csv")
OUT_PDF = os.path.join(BASE, "IBD_metrics_recall_accuracy_f1_by_length.pdf")

BIN_WIDTH_CM = 0.5

CHR_SIZES = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
}


def clean_id(value):
    return value.replace("_1X", "").replace("_30X", "")


def pair_key(chrom, id1, id2):
    a, b = sorted([clean_id(id1), clean_id(id2)])
    return chrom, a, b


def length_bin(length_cm):
    start = math.floor(length_cm / BIN_WIDTH_CM) * BIN_WIDTH_CM
    end = start + BIN_WIDTH_CM
    mid = start + BIN_WIDTH_CM / 2
    return round(start, 6), round(end, 6), round(mid, 6)


def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [list(intervals[0])]

    for start, end in intervals[1:]:
        last = merged[-1]
        if start <= last[1]:
            last[1] = max(last[1], end)
        else:
            merged.append([start, end])

    return [(start, end) for start, end in merged]


def interval_sum(intervals):
    return sum(end - start for start, end in intervals)


def intersection_length(a_intervals, b_intervals):
    total = 0
    i = 0
    j = 0

    while i < len(a_intervals) and j < len(b_intervals):
        a_start, a_end = a_intervals[i]
        b_start, b_end = b_intervals[j]

        ov_start = max(a_start, b_start)
        ov_end = min(a_end, b_end)
        if ov_end > ov_start:
            total += ov_end - ov_start

        if a_end < b_end:
            i += 1
        else:
            j += 1

    return total


def union_length(a_intervals, b_intervals):
    return interval_sum(merge_intervals(list(a_intervals) + list(b_intervals)))


def safe_div(num, den):
    return num / den if den else ""


def load_ibd(pattern, label):
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found for pattern: {pattern}")

    rows = []
    total_rows = 0
    skipped_self = 0
    skipped_bad = 0

    print(f"\nLoading {label}: {len(files)} chromosome files")

    for path in files:
        with open(path, "rt") as handle:
            for line in handle:
                total_rows += 1
                fields = line.strip().split()

                if len(fields) < 9:
                    skipped_bad += 1
                    continue

                id1 = fields[0]
                id2 = fields[2]
                chrom = fields[4]

                if clean_id(id1) == clean_id(id2):
                    skipped_self += 1
                    continue

                try:
                    start = int(fields[5])
                    end = int(fields[6])
                    length_cm = float(fields[8])
                except ValueError:
                    skipped_bad += 1
                    continue

                if chrom not in CHR_SIZES or end <= start:
                    skipped_bad += 1
                    continue

                bin_start, bin_end, bin_mid = length_bin(length_cm)
                rows.append({
                    "key": pair_key(chrom, id1, id2),
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "length_bp": end - start,
                    "length_cm": length_cm,
                    "bin_start_cm": bin_start,
                    "bin_end_cm": bin_end,
                    "bin_mid_cm": bin_mid,
                })

    print(f"  rows read: {total_rows}")
    print(f"  rows kept: {len(rows)}")
    print(f"  same-individual skipped: {skipped_self}")
    print(f"  malformed skipped: {skipped_bad}")

    return rows


def build_index(rows):
    raw = defaultdict(list)
    for row in rows:
        raw[row["key"]].append((row["start"], row["end"]))
    return {key: merge_intervals(intervals) for key, intervals in raw.items()}


def build_bin_index(rows):
    raw = defaultdict(list)
    for row in rows:
        bin_key = (row["bin_start_cm"], row["bin_end_cm"], row["bin_mid_cm"])
        raw[(bin_key, row["key"])].append((row["start"], row["end"]))
    return {key: merge_intervals(intervals) for key, intervals in raw.items()}


def metrics_from_counts(tp, fn, fp, tn):
    # In this analysis, "accuracy" is the fraction of LC-positive sequence that overlaps HC.
    accuracy = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    f1 = safe_div(2 * tp, 2 * tp + fp + fn)
    fnr = safe_div(fn, fn + tp)
    fpr = safe_div(fp, fp + tn)
    tnr = safe_div(tn, tn + fp)

    return {
        "accuracy": accuracy,
        "recall": recall,
        "f1": f1,
        "FNR": fnr,
        "FPR": fpr,
        "TNR": tnr,
    }


def compute_global_metrics(hc_index, lc_index):
    rows = []
    totals = {"TP_bp": 0, "FN_bp": 0, "FP_bp": 0, "TN_bp": 0}

    for key in sorted(set(hc_index) | set(lc_index)):
        chrom, ind1, ind2 = key
        chrom_size = CHR_SIZES[chrom]
        hc_intervals = hc_index.get(key, [])
        lc_intervals = lc_index.get(key, [])

        tp = intersection_length(hc_intervals, lc_intervals)
        hc_positive = interval_sum(hc_intervals)
        lc_positive = interval_sum(lc_intervals)
        fn = hc_positive - tp
        fp = lc_positive - tp
        tn = max(0, chrom_size - union_length(hc_intervals, lc_intervals))

        row = {
            "chr": chrom,
            "ind1": ind1,
            "ind2": ind2,
            "TP_bp": tp,
            "FN_bp": fn,
            "FP_bp": fp,
            "TN_bp": tn,
            "HC_positive_bp": hc_positive,
            "LC_positive_bp": lc_positive,
            "chr_size_bp": chrom_size,
            **metrics_from_counts(tp, fn, fp, tn),
        }
        rows.append(row)

        for field in totals:
            totals[field] += row[field]

    global_row = {
        **totals,
        **metrics_from_counts(totals["TP_bp"], totals["FN_bp"], totals["FP_bp"], totals["TN_bp"]),
    }
    return rows, global_row


def compute_by_length_metrics(hc_bin_index, lc_bin_index):
    by_bin = defaultdict(lambda: {"TP_bp": 0, "FN_bp": 0, "FP_bp": 0, "TN_bp": 0})
    all_bin_keys = sorted(set(key[0] for key in hc_bin_index) | set(key[0] for key in lc_bin_index))

    for bin_key in all_bin_keys:
        keys_for_bin = {
            key for current_bin, key in hc_bin_index if current_bin == bin_key
        } | {
            key for current_bin, key in lc_bin_index if current_bin == bin_key
        }

        for key in keys_for_bin:
            chrom = key[0]
            chrom_size = CHR_SIZES[chrom]
            hc_intervals = hc_bin_index.get((bin_key, key), [])
            lc_intervals = lc_bin_index.get((bin_key, key), [])

            tp = intersection_length(hc_intervals, lc_intervals)
            fn = interval_sum(hc_intervals) - tp
            fp = interval_sum(lc_intervals) - tp
            tn = max(0, chrom_size - union_length(hc_intervals, lc_intervals))

            by_bin[bin_key]["TP_bp"] += tp
            by_bin[bin_key]["FN_bp"] += fn
            by_bin[bin_key]["FP_bp"] += fp
            by_bin[bin_key]["TN_bp"] += tn

    rows = []
    for bin_key, counts in sorted(by_bin.items()):
        bin_start, bin_end, bin_mid = bin_key
        rows.append({
            "bin_start_cm": bin_start,
            "bin_end_cm": bin_end,
            "length_cm": bin_mid,
            **counts,
            **metrics_from_counts(counts["TP_bp"], counts["FN_bp"], counts["FP_bp"], counts["TN_bp"]),
        })

    return rows


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def numeric(value):
    return value if isinstance(value, (int, float)) else float("nan")


def plot_by_length(rows):
    rows = [row for row in rows if row["TP_bp"] + row["FN_bp"] + row["FP_bp"] > 0]
    x = [row["length_cm"] for row in rows]
    series = {
        "Recall": [numeric(row["recall"]) for row in rows],
        "Accuracy": [numeric(row["accuracy"]) for row in rows],
        "F1": [numeric(row["f1"]) for row in rows],
    }
    markers = {"Recall": "o", "Accuracy": "x", "F1": "s"}
    finite_values = [
        value for values in series.values() for value in values
        if isinstance(value, (int, float)) and not math.isnan(value)
    ]
    zoom_top = min(1.0, max(finite_values) * 1.15) if finite_values else 1.0
    zoom_top = max(zoom_top, 0.05)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=True)
    for ax in axes:
        for label, values in series.items():
            ax.scatter(x, values, label=label, alpha=0.8, marker=markers[label])
        ax.set_xlabel("IBD segment length (cM)")
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("Metric")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title("Full scale")
    axes[1].set_ylim(0, zoom_top)
    axes[1].set_title("Zoomed scale")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3)
    fig.suptitle("LC vs HC IBD metrics by segment length")
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    fig.savefig(OUT_PDF)


def main():
    hc_rows = load_ibd(FILE_HC, "HC reference")
    lc_rows = load_ibd(FILE_LC, "LC prediction")

    hc_index = build_index(hc_rows)
    lc_index = build_index(lc_rows)
    hc_bin_index = build_bin_index(hc_rows)
    lc_bin_index = build_bin_index(lc_rows)

    per_pair_rows, global_row = compute_global_metrics(hc_index, lc_index)
    by_length_rows = compute_by_length_metrics(hc_bin_index, lc_bin_index)

    per_pair_fields = [
        "chr", "ind1", "ind2", "TP_bp", "FN_bp", "FP_bp", "TN_bp",
        "HC_positive_bp", "LC_positive_bp", "chr_size_bp",
        "accuracy", "recall", "f1", "FNR", "FPR", "TNR",
    ]
    by_length_fields = [
        "bin_start_cm", "bin_end_cm", "length_cm", "TP_bp", "FN_bp", "FP_bp", "TN_bp",
        "accuracy", "recall", "f1", "FNR", "FPR", "TNR",
    ]

    write_csv(OUT_PER_PAIR, per_pair_rows, per_pair_fields)
    write_csv(OUT_GLOBAL, [global_row], list(global_row.keys()))
    write_csv(OUT_BY_LENGTH, by_length_rows, by_length_fields)
    plot_by_length(by_length_rows)

    print("\nGlobal metrics, HC as reference and LC as prediction")
    for key in ["TP_bp", "FN_bp", "FP_bp", "TN_bp"]:
        print(f"  {key}: {global_row[key]}")
    for key in ["recall", "accuracy", "f1"]:
        print(f"  {key}: {global_row[key]:.6f}")

    print("\nDONE")
    print(f"  Global CSV: {OUT_GLOBAL}")
    print(f"  Per-pair CSV: {OUT_PER_PAIR}")
    print(f"  By-length CSV: {OUT_BY_LENGTH}")
    print(f"  Plot: {OUT_PDF}")


if __name__ == "__main__":
    main()
