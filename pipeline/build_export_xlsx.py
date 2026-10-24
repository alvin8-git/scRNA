#!/usr/bin/env python3
"""Pack a run's export CSVs into one workbook for the customer.

08f_customer_export.R writes the derived tables as CSVs; this collects them
into <run>/export/<run>_analysis.xlsx. Kept in Python because this conda env
has no R xlsx writer (openxlsx, writexl and xlsx are all absent) but does have
xlsxwriter.

Sheet order is deliberate: provenance first, then what the samples are, then
what was found, then the evidence, then raw QC. A customer opening this should
meet the run's parameters before any numbers.

Usage: python3 pipeline/build_export_xlsx.py <run_dir>
"""
import csv, os, sys, glob

ORDER = [
    ("run_parameters",    "Run_parameters",    "How this run was configured and which software produced it"),
    ("sample_overview",   "Sample_overview",   "Cells and quality metrics per sample"),
    ("qc_cell_fate",      "QC_cell_fate",      "How many cells each filter removed, per sample"),
    ("qc_summary",        "QC_summary",        "Per-sample QC summary as the pipeline recorded it"),
    ("ambient_rho",       "Ambient_rho",       "SoupX ambient RNA fraction per sample"),
    ("celltype_counts",   "CellType_counts",   "Cells and percent per cell type, per sample"),
    ("cluster_evidence",  "Cluster_evidence",  "Why each cluster received its label"),
    ("annotation_detail", "Annotation_detail", "Cluster table as written by the annotation step"),
    ("marker_panel",      "Marker_panel",      "The marker panel used, and which genes were present"),
    ("cluster_markers",   "Cluster_markers",   "Top differentially expressed genes per cluster"),
]

def main(run_dir):
    exp = os.path.join(run_dir, "export")
    if not os.path.isdir(exp):
        sys.exit(f"no export directory at {exp} — run 08f_customer_export.R first")
    import xlsxwriter
    run = os.path.basename(run_dir.rstrip("/"))
    out = os.path.join(exp, f"{run}_analysis.xlsx")
    wb = xlsxwriter.Workbook(out, {"constant_memory": True})
    hdr  = wb.add_format({"bold": True, "bg_color": "#1a4f7a", "font_color": "white",
                          "border": 1, "text_wrap": True, "valign": "top"})
    bold = wb.add_format({"bold": True})
    wrap = wb.add_format({"text_wrap": True, "valign": "top"})
    num  = wb.add_format({"num_format": "0.00"})
    intf = wb.add_format({"num_format": "#,##0"})

    contents = []
    for stem, sheet, desc in ORDER:
        path = os.path.join(exp, stem + ".csv")
        if not os.path.exists(path):
            continue
        ws = wb.add_worksheet(sheet)
        with open(path, newline="", encoding="utf8") as fh:
            rows = csv.reader(fh)
            head = next(rows)
            for c, h in enumerate(head):
                ws.write(0, c, h, hdr)
            # width from the header, widened for free-text marker columns
            for c, h in enumerate(head):
                w = 14 if len(h) < 14 else min(len(h) + 2, 30)
                if "marker" in h.lower() or "enriched" in h.lower() or h in ("Value", "Parameter"):
                    w = 52
                ws.set_column(c, c, w, wrap if w >= 52 else None)
            ws.freeze_panes(1, 0)
            ws.autofilter(0, 0, 0, len(head) - 1)
            n = 0
            for i, row in enumerate(rows, start=1):
                for c, v in enumerate(row):
                    try:
                        f = float(v)
                        ws.write_number(i, c, f, intf if f == int(f) and abs(f) >= 1000 else num)
                    except ValueError:
                        ws.write_string(i, c, v)
                n = i
        contents.append((sheet, desc, n))
        print(f"  {sheet:<20} {n:>6} rows")

    # contents sheet, written last but moved to the front by index
    ws = wb.add_worksheet("Contents")
    ws.set_column(0, 0, 24); ws.set_column(1, 1, 64); ws.set_column(2, 2, 10)
    ws.write(0, 0, f"{run} — analysis data", bold)
    ws.write(2, 0, "Sheet", hdr); ws.write(2, 1, "Contents", hdr); ws.write(2, 2, "Rows", hdr)
    for i, (sheet, desc, n) in enumerate(contents, start=3):
        ws.write(i, 0, sheet); ws.write(i, 1, desc, wrap); ws.write_number(i, 2, n)
    r = len(contents) + 4
    ws.write(r, 0, "Also in this folder", bold)
    for f, d in [("per_cell_assignments.csv.gz",
                  "One row per cell: sample, cluster, cell type, UMIs, genes, mito%, UMAP position. "
                  "Every figure in the report can be reproduced from this file."),
                 ("pseudobulk_cpm.csv.gz",
                  "Expression per gene per sample, counts-per-million normalised so samples compare directly.")]:
        p = os.path.join(exp, f)
        if os.path.exists(p):
            r += 1
            ws.write(r, 0, f); ws.write(r, 1, d, wrap)
            ws.write_number(r, 2, round(os.path.getsize(p) / 1e6, 2))
    r += 2
    ws.write(r, 0, "Note", bold)
    ws.write(r, 1, "Cell-type labels are a computational call, not a measurement. Cluster_evidence "
                   "shows the marker support for each one; a cluster labelled Unknown is one where "
                   "the evidence did not support a confident call and is reported rather than guessed.",
             wrap)
    wb.close()
    # Contents was added last; xlsxwriter has no reorder, so state the intent
    print(f"\nwrote {out} ({os.path.getsize(out)/1e6:.2f} MB, {len(contents)+1} sheets)")
    print("XLSX_EXPORT_DONE")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: build_export_xlsx.py <run_dir>")
    main(sys.argv[1])
