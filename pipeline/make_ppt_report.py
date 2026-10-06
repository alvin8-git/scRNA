#!/usr/bin/env python3
"""Generate a PowerPoint report for a finished scRNA run, on the MGI template.

Slide order follows the flowcell metrics spreadsheet, NOT alphabetical order:
  1  title
  2  sequencing metrics table (read from the xlsx at runtime)
  3+ cell-type UMAPs, 2 samples per slide, in spreadsheet column order
     contamination
     cell-type proportions
     canonical markers
  last  "Thank you" (kept from the template)

Nothing is sample-specific: point it at another run + that flowcell's xlsx and it
produces the matching deck.

  python3 pipeline/make_ppt_report.py \
      --run-dir  Results/results_PBMC5KAQL_6samples_filtered \
      --xls      Mirxes/Mirxes_PBMC_metrics.xlsx \
      --template "Mirxes/MGI PPT templates_2026_v1.pptx"
"""
import argparse, datetime, re, subprocess, sys
from pathlib import Path

from openpyxl import load_workbook
from pptx import Presentation
from pptx.util import Inches, Pt

TMP = Path("/data/alvin/tmp/ppt_figs")          # project rule: never /tmp
LAYOUT_TITLE_ONLY = 3                            # 仅标题 — what the reference deck uses
DPI = 200

# figure -> (relative path, slide title(s)). A list of titles names each page of a multi-page PDF.
# celltype_proportions_bar.pdf is deliberately absent: it is the "Labels: de-novo" variant of the
# same % chart, and it is not ordered by the sample list.
FIGURES = [
    ("annotation/contamination_summary.pdf",         "Contamination"),
    ("integrated/celltype_composition_combined.pdf",
     ["Cell Type – % of Sample", "Cell Type – Number of Cells"]),
    ("annotation/canonical_markers_dotplot.pdf",      "Cell Markers"),
]
UMAP_PDF = "integrated/umap_split_by_sample.pdf"

# Metrics rows to drop so the table fits one slide. Everything else in the xlsx is kept, so a
# flowcell with extra rows still reports them.
DROP_EXACT = {
    "reads mapped confidently to genome",
    "reads mapped confidently to transcriptome",
    "include introns",
    "sample name",            # repeats the header row
}
# of the cDNA_/oligo_ block keep only these two measures
KEEP_SEQ_SUFFIX = ("number of reads", "valid barcodes")


def keep_metric(label: str, seen_before: bool) -> bool:
    low = label.strip().lower()
    if low in DROP_EXACT:
        return False
    if low.startswith(("cdna_", "oligo_")):
        return low.split("_", 1)[1].strip() in KEEP_SEQ_SUFFIX
    # the sheet repeats "Estimated number of cells" and "Median UMI counts per cell" at the end
    # (the second one spelled "Estimated number ofcells") — keep only the first occurrence
    if seen_before:
        return False
    return True


def sample_order_from_xls(xls: Path):
    """Sample names in spreadsheet column order.

    The metrics sheet is transposed: row 1 holds the sample names across the columns,
    and each later row is one metric. That row order IS the intended report order.
    """
    ws = load_workbook(xls, data_only=True).worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    header = [str(c).strip() for c in rows[0] if c is not None and str(c).strip()]
    # some sheets repeat the names in a "Sample name" row; the first row wins
    samples = [h for h in header if h.lower() not in ("metric", "metrics", "sample name")]
    metrics, seen = [], set()
    for r in rows[1:]:
        if not r or r[0] is None:
            continue
        label = str(r[0]).strip()
        if not label:
            continue
        # "Estimated number of cells" reappears as "Estimated number ofcells" — compare with all
        # whitespace removed so the duplicate is caught despite the typo.
        norm = "".join(label.lower().split())
        if not keep_metric(label, norm in seen):
            seen.add(norm)
            continue
        seen.add(norm)
        metrics.append((label, list(r[1:1 + len(samples)])))
    return samples, metrics


def pdf_to_pngs(pdf: Path, tag: str):
    """Render every page of a PDF to PNG; returns the page images in order."""
    TMP.mkdir(parents=True, exist_ok=True)
    out = TMP / tag
    for stale in TMP.glob(f"{tag}-*.png"):
        stale.unlink()
    subprocess.run(["pdftoppm", "-r", str(DPI), "-png", str(pdf), str(out)],
                   check=True, capture_output=True)
    return sorted(TMP.glob(f"{tag}-*.png"))


def add_slide(prs, title):
    s = prs.slides.add_slide(prs.slide_layouts[LAYOUT_TITLE_ONLY])
    if s.shapes.title is not None:
        s.shapes.title.text = title
    return s


def place_image(prs, slide, img: Path, top_in=1.35, margin_in=0.45):
    """Fit the image into the area under the title, preserving aspect ratio."""
    from PIL import Image
    avail_w = prs.slide_width / 914400 - 2 * margin_in
    avail_h = prs.slide_height / 914400 - top_in - 0.35
    with Image.open(img) as im:
        iw, ih = im.size
    scale = min(avail_w / (iw / DPI), avail_h / (ih / DPI))
    w, h = (iw / DPI) * scale, (ih / DPI) * scale
    slide.shapes.add_picture(str(img), Inches((prs.slide_width / 914400 - w) / 2),
                             Inches(top_in + (avail_h - h) / 2), Inches(w), Inches(h))


def placeholder(slide, msg):
    tb = slide.shapes.add_textbox(Inches(1), Inches(3), Inches(11), Inches(1))
    tb.text_frame.text = msg
    tb.text_frame.paragraphs[0].runs[0].font.size = Pt(18)
    print(f"  MISSING: {msg}")


def metrics_slide(prs, samples, metrics, title):
    s = add_slide(prs, title)
    rows, cols = len(metrics) + 1, len(samples) + 1
    # LibreOffice/PowerPoint grow a row when its label wraps, so the table must be sized to keep
    # every label on one line: wide first column, small type, and a row height that leaves the
    # whole table inside the 7.5in slide.
    top, avail_h = 1.12, 6.1
    row_h = min(0.26, avail_h / rows)
    tbl = s.shapes.add_table(rows, cols, Inches(0.3), Inches(top),
                             Inches(12.7), Inches(row_h * rows)).table
    tbl.columns[0].width = Inches(3.0)
    for c in range(1, cols):
        tbl.columns[c].width = Inches(9.7 / (cols - 1))
    tbl.cell(0, 0).text = "Metric"
    for j, nm in enumerate(samples, start=1):
        tbl.cell(0, j).text = str(nm)
    for i, (label, vals) in enumerate(metrics, start=1):
        tbl.cell(i, 0).text = label
        for j, v in enumerate(vals, start=1):
            if isinstance(v, float) and 0 <= v <= 1:
                txt = f"{v*100:.2f}%"
            elif isinstance(v, (int, float)):
                txt = f"{v:,.0f}"
            else:
                txt = "" if v is None else str(v)
            tbl.cell(i, j).text = txt
    for row in tbl.rows:
        row.height = Inches(row_h)
        for c in row.cells:
            c.margin_top = c.margin_bottom = Inches(0.01)
            for p in c.text_frame.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(7.5)
    return s


def summary_slide(prs, md: Path):
    """Render a markdown bullets (+ optional pipe table) file onto one slide.

    Conclusions are interpretive and run-specific, so they live in a file rather than in this
    script — another flowcell supplies its own, or omits --summary entirely.
    """
    title, bullets, rows = "Summary", [], []
    for line in md.read_text(encoding="utf-8").splitlines():
        t = line.strip()
        if t.startswith("#"):
            if not bullets and not rows:
                title = t.lstrip("#").strip()
        elif t.startswith(("-", "*")):
            bullets.append(re.sub(r"\*\*(.+?)\*\*", r"\1", t[1:].strip()))
        elif t.startswith("|"):
            cells = [c.strip() for c in t.strip("|").split("|")]
            if not all(set(c) <= set("-: ") for c in cells):      # skip the --- separator row
                rows.append(cells)
    s = add_slide(prs, title)
    tb = s.shapes.add_textbox(Inches(0.45), Inches(1.25), Inches(12.4), Inches(3.0))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, b in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = "• " + b
        p.font.size = Pt(11)
        p.space_after = Pt(3)
    if rows:
        tbl = s.shapes.add_table(len(rows), len(rows[0]), Inches(0.45),
                                 Inches(4.45), Inches(12.4),
                                 Inches(0.26 * len(rows))).table
        for r, row in enumerate(rows):
            for c, val in enumerate(row):
                tbl.cell(r, c).text = val
        for row in tbl.rows:
            row.height = Inches(0.26)
            for c in row.cells:
                c.margin_top = c.margin_bottom = Inches(0.01)
                for p in c.text_frame.paragraphs:
                    for run in p.runs:
                        run.font.size = Pt(9)
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--xls", required=True, type=Path)
    ap.add_argument("--template", required=True, type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--title", default=None)
    ap.add_argument("--summary", type=Path,
                    help="markdown file of findings bullets (+ optional table) for a summary slide")
    a = ap.parse_args()

    run = a.run_dir.resolve()
    out = a.out or run / "reports" / f"{run.name.replace('results_', '')}_report.pptx"
    out.parent.mkdir(parents=True, exist_ok=True)

    samples, metrics = sample_order_from_xls(a.xls)
    print(f"sample order from {a.xls.name}: {', '.join(samples)}")

    prs = Presentation(a.template)
    # The template ships a cover, an empty body slide and a closing "Thank you." slide. Keep the
    # cover, drop the other two, and re-add the closing slide at the end after the content.
    # Reordering via sldIdLst remove/append duplicated a slide part instead of moving it, so build
    # the deck in order rather than rearranging it afterwards.
    # Dropping the sldId alone leaves the slide part orphaned in the package, and python-pptx then
    # writes it alongside a new slide of the same partname ("Duplicate name: ppt/slides/slide2.xml").
    # The relationship must go too.
    xml_slides = prs.slides._sldIdLst
    for sld in list(xml_slides)[1:]:
        prs.part.drop_rel(sld.rId)
        xml_slides.remove(sld)

    title_txt = a.title or f"{run.name.replace('results_', '').replace('_filtered', '')} scRNA analysis"
    cover = list(prs.slides)[0]
    for sh in cover.shapes:
        if not sh.has_text_frame:
            continue
        t = sh.text_frame.text
        if "Title" in t:
            sh.text_frame.text = title_txt
        elif t.strip().upper().startswith("DATE"):
            sh.text_frame.text = datetime.date.today().strftime("%d %B %Y")

    metrics_slide(prs, samples, metrics, "Sequencing metrics")

    if a.summary:
        if a.summary.exists():
            summary_slide(prs, a.summary)
        else:
            print(f"  MISSING: summary file not found: {a.summary}")

    # UMAPs: one PDF page per sample pair, already in run order once step 06 orders by SAMPLE_NAMES
    upath = run / UMAP_PDF
    if upath.exists():
        pages = pdf_to_pngs(upath, "umap")
        for i, pg in enumerate(pages):
            pair = samples[i * 2:i * 2 + 2]
            s = add_slide(prs, " and ".join(pair) if pair else "Cell Types – Split by Sample")
            place_image(prs, s, pg)
    else:
        placeholder(add_slide(prs, "Cell Types – Split by Sample"), f"figure not found: {upath}")

    for rel, title in FIGURES:
        titles = title if isinstance(title, list) else None
        p = run / rel
        if not p.exists():
            placeholder(add_slide(prs, titles[0] if titles else title), f"figure not found: {p}")
            continue
        for i, pg in enumerate(pdf_to_pngs(p, re.sub(r"\W+", "_", rel))):
            t = titles[i] if titles and i < len(titles) else (titles[-1] if titles else title)
            place_image(prs, add_slide(prs, t), pg)

    # closing slide, on the template's blank layout
    closing = prs.slides.add_slide(prs.slide_layouts[4])
    tb = closing.shapes.add_textbox(Inches(5.2), Inches(3.2), Inches(3), Inches(1))
    tb.text_frame.text = "Thank you."
    tb.text_frame.paragraphs[0].runs[0].font.size = Pt(28)

    prs.save(out)
    print(f"wrote {out} ({len(list(prs.slides))} slides)")


if __name__ == "__main__":
    sys.exit(main())
