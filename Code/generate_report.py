import csv
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import mm, cm
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                  TableStyle, HRFlowable, KeepTogether, PageBreak)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT, TA_JUSTIFY
from reportlab.platypus import Flowable
from reportlab.lib.colors import HexColor

# ══════════════════════════════════════════════════════════════════════════════
# COLOR PALETTE
# ══════════════════════════════════════════════════════════════════════════════
C_NAVY      = HexColor("#0d1b2a")
C_BLUE      = HexColor("#1b4f72")
C_ACCENT    = HexColor("#2e86ab")
C_ACCENT2   = HexColor("#a23b72")
C_ORANGE    = HexColor("#e07b39")
C_GREEN     = HexColor("#1b7a4a")
C_GREEN_L   = HexColor("#d4edda")
C_PURPLE    = HexColor("#6a2d6a")
C_LIGHT     = HexColor("#f4f7fb")
C_LIGHTER   = HexColor("#eaf1f8")
C_MID       = HexColor("#c8d8e8")
C_BORDER    = HexColor("#b0c4d8")
C_TEXT      = HexColor("#1a2332")
C_SUBTEXT   = HexColor("#4a5568")
C_MUTED     = HexColor("#8898a8")
C_WHITE     = colors.white
C_RED_L     = HexColor("#fdecea")
C_RED       = HexColor("#c0392b")

PAGE_W = A4[0] - 2*cm

# ══════════════════════════════════════════════════════════════════════════════
# STYLES
# ══════════════════════════════════════════════════════════════════════════════
def S(name):
    styles = {
        "h1": ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=16,
                             textColor=C_NAVY, leading=20, spaceAfter=4),
        "h2": ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=12,
                             textColor=C_BLUE, leading=16, spaceBefore=2),
        "h3": ParagraphStyle("h3", fontName="Helvetica-Bold", fontSize=10,
                             textColor=C_ACCENT, leading=14),
        "body": ParagraphStyle("body", fontName="Helvetica", fontSize=9,
                               textColor=C_TEXT, leading=14, spaceAfter=4,
                               alignment=TA_JUSTIFY),
        "body_c": ParagraphStyle("body_c", fontName="Helvetica", fontSize=9,
                                 textColor=C_TEXT, leading=14, alignment=TA_CENTER),
        "small": ParagraphStyle("small", fontName="Helvetica", fontSize=8,
                                textColor=C_MUTED, leading=11),
        "small_i": ParagraphStyle("small_i", fontName="Helvetica-Oblique", fontSize=8,
                                  textColor=C_MUTED, leading=11),
        "th": ParagraphStyle("th", fontName="Helvetica-Bold", fontSize=8.5,
                             textColor=C_WHITE, leading=12),
        "th_dark": ParagraphStyle("th_dark", fontName="Helvetica-Bold", fontSize=8.5,
                                  textColor=C_NAVY, leading=12),
        "td": ParagraphStyle("td", fontName="Helvetica", fontSize=8.5,
                             textColor=C_TEXT, leading=12),
        "td_b": ParagraphStyle("td_b", fontName="Helvetica-Bold", fontSize=8.5,
                               textColor=C_NAVY, leading=12),
        "td_accent": ParagraphStyle("td_accent", fontName="Helvetica-Bold", fontSize=8.5,
                                    textColor=C_ACCENT, leading=12),
        "td_green": ParagraphStyle("td_green", fontName="Helvetica-Bold", fontSize=8.5,
                                   textColor=C_GREEN, leading=12),
        "td_muted": ParagraphStyle("td_muted", fontName="Helvetica", fontSize=8,
                                   textColor=C_MUTED, leading=11),
        "badge": ParagraphStyle("badge", fontName="Helvetica-Bold", fontSize=7.5,
                                textColor=C_WHITE, leading=10, alignment=TA_CENTER),
        "footer": ParagraphStyle("footer", fontName="Helvetica", fontSize=7.5,
                                 textColor=C_MUTED, alignment=TA_CENTER),
        "caption": ParagraphStyle("caption", fontName="Helvetica-Oblique", fontSize=8,
                                  textColor=C_SUBTEXT, leading=12, spaceAfter=6),
    }
    return styles[name]

# ══════════════════════════════════════════════════════════════════════════════
# CUSTOM FLOWABLES
# ══════════════════════════════════════════════════════════════════════════════

class CoverBanner(Flowable):
    """Full-width cover header with gradient effect"""
    def __init__(self, title, subtitle, meta_line, width=None):
        super().__init__()
        self.title     = title
        self.subtitle  = subtitle
        self.meta_line = meta_line
        self.width     = width or PAGE_W
        self.height    = 72

    def draw(self):
        c = self.canv
        w, h = self.width, self.height
        # Background layers for depth
        c.setFillColor(C_NAVY)
        c.roundRect(0, 0, w, h, 5, fill=1, stroke=0)
        c.setFillColor(HexColor("#142035"))
        c.roundRect(w*0.55, 0, w*0.45, h, 5, fill=1, stroke=0)
        # Accent stripe top
        c.setFillColor(C_ACCENT)
        c.rect(0, h-4, w*0.45, 4, fill=1, stroke=0)
        c.setFillColor(C_ACCENT2)
        c.rect(w*0.45, h-4, w*0.25, 4, fill=1, stroke=0)
        c.setFillColor(C_ORANGE)
        c.rect(w*0.7, h-4, w*0.3, 4, fill=1, stroke=0)
        # Title
        c.setFillColor(C_WHITE)
        c.setFont("Helvetica-Bold", 21)
        c.drawString(16, h - 30, self.title)
        # Subtitle
        c.setFillColor(HexColor("#a8c4d8"))
        c.setFont("Helvetica", 10)
        c.drawString(16, h - 45, self.subtitle)
        # Meta line
        c.setFillColor(HexColor("#6a8fa8"))
        c.setFont("Helvetica", 8)
        c.drawString(16, h - 58, self.meta_line)
        # Cell-Hub badge
        bx, by, bw, bh = w-90, h-26, 78, 18
        c.setFillColor(C_ACCENT)
        c.roundRect(bx, by, bw, bh, 4, fill=1, stroke=0)
        c.setFillColor(C_WHITE)
        c.setFont("Helvetica-Bold", 9)
        c.drawCentredString(bx + bw/2, by + 5, "Cell-Hub  v1.0")
        # Bottom meta
        c.setFillColor(HexColor("#4a6a82"))
        c.setFont("Helvetica", 7.5)
        c.drawRightString(w - 12, 5, self.meta_line)

    def wrap(self, aw, ah):
        self.width = aw
        return self.width, self.height


class SectionBanner(Flowable):
    """Colored section header"""
    def __init__(self, number, title, color=C_ACCENT, width=None):
        super().__init__()
        self.number = number
        self.title  = title
        self.color  = color
        self.width  = width or PAGE_W
        self.height = 24

    def draw(self):
        c = self.canv
        w, h = self.width, self.height
        c.setFillColor(self.color)
        c.roundRect(0, 0, w, h, 4, fill=1, stroke=0)
        # Number badge
        c.setFillColor(colors.white)
        c.setFillAlpha(0.18)
        c.circle(16, h/2, 10, fill=1, stroke=0)
        c.setFillAlpha(1)
        c.setFillColor(C_WHITE)
        c.setFont("Helvetica-Bold", 10)
        c.drawCentredString(16, h/2 - 4, str(self.number))
        # Title
        c.setFont("Helvetica-Bold", 11)
        c.drawString(34, h/2 - 4, self.title)

    def wrap(self, aw, ah):
        self.width = aw
        return self.width, self.height


class StatCard(Flowable):
    """Metric summary card"""
    def __init__(self, value, label, sublabel="", color=C_ACCENT, w=38*mm, h=20*mm):
        super().__init__()
        self.value    = str(value)
        self.label    = label
        self.sublabel = sublabel
        self.color    = color
        self.width    = w
        self.height   = h

    def draw(self):
        c = self.canv
        w, h = self.width, self.height
        # Card background
        c.setFillColor(C_LIGHT)
        c.setStrokeColor(self.color)
        c.setLineWidth(1.2)
        c.roundRect(0, 0, w, h, 4, fill=1, stroke=1)
        # Top color bar
        c.setFillColor(self.color)
        c.roundRect(0, h-5, w, 5, 3, fill=1, stroke=0)
        # Value
        c.setFillColor(self.color)
        font_size = 15 if len(self.value) <= 6 else 11
        c.setFont("Helvetica-Bold", font_size)
        c.drawCentredString(w/2, h/2 - 1, self.value)
        # Label
        c.setFillColor(C_SUBTEXT)
        c.setFont("Helvetica-Bold", 7)
        c.drawCentredString(w/2, 8, self.label)
        # Sublabel
        if self.sublabel:
            c.setFillColor(C_MUTED)
            c.setFont("Helvetica", 6)
            c.drawCentredString(w/2, 2, self.sublabel)

    def wrap(self, aw, ah):
        return self.width, self.height


class InlineBadge(Flowable):
    """Small colored badge for values like species, method"""
    def __init__(self, text, color=C_ACCENT, w=None, h=12):
        super().__init__()
        self.text   = text
        self.color  = color
        self.width  = w or (len(text) * 5.5 + 12)
        self.height = h

    def draw(self):
        c = self.canv
        c.setFillColor(self.color)
        c.roundRect(0, 0, self.width, self.height, 3, fill=1, stroke=0)
        c.setFillColor(C_WHITE)
        c.setFont("Helvetica-Bold", 7)
        c.drawCentredString(self.width/2, 3, self.text)

    def wrap(self, aw, ah):
        return self.width, self.height


class DividerLine(Flowable):
    def __init__(self, color=C_BORDER, thickness=0.5):
        super().__init__()
        self.color     = color
        self.thickness = thickness
        self.height    = 4

    def draw(self):
        self.canv.setStrokeColor(self.color)
        self.canv.setLineWidth(self.thickness)
        self.canv.line(0, self.height/2, self._availableWidth, self.height/2)

    def wrap(self, aw, ah):
        self._availableWidth = aw
        return aw, self.height


# ══════════════════════════════════════════════════════════════════════════════
# TABLE BUILDERS
# ══════════════════════════════════════════════════════════════════════════════

def param_table_2col(rows, widths=None, header_color=C_BLUE):
    """Simple 2-column Parameter / Value table"""
    if widths is None:
        widths = [90*mm, 80*mm]
    data = [[Paragraph("<b>Parameter</b>", S("th")),
             Paragraph("<b>Value</b>", S("th"))]]
    for param, value, *_ in rows:
        is_na = str(value).strip() in ("NA", "N/A", "unknown", "not_run", "None", "")
        val_style = S("td_muted") if is_na else S("td_accent")
        data.append([
            Paragraph(str(param), S("td_b")),
            Paragraph(str(value), val_style),
        ])
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND",    (0,0), (-1,0), header_color),
        ("ROWBACKGROUNDS",(0,1), (-1,-1), [C_WHITE, C_LIGHTER]),
        ("GRID",          (0,0), (-1,-1), 0.3, C_BORDER),
        ("TOPPADDING",    (0,0), (-1,-1), 5),
        ("BOTTOMPADDING", (0,0), (-1,-1), 5),
        ("LEFTPADDING",   (0,0), (-1,-1), 8),
        ("RIGHTPADDING",  (0,0), (-1,-1), 8),
        ("VALIGN",        (0,0), (-1,-1), "MIDDLE"),
        ("LINEAFTER",     (0,0), (0,-1), 0.5, C_ACCENT),
        ("ROUNDEDCORNERS",[3]),
    ]))
    return t


def param_table_3col(rows, widths=None, header_color=C_BLUE):
    """3-column Parameter / Value / Description table"""
    if widths is None:
        widths = [68*mm, 32*mm, 72*mm]
    data = [[Paragraph("<b>Parameter</b>", S("th")),
             Paragraph("<b>Value</b>", S("th")),
             Paragraph("<b>Note</b>", S("th"))]]
    for param, value, desc in rows:
        is_na = str(value).strip() in ("NA", "N/A", "unknown", "not_run", "None", "")
        val_style = S("td_muted") if is_na else S("td_accent")
        data.append([
            Paragraph(str(param), S("td_b")),
            Paragraph(str(value), val_style),
            Paragraph(str(desc), S("td_muted")),
        ])
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND",    (0,0), (-1,0), header_color),
        ("ROWBACKGROUNDS",(0,1), (-1,-1), [C_WHITE, C_LIGHTER]),
        ("GRID",          (0,0), (-1,-1), 0.3, C_BORDER),
        ("TOPPADDING",    (0,0), (-1,-1), 5),
        ("BOTTOMPADDING", (0,0), (-1,-1), 5),
        ("LEFTPADDING",   (0,0), (-1,-1), 8),
        ("RIGHTPADDING",  (0,0), (-1,-1), 8),
        ("VALIGN",        (0,0), (-1,-1), "MIDDLE"),
        ("LINEAFTER",     (0,0), (0,-1), 0.5, C_BORDER),
        ("LINEAFTER",     (1,0), (1,-1), 0.5, C_BORDER),
    ]))
    return t


def cluster_composition_table(data):
    """Rich cluster × dataset composition table with QC medians"""
    cluster_rows = [(r["parameter"].replace("Cluster: ", ""), r["value"], r["description"])
                    for r in data.get("Cluster Composition", [])
                    if r["parameter"].startswith("Cluster:")]
    
    dataset_cats  = sorted([k for k in data.keys() if k.startswith("Dataset Composition")])
    dataset_names = [c.replace("Dataset Composition - ", "") for c in dataset_cats]
    total         = next((r["value"] for r in data.get("Cluster Composition", [])
                          if r["parameter"] == "Total nuclei (post-clustering)"), "?")
    
    # Build header
    col_labels = ["Cluster", "N", "%"] + dataset_names + ["Median genes", "Median %mito"]
    col_w = [32*mm, 14*mm, 12*mm] + [22*mm]*len(dataset_names) + [24*mm, 22*mm]
    
    header_row = [Paragraph(f"<b>{c}</b>", S("th")) for c in col_labels]
    table_data = [header_row]
    
    for cl_name, n_val, desc in cluster_rows:
        n_int = int(n_val)
        pct   = round(n_int / int(total) * 100, 1)
        
        # Extract QC medians from description field
        med_genes = med_mt = "—"
        if "median genes:" in desc:
            med_genes = desc.split("median genes:")[1].split("|")[0].strip()
        if "median %mito:" in desc:
            med_mt = desc.split("median %mito:")[1].split("|")[0].strip()
        
        row = [
            Paragraph(f"<b>{cl_name}</b>",
                      ParagraphStyle("cn", fontName="Helvetica-Bold", fontSize=9,
                                     textColor=C_ACCENT2, leading=12)),
            Paragraph(str(n_int), S("td_b")),
            Paragraph(f"{pct}%", S("td")),
        ]
        for ds_cat, ds_name in zip(dataset_cats, dataset_names):
            n_ds = next((r["value"] for r in data.get(ds_cat, [])
                         if r["parameter"] == f"[{ds_name}] Cluster: {cl_name}"), "0")
            pct_ds = round(int(n_ds) / n_int * 100, 1) if n_int > 0 else 0
            if pct_ds > 60:
                color_hex = "#1b7a4a"
            elif pct_ds > 20:
                color_hex = "#2e86ab"
            else:
                color_hex = "#8898a8"
            row.append(Paragraph(
                f"<b>{n_ds}</b><br/><font size='7' color='{color_hex}'>{pct_ds}%</font>",
                S("td")))
        
        row += [Paragraph(str(med_genes), S("td")),
                Paragraph(str(med_mt), S("td"))]
        table_data.append(row)
    
    # Total row
    total_row = [Paragraph("<b>TOTAL</b>", S("td_b")),
                 Paragraph(f"<b>{total}</b>", S("td_b")),
                 Paragraph("<b>100%</b>", S("td_b"))]
    for ds_cat, ds_name in zip(dataset_cats, dataset_names):
        n_ds = next((r["value"] for r in data.get(ds_cat, [])
                     if r["parameter"] == f"[{ds_name}] Total nuclei"), "?")
        pct_ds = round(int(n_ds)/int(total)*100, 1) if total != "?" else "?"
        total_row.append(Paragraph(f"<b>{n_ds}</b><br/><font size='7'>{pct_ds}%</font>",
                                   S("td_b")))
    total_row += [Paragraph("—", S("td")), Paragraph("—", S("td"))]
    table_data.append(total_row)
    
    t = Table(table_data, colWidths=col_w, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND",    (0,0), (-1,0),  C_NAVY),
        ("ROWBACKGROUNDS",(0,1), (-1,-2), [C_WHITE, C_LIGHTER]),
        ("BACKGROUND",    (0,-1),(-1,-1), HexColor("#e8f0fe")),
        ("FONTNAME",      (0,-1),(-1,-1), "Helvetica-Bold"),
        ("GRID",          (0,0), (-1,-1), 0.3, C_BORDER),
        ("TOPPADDING",    (0,0), (-1,-1), 5),
        ("BOTTOMPADDING", (0,0), (-1,-1), 5),
        ("LEFTPADDING",   (0,0), (-1,-1), 6),
        ("RIGHTPADDING",  (0,0), (-1,-1), 6),
        ("VALIGN",        (0,0), (-1,-1), "MIDDLE"),
        ("LINEAFTER",     (2,0), (2,-1),  1.0, C_ACCENT),
    ]))
    return t


def stat_cards_row(cards, n_cols=None):
    """Row of StatCard flowables"""
    if n_cols is None:
        n_cols = len(cards)
    card_w = (PAGE_W - (n_cols-1)*3*mm) / n_cols
    card_h = 20*mm
    row    = [[StatCard(v, l, s, c, card_w, card_h) for v, l, s, c in cards]]
    t = Table(row, colWidths=[card_w]*n_cols, rowHeights=[card_h])
    t.setStyle(TableStyle([
        ("LEFTPADDING",  (0,0),(-1,-1), 1),
        ("RIGHTPADDING", (0,0),(-1,-1), 1),
        ("TOPPADDING",   (0,0),(-1,-1), 0),
        ("BOTTOMPADDING",(0,0),(-1,-1), 0),
    ]))
    return t


# ══════════════════════════════════════════════════════════════════════════════
# DATA HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def load_csv(path):
    data = {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            cat = row["category"]
            data.setdefault(cat, []).append(row)
    return data


def g(data, category, parameter, default="N/A"):
    for r in data.get(category, []):
        if r["parameter"] == parameter:
            v = r["value"]
            return v if v and v.strip() not in ("NA", "") else default
    return default


# ══════════════════════════════════════════════════════════════════════════════
# PAGE TEMPLATE  — now a closure to capture seed dynamically
# ══════════════════════════════════════════════════════════════════════════════
def make_on_page(seed):
    """Return a page callback that embeds the real seed value in the footer."""
    def on_page(canvas, doc):
        canvas.saveState()
        w, h = A4
        canvas.setFillColor(HexColor("#eaf1f8"))
        canvas.rect(cm, 0, w - 2*cm, 10*mm, fill=1, stroke=0)
        canvas.setStrokeColor(C_BORDER)
        canvas.setLineWidth(0.5)
        canvas.line(cm, 10*mm, w - cm, 10*mm)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(C_MUTED)
        canvas.drawString(cm + 2, 3.5*mm,
                          f"Cell-Hub Analysis Report  •  Generated automatically  •  Seed = {seed}")
        canvas.drawRightString(w - cm - 2, 3.5*mm, f"Page {doc.page}")
        canvas.setStrokeColor(C_ACCENT)
        canvas.setLineWidth(1.5)
        canvas.line(cm, h - 8*mm, cm + (w - 2*cm)*0.45, h - 8*mm)
        canvas.setStrokeColor(C_ACCENT2)
        canvas.line(cm + (w - 2*cm)*0.45, h - 8*mm, cm + (w - 2*cm)*0.7, h - 8*mm)
        canvas.setStrokeColor(C_ORANGE)
        canvas.line(cm + (w - 2*cm)*0.7, h - 8*mm, w - cm, h - 8*mm)
        canvas.restoreState()
    return on_page


# ══════════════════════════════════════════════════════════════════════════════
# REPORT BUILDER
# ══════════════════════════════════════════════════════════════════════════════
def build_report(csv_path, output_path, module="multiple"):
    data  = load_csv(csv_path)
    story = []

    export_ts    = g(data, "Reproducibility",         "Export timestamp")
    n_total      = g(data, "Cluster Composition",     "Total nuclei (post-clustering)")
    n_clusters   = g(data, "Cluster Composition",     "Number of clusters")
    n_datasets   = g(data, "Integration",             "Number of datasets integrated")
    n_pcs        = g(data, "Dimensionality Reduction","Number of PCs (FindNeighbors / UMAP)")
    resolution   = g(data, "Clustering",              "Clustering resolution")
    k_param      = g(data, "Dimensionality Reduction","k neighbors (FindNeighbors k.param)")
    algo         = g(data, "Clustering",              "Clustering algorithm")
    norm_method  = g(data, "Normalization",           "Normalization method")
    n_var_feat   = g(data, "Normalization",           "Variable features (n)")
    integ_method = g(data, "Integration",             "Integration method")
    ds_names     = g(data, "Integration",             "Dataset names")
    seed         = g(data, "Reproducibility",         "Random seed", "42")

    # Cover subtitle adapts to module
    if module == "single":
        cover_subtitle = "Single Dataset Analysis  •  Cell-Hub Pipeline"
    else:
        cover_subtitle = "Multi-dataset Integration  •  Cell-Hub Pipeline"

    # ════════════════════════════════════════════════════════════════════════
    # PAGE 1 — COVER + SUMMARY
    # ════════════════════════════════════════════════════════════════════════
    story.append(CoverBanner(
        title     = "snRNA-seq Analysis Report",
        subtitle  = cover_subtitle,
        meta_line = f"Exported: {export_ts}  •  Seed: {seed}",
    ))
    story.append(Spacer(1, 5*mm))

    # Summary stat cards — hide n_datasets card for single module
    if module == "single":
        cards = [
            (n_total,    "Total nuclei",   "post-QC & clustering", C_ACCENT),
            (n_clusters, "Clusters",       "identified",           C_ACCENT2),
            (n_pcs,      "PCs used",       "FindNeighbors/UMAP",   C_GREEN),
            (resolution, "Resolution",     "FindClusters",         C_PURPLE),
            (k_param,    "k neighbors",    "SNN graph",            C_BLUE),
        ]
    else:
        cards = [
            (n_total,    "Total nuclei",   "post-QC & clustering", C_ACCENT),
            (n_clusters, "Clusters",       "identified",           C_ACCENT2),
            (n_datasets, "Datasets",       "integrated",           C_ORANGE),
            (n_pcs,      "PCs used",       "FindNeighbors/UMAP",   C_GREEN),
            (resolution, "Resolution",     "FindClusters",         C_PURPLE),
            (k_param,    "k neighbors",    "SNN graph",            C_BLUE),
        ]
    story.append(stat_cards_row(cards))
    story.append(Spacer(1, 5*mm))

    story.append(Paragraph(
        "This report summarizes the key analytical parameters, quality control metrics, "
        "and cluster composition of the snRNA-seq dataset processed through <b>Cell-Hub</b>. "
        "All stochastic steps (PCA, UMAP, neighbor graph, clustering) were performed with "
        f"a fixed random seed of <b>{seed}</b> to ensure full reproducibility. "
        "Parameters are recorded at the time of export and reflect the current state of the "
        "Seurat object.",
        S("body")
    ))
    story.append(Spacer(1, 4*mm))

    pipeline_text = [
        ["Step", "Method / Value"],
        ["Normalization",           norm_method],
        ["Variable features",       n_var_feat],
        ["Integration",             integ_method if module == "multiple" else "N/A (single dataset)"],
        ["Dimensions (PCA → UMAP)", n_pcs],
        ["k neighbors (SNN)",       k_param],
        ["Clustering algorithm",    algo],
        ["Clustering resolution",   resolution],
    ]
    pipeline_table_data = []
    for i, (step, val) in enumerate(pipeline_text):
        if i == 0:
            pipeline_table_data.append([
                Paragraph(f"<b>{step}</b>", S("th")),
                Paragraph(f"<b>{val}</b>",  S("th"))
            ])
        else:
            is_na = str(val).strip() in ("NA", "N/A", "unknown", "not_run", "None", "")
            pipeline_table_data.append([
                Paragraph(step, S("td_b")),
                Paragraph(val,  S("td_muted") if is_na else S("td_accent")),
            ])
    pt = Table(pipeline_table_data, colWidths=[85*mm, 85*mm], repeatRows=1)
    pt.setStyle(TableStyle([
        ("BACKGROUND",    (0,0), (-1,0), C_BLUE),
        ("ROWBACKGROUNDS",(0,1), (-1,-1), [C_WHITE, C_LIGHTER]),
        ("GRID",          (0,0), (-1,-1), 0.3, C_BORDER),
        ("TOPPADDING",    (0,0), (-1,-1), 5),
        ("BOTTOMPADDING", (0,0), (-1,-1), 5),
        ("LEFTPADDING",   (0,0), (-1,-1), 10),
        ("RIGHTPADDING",  (0,0), (-1,-1), 10),
        ("VALIGN",        (0,0), (-1,-1), "MIDDLE"),
        ("LINEAFTER",     (0,0), (0,-1), 0.5, C_BORDER),
    ]))
    story.append(KeepTogether([
        Paragraph("Pipeline Summary", S("h2")),
        Spacer(1, 2*mm),
        pt,
    ]))

    # ════════════════════════════════════════════════════════════════════════
    # PAGE 2 — QC + INTEGRATION
    # ════════════════════════════════════════════════════════════════════════
    story.append(SectionBanner(1, "Quality Control", C_ACCENT))
    story.append(Spacer(1, 3*mm))
    story.append(Paragraph(
        "Quality control was applied to remove low-quality nuclei based on the number "
        "of detected genes (<i>nFeature_RNA</i>), total UMI counts (<i>nCount_RNA</i>), "
        "and the fraction of mitochondrial transcripts (<i>percent.mt</i>). "
        "Doublets were detected using <b>DoubletFinder</b>.",
        S("body")
    ))
    story.append(Spacer(1, 2*mm))
    warning_data = [[Paragraph(
        "<b>⚠ Note:</b> QC threshold values below reflect the state of the interface sliders "
        "at the time of export. If parameters were adjusted after running QC, these values "
        "may not correspond to the thresholds actually applied to the data. "
        "Refer to the QC Statistics section for the actual pre/post-QC cell counts.",
        ParagraphStyle("warn", fontName="Helvetica", fontSize=8,
                       textColor=HexColor("#7d4e00"), leading=12)
    )]]
    warning_table = Table(warning_data, colWidths=[PAGE_W])
    warning_table.setStyle(TableStyle([
        ("BACKGROUND",    (0,0), (-1,-1), HexColor("#fff8e1")),
        ("LEFTPADDING",   (0,0), (-1,-1), 10),
        ("RIGHTPADDING",  (0,0), (-1,-1), 10),
        ("TOPPADDING",    (0,0), (-1,-1), 7),
        ("BOTTOMPADDING", (0,0), (-1,-1), 7),
        ("LINEBEFORE",    (0,0), (0,-1),  3, HexColor("#f0a500")),
        ("ROUNDEDCORNERS",[3]),
    ]))
    story.append(warning_table)
    story.append(Spacer(1, 3*mm))
    
    qc_thresh_rows = [
        ("Min nCount_RNA",      g(data,"QC Thresholds","Min nCount_RNA"),      "Minimum total UMI counts per nucleus"),
        ("Max nCount_RNA",      g(data,"QC Thresholds","Max nCount_RNA"),      "Maximum total UMI counts per nucleus"),
        ("Min nFeature_RNA",    g(data,"QC Thresholds","Min nFeature_RNA"),    "Minimum unique genes per nucleus"),
        ("Max nFeature_RNA",    g(data,"QC Thresholds","Max nFeature_RNA"),    "Maximum unique genes per nucleus"),
        ("Max mitochondrial %", g(data,"QC Thresholds","Max mitochondrial %"), "Upper limit for mitochondrial fraction"),
    ]

    qc_stat_cats = [k for k in data.keys() if k.startswith("QC Statistics")]

    story.append(KeepTogether([
        Paragraph("QC Thresholds Applied", S("h3")),
        Spacer(1, 2*mm),
        param_table_3col(qc_thresh_rows, widths=[68*mm, 28*mm, 76*mm], header_color=C_ACCENT),
    ]))
    story.append(Spacer(1, 4*mm))

    if qc_stat_cats:
        for cat in qc_stat_cats:
            ds_label = cat.replace("QC Statistics - ", "").replace("QC Statistics", "Global")
            rows = [(r["parameter"].replace(f"[{ds_label}] ", ""), r["value"], r["description"])
                    for r in data[cat]]
            story.append(KeepTogether([
                Paragraph(f"QC Statistics — {ds_label}", S("h3")),
                Spacer(1, 2*mm),
                param_table_3col(rows, widths=[72*mm, 28*mm, 72*mm], header_color=HexColor("#1a5276")),
                Spacer(1, 3*mm),
            ]))

    df_rows = [(r["parameter"], r["value"], r["description"])
               for r in data.get("DoubletFinder", [])]
    if df_rows:
        story.append(Spacer(1, 2*mm))
        story.append(KeepTogether([
            Paragraph("DoubletFinder Parameters", S("h3")),
            Spacer(1, 2*mm),
            param_table_3col(df_rows, widths=[72*mm, 28*mm, 72*mm], header_color=HexColor("#7d3c98")),
        ]))
        story.append(Spacer(1, 4*mm))

    # Integration section — only for multiple datasets
    if module == "multiple":
        story.append(SectionBanner(2, "Integration & Normalisation", C_ACCENT2))
        story.append(Spacer(1, 3*mm))
        story.append(Paragraph(
            f"A total of <b>{n_datasets}</b> dataset(s) were processed and integrated: "
            f"<b>{ds_names}</b>. "
            "Each dataset was independently normalised and highly variable features were identified "
            "before integration. The integration strategy is indicated below.",
            S("body")
        ))
        story.append(Spacer(1, 3*mm))
        int_rows = [
            ("Number of datasets",    n_datasets,   "Total datasets merged"),
            ("Dataset names",         ds_names,     ""),
            ("Integration method",    integ_method, "Batch correction strategy"),
            ("Normalization method",  norm_method,  "Applied per dataset"),
            ("Variable features (n)", n_var_feat,   "HVG selected for PCA"),
            ("Scale factor",          g(data,"Normalization","Scale factor"), "LogNormalize only"),
        ]
        harmony_vars = g(data, "Integration", "Harmony batch correction variable", None)
        if harmony_vars and harmony_vars != "N/A":
            int_rows.append(("Harmony variable", harmony_vars, "Metadata used for batch correction"))
        story.append(param_table_3col(int_rows, widths=[68*mm, 40*mm, 64*mm], header_color=C_ACCENT2))
    else:
        # For single dataset: just show normalization params
        story.append(SectionBanner(2, "Normalisation", C_ACCENT2))
        story.append(Spacer(1, 3*mm))
        norm_rows = [
            ("Normalization method",  norm_method, "Applied to raw counts"),
            ("Variable features (n)", n_var_feat,  "HVG selected for PCA"),
            ("Scale factor",          g(data,"Normalization","Scale factor"), "LogNormalize only"),
        ]
        story.append(param_table_3col(norm_rows, widths=[68*mm, 40*mm, 64*mm], header_color=C_ACCENT2))

    # ════════════════════════════════════════════════════════════════════════
    # PAGE 3 — CLUSTERING + COMPOSITION
    # ════════════════════════════════════════════════════════════════════════
    story.append(PageBreak())
    story.append(SectionBanner(3, "Dimensionality Reduction & Clustering", C_GREEN))
    story.append(Spacer(1, 3*mm))
    story.append(Paragraph(
        f"Principal component analysis was performed, and the top <b>{n_pcs} PCs</b> were "
        f"retained for downstream neighbour graph construction (<b>k = {k_param}</b>) and "
        "UMAP embedding. Clustering was performed using the "
        f"<b>{algo}</b> algorithm at resolution <b>{resolution}</b>.",
        S("body")
    ))
    story.append(Spacer(1, 3*mm))

    umap_nn = g(data,"Dimensionality Reduction","UMAP n.neighbors","30")
    umap_md = g(data,"Dimensionality Reduction","UMAP min.dist","0.3")
    dr_rows = [
        ("PCs used (FindNeighbors / UMAP)", n_pcs,      "Top PCs retained"),
        ("k neighbors (SNN graph)",         k_param,    "k.param in FindNeighbors"),
        ("UMAP n.neighbors",                umap_nn,    "Local neighborhood size"),
        ("UMAP min.dist",                   umap_md,    "Minimum distance in UMAP space"),
        ("Clustering algorithm",            algo,       "Graph partitioning method"),
        ("Clustering resolution",           resolution, "Higher = more clusters"),
        ("Random seed",                     seed,       "Fixed for reproducibility"),
    ]
    story.append(param_table_3col(dr_rows, widths=[72*mm, 30*mm, 70*mm], header_color=C_GREEN))
    story.append(Spacer(1, 5*mm))

    story.append(SectionBanner(4, "Cluster Composition", C_PURPLE))
    story.append(Spacer(1, 3*mm))
    cl_names_str = g(data, "Cluster Composition", "Cluster names (all)", "")
    story.append(Paragraph(
        f"A total of <b>{n_clusters} clusters</b> were identified across "
        f"<b>{n_total} nuclei</b>. "
        f"Cluster labels: <b>{cl_names_str}</b>. "
        "The table below shows the composition of each cluster, including "
        "the contribution per dataset and median QC metrics.",
        S("body")
    ))
    story.append(Spacer(1, 3*mm))
    story.append(Paragraph(
        "Values in the dataset columns indicate absolute nucleus count with percentage "
        "of the cluster in parenthesis. Green = majority contributor (&gt;60%).",
        S("caption")
    ))
    story.append(cluster_composition_table(data))

    # ════════════════════════════════════════════════════════════════════════
    # PAGE 4 — REPRODUCIBILITY + NOTES
    # ════════════════════════════════════════════════════════════════════════
    story.append(PageBreak())
    story.append(SectionBanner(5, "Reproducibility & Methods Notes", C_NAVY))
    story.append(Spacer(1, 4*mm))
    story.append(Paragraph("Reproducibility", S("h2")))
    story.append(Paragraph(
        f"All stochastic steps were performed with a fixed random seed of <b>{seed}</b> "
        "to ensure full reproducibility of results. This includes PCA "
        "(<code>RunPCA</code>), UMAP embedding (<code>RunUMAP</code>), "
        "neighbour graph construction (<code>FindNeighbors</code>), and "
        "graph-based clustering (<code>FindClusters</code>).",
        S("body")
    ))
    story.append(Spacer(1, 4*mm))
    story.append(Paragraph("Software & References", S("h2")))
    refs = [
        ("Seurat v5",     "Hao et al., Cell 2024. DOI: 10.1016/j.cell.2021.04.048"),
        ("Harmony",       "Korsunsky et al., Nat Methods 2019. DOI: 10.1038/s41592-019-0619-0"),
        ("DoubletFinder", "McGinnis et al., Cell Systems 2019. DOI: 10.1016/j.cels.2019.03.003"),
        ("Monocle 3",     "Cao et al., Nature 2019. DOI: 10.1038/s41586-019-0969-x"),
        ("Cell-Hub",      "Cell-Hub GitHub repository — GaspardMacaux/CellHub"),
    ]
    for tool, ref in refs:
        story.append(Paragraph(f"<b>{tool}:</b> {ref}", S("body")))
    story.append(Spacer(1, 4*mm))
    story.append(DividerLine())
    story.append(Spacer(1, 3*mm))
    story.append(Paragraph(
        f"Report generated by Cell-Hub on {export_ts}. "
        "Parameter values reflect the state of the Seurat object at export time. "
        "For questions or issues, refer to the Cell-Hub documentation on Read the Docs.",
        S("small_i")
    ))

    doc = SimpleDocTemplate(
        output_path, pagesize=A4,
        leftMargin=cm, rightMargin=cm,
        topMargin=1.2*cm, bottomMargin=1.5*cm,
        title="Cell-Hub snRNA-seq Analysis Report",
        author="Cell-Hub",
    )
    on_page = make_on_page(seed)
    doc.build(story, onFirstPage=on_page, onLaterPages=on_page)
    print(f"PDF generated: {output_path}")


# ══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--module", default="single")
    args = parser.parse_args()
    build_report(args.csv, args.output)


