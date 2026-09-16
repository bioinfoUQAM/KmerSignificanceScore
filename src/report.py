"""Automated visual report of a KSS run (Reviewer 4.5).

Turns the raw JSON scores into a printable PDF: a header page with the run summary and the
top-N positions table (every variant of each position, in the S4-Table layout), then, per gene
and laid out like Fig. 3, the KSS track with the top positions labelled above the discriminative
and mutational scores mirrored about a shared zero, with the protein's UniProt link. Built from
a run's <gene>/results/<gene>_compiled_results.json, so it works on any dataset.
"""

from __future__ import annotations

import json
import os
import re
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea

try:
    from . import kss, protein_score
except ImportError:
    import kss
    import protein_score

# Hues matching Figure 3 of the manuscript.
COLORS = {"discriminative": "#E07B54", "mutational": "#5BA85B", "kss": "#2C6FBB"}
UNIPROT_URL = "https://www.uniprot.org/uniprotkb/{}"
ROWS_PAGE1 = 20   # variant rows on the first table page (shares room with the legend)
ROWS_CONT = 32    # variant rows on a continuation page
TABLE_COLS = ["Rank", "Gene", "Pos", "KSS", "Disc", "Mut", "Prot", "Ref k-mer", "Variant",
              "AA change", "Top class"]
# Three columns hold free text and the rest hold short numbers, so the width goes where the
# text is: a k-mer lengthened by an insertion, the changes with their impacts, and the class
# distribution. The numeric columns are trimmed to what "0.997" and a gene name need, and the
# widths now sum to the full axes rather than leaving six percent of it unused.
TABLE_WIDTHS = [0.026, 0.038, 0.038, 0.042, 0.042, 0.042, 0.042, 0.088, 0.135, 0.280, 0.227]
TABLE_FONT = 6.0
TABLE_FIGURE_INCHES = 11.0     # figure width the table pages are drawn at
TABLE_WIDTH_FRACTION = 0.97    # share of that width the table occupies
CHARACTER_WIDTH_RATIO = 0.58   # average glyph width as a fraction of the font size, sans-serif
VARIANT_COL = 8  # index of the "Variant" column, whose changed nucleotides are coloured
DIFF_COLOR = "#C0392B"
CITATION = (
    "How to cite: Lebatteux D, Corso F, Soudeyns H, Boucoiran I, Gantt S, Diallo AB. "
    "KmerSignificance Score: a discriminative and biologically informed\n"
    "framework for viral k-mer prioritization.  "
    "https://github.com/bioinfoUQAM/KmerSignificanceScore")


# --------------------------------------------------------------------------------------
# Amino acid change labels (same conventions as Fig. 3: Delta for indels, grouped by position)
# --------------------------------------------------------------------------------------
def _has_indel(change: str) -> bool:
    """True when a change label carries a gap, which is how indels are written here."""
    return "-" in change


def _position_number(change: str) -> int:
    """The residue number inside a change label, 0 when the label carries none."""
    match = re.search(r"\d+", change)
    return int(match.group()) if match else 0


def _compact_label(changes: Dict[str, float], with_impacts: bool = True) -> str:
    """"DY144- (1.000) | E484A/K/Q (0.55/0.42/0.31)"; "Synonymous" when there are none.

    Each change carries its mutational impact, in the format the supplementary tables use, so
    that the column can be read against the position's mutational score. That score is the
    maximum over these changes, and without the values there is no way to see which change it
    came from: a position showing two changes told the reader nothing about which one drove it.
    """
    if not changes:
        return "Synonymous"
    scores = changes
    changes = list(changes)
    groups: "OrderedDict[int, List[str]]" = OrderedDict()
    for change in sorted(changes, key=_position_number):
        groups.setdefault(_position_number(change), []).append(change)
    parts = []
    for variants in groups.values():
        prefix = "Δ" if any(_has_indel(c) for c in variants) else ""
        # One value when they all share it, which indels routinely do: seven repetitions of
        # 1.00 carry no more information than one and push the column over its neighbour.
        values = [f"{scores.get(v, 0.0):.2f}" for v in variants]
        collapsed = values[0] if len(set(values)) == 1 else "/".join(values)
        # On the gene pages the bars are the scores, so the label only has to say which change
        # it is; printing the number beside a bar of that height says it twice.
        impacts = f" ({collapsed})" if with_impacts else ""
        if len(variants) == 1:
            parts.append(f"{prefix}{variants[0]}{impacts}")
        else:
            m = re.match(r"([A-Z*-])(\d+)(.+)", variants[0])
            if m:
                ref, pos, first_new = m.groups()
                news = [first_new] + [re.match(r"[A-Z*-]\d+(.+)", c).group(1)
                                      if re.match(r"[A-Z*-]\d+(.+)", c) else c for c in variants[1:]]
                parts.append(f"{prefix}{ref}{pos}{'/'.join(news)}{impacts}")
            else:
                parts.append(f"{prefix}{'/'.join(variants)}{impacts}")
    return " | ".join(parts)


def _distribution(class_counts: Dict[str, int], top: int = 3) -> str:
    """"A1 (60%), C (25%), D (15%), +2 others": per-class share among this variant's carriers."""
    counts = {c: n for c, n in class_counts.items() if n}
    total = sum(counts.values())
    if not total:
        return ""
    ordered = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    parts = []
    for c, n in ordered[:top]:
        pct = f"{n / total:.0%}"
        parts.append(f"{c} ({'<1%' if pct == '0%' else pct})")
    text = ", ".join(parts)
    if len(ordered) > top:
        text += f", +{len(ordered) - top} others"
    return text


# --------------------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------------------
def _load_gene_scores(results_dir: str, gene: str) -> Optional[Dict[str, Any]]:
    """The compiled results of one gene, or None when that gene has not been scored."""
    path = os.path.join(results_dir, gene, "results", f"{gene}_compiled_results.json")
    if not os.path.exists(path):
        return None
    return json.load(open(path, encoding="utf-8")).get(gene)


def _variant_rows(rec: dict) -> List[dict]:
    """One row per variant k-mer of a position, most frequent first, reference last."""
    ref = rec.get("ref", "")
    rows = []
    for kmer, data in rec.get("alts", {}).items():
        # The mapping, not its keys: the impact of each change is what the label reports.
        changes = data.get("amino_acid_changes", {})
        is_ref = kmer == ref
        rows.append({
            "variant": kmer, "is_ref": is_ref,
            "count": sum(data.get("class_counts", {}).values()),
            # The reference k-mer carries no change (kept for its class frequencies); a
            # different k-mer with no change is a genuine synonymous variant.
            "label": "reference" if is_ref else _compact_label(changes),
            # The same label without the impacts, for the gene pages where the bars carry them.
            "label_plain": "reference" if is_ref else _compact_label(changes, with_impacts=False),
            "indel": "yes" if any(_has_indel(c) for c in changes) else "",
            "distribution": _distribution(data.get("class_counts", {})),
        })
    rows.sort(key=lambda r: (r["is_ref"], -r["count"]))
    return rows


def _top_positions(scores_by_gene: Dict[str, dict], top_n: int) -> List[dict]:
    """Top-N positions across all genes by KSS, each with its variant rows."""
    ranked = []
    for gene, positions in scores_by_gene.items():
        for pos, rec in positions.items():
            ranked.append({
                "Gene": gene, "Position": int(pos), "KSS": rec.get("kss", 0.0),
                "Disc": rec.get("discriminative_score", 0.0),
                "Mut": rec.get("mutational_score", 0.0),
                "Prot": rec.get("protein_score", 0.0),
                "Ref": rec.get("ref", ""), "variants": _variant_rows(rec),
            })
    ranked.sort(key=lambda r: r["KSS"], reverse=True)
    return ranked[:top_n]


def _protein_info(input_folder: str, gene: str) -> Dict[str, str]:
    """UniProt name/accession/link and the GenBank reference for a gene (offline if cached)."""
    info = {"name": "protein annotation unavailable", "accession": "", "link": "",
            "genbank": "", "genbank_link": ""}
    try:
        gb_dir = os.path.join(input_folder, gene)
        gb = next(f for f in os.listdir(gb_dir) if f.lower().endswith(".gb"))
        accession = os.path.splitext(gb)[0]
        info["genbank"] = accession
        info["genbank_link"] = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
        entry = protein_score.fetch_uniprot_data(kss.get_taxon_id(os.path.join(gb_dir, gb)), gene)
        if not entry or entry.startswith(("No", "Error")):
            info["name"] = "UniProt entry not found"
            return info
        cache = os.path.join(os.path.dirname(os.path.abspath(protein_score.__file__)),
                             "uniprot", f"{entry}.json")
        data = json.load(open(cache, encoding="utf-8"))
        record = data["results"][0] if "results" in data else data
        info["name"] = (record.get("proteinDescription", {}).get("recommendedName", {})
                        .get("fullName", {}).get("value", "unnamed protein"))
        info["accession"] = entry
        info["link"] = UNIPROT_URL.format(entry)
        return info
    except Exception:
        return info


# --------------------------------------------------------------------------------------
# Pages
# --------------------------------------------------------------------------------------
def _paginate(top: List[dict]) -> List[List[dict]]:
    """Split positions into pages by variant-row budget, never splitting a position."""
    pages, page, used, budget = [], [], 0, ROWS_PAGE1
    for r in top:
        n = len(r["variants"])
        if page and used + n > budget:
            pages.append(page); page, used, budget = [], 0, ROWS_CONT
        page.append(r); used += n
    if page:
        pages.append(page)
    return pages


def _colour_variant(ax, cell, variant: str, ref: str, budget: int) -> None:
    """Redraw a variant k-mer in the cell, nucleotides that differ from the reference in red.

    This overlay replaces the cell's own text, so it is the one place where the k-mer's length
    has to be respected: an insertion can produce a k-mer of 27 nucleotides, and drawn in full
    it runs straight over the next column. Beyond the budget the middle is dropped rather than
    the tail, both ends of a k-mer carrying information, and each character kept is still
    compared against the reference at the position it actually occupies.
    """
    cell.get_text().set_text("")
    if len(variant) <= budget:
        kept = list(enumerate(variant))
    else:
        head = (budget - 1) // 2
        tail = budget - 1 - head
        kept = ([(i, variant[i]) for i in range(head)]
                + [(None, "…")]
                + [(i, variant[i]) for i in range(len(variant) - tail, len(variant))])

    chars = [TextArea(ch, textprops=dict(
        family="monospace", fontsize=6,
        color="#777777" if i is None
        else (DIFF_COLOR if (i >= len(ref) or ch != ref[i]) else "#111111")))
        for i, ch in kept]
    packed = HPacker(children=chars, align="center", pad=0, sep=0)
    ax.add_artist(AnnotationBbox(
        packed, (cell.get_x() + cell.get_width() / 2, cell.get_y() + cell.get_height() / 2),
        xycoords="axes fraction", frameon=False, pad=0, box_alignment=(0.5, 0.5)))


def _shorten(text: str, budget: int, middle: bool = False) -> str:
    """Bring a cell's text inside its column, marking the cut with an ellipsis.

    A matplotlib table neither wraps nor clips: text longer than its column is drawn straight
    over the neighbouring one, which is how a report looks right on one dataset and overruns on
    the next. Measuring the drawn text does not help here, because a table's cells are laid out
    after the text is set, so the budget is taken from the column width instead and applied to
    the string itself: deterministic, and independent of when the figure happens to be drawn.

    `middle` keeps both ends, which is what a k-mer needs: an insertion can make one 27
    nucleotides long, and its first and last codons say more than its first fourteen characters.
    """
    if len(text) <= budget:
        return text
    if not middle:
        return text[:budget - 1] + "…"
    head = (budget - 1) // 2
    tail = budget - 1 - head
    return text[:head] + "…" + text[len(text) - tail:]


def _budgets() -> List[int]:
    """Characters each column can hold, from its width at the table's font size."""
    inches = TABLE_WIDTH_FRACTION * TABLE_FIGURE_INCHES
    per_character = TABLE_FONT * CHARACTER_WIDTH_RATIO / 72.0
    return [max(4, int(width * inches / per_character)) for width in TABLE_WIDTHS]


def _render_table(ax, positions: List[dict], start_rank: int) -> None:
    """Draw the variant table with the shared columns truly merged per position."""
    cells, blocks = [], []          # blocks: (row_start, row_end, position, rank)
    budgets = _budgets()
    for offset, r in enumerate(positions):
        row_start = len(cells)
        for var in r["variants"]:
            cells.append(["", "", "", "", "", "", "", "",
                          _shorten(var["variant"], budgets[8], middle=True),
                          _shorten(var["label"], budgets[9]),
                          _shorten(var["distribution"], budgets[10])])
        blocks.append((row_start, len(cells) - 1, r, start_rank + offset))

    table = ax.table(cellText=cells, colLabels=TABLE_COLS, loc="upper center",
                     cellLoc="center", colWidths=TABLE_WIDTHS)
    table.auto_set_font_size(False); table.set_fontsize(TABLE_FONT); table.scale(1, 1.25)
    for j in range(len(TABLE_COLS)):
        table[0, j].set_facecolor("#2A2A2A"); table[0, j].set_text_props(color="white")

    # Merge the shared columns: value on the block's middle row, internal borders removed.
    for row_start, row_end, r, rank in blocks:
        shared = {0: rank, 1: r["Gene"], 2: r["Position"], 3: f"{r['KSS']:.3f}",
                  4: f"{r['Disc']:.3f}", 5: f"{r['Mut']:.3f}", 6: f"{r['Prot']:.3f}",
                  7: _shorten(r["Ref"], budgets[7], middle=True)}
        middle = (row_start + row_end) // 2
        for j, value in shared.items():
            for row in range(row_start, row_end + 1):
                cell = table[row + 1, j]     # +1 skips the header row
                cell.get_text().set_text(str(value) if row == middle else "")
                if row_start == row_end:
                    edges = "LRTB"
                elif row == row_start:
                    edges = "LRT"
                elif row == row_end:
                    edges = "LRB"
                else:
                    edges = "LR"
                cell.visible_edges = edges


    # Cell positions are only computed at draw time; force it, then colour the variant cells.
    ax.figure.canvas.draw()
    for row_start, _row_end, r, _rank in blocks:
        for k, var in enumerate(r["variants"]):
            _colour_variant(ax, table[row_start + 1 + k, VARIANT_COL], var["variant"],
                            r["Ref"], budgets[VARIANT_COL])


def _table_pages(pdf, dataset_name, weights, scores_by_gene, top):
    """Page 1 (title, run summary, component legend, first table chunk) then continuations."""
    pages = _paginate(top)
    rank = 1
    for page_index, chunk in enumerate(pages):
        fig = plt.figure(figsize=(11, 8.5))
        if page_index == 0:
            fig.suptitle(f"KSS report - {dataset_name}", fontsize=16, fontweight="bold", y=0.97)
            n_pos = sum(len(p) for p in scores_by_gene.values())
            fig.text(0.06, 0.90,
                     f"{len(scores_by_gene)} gene(s), {n_pos} scored positions.  "
                     f"Weights: discriminative {weights['discriminative']}, "
                     f"mutational {weights['mutational']}, protein {weights['protein']}.",
                     fontsize=10)
            fig.text(0.06, 0.87,
                     "KSS = weighted mean of three components, each in [0, 1]:\n"
                     "  • Discriminative - class separation, from information theory "
                     "(normalized mutual information × k-mer purity).\n"
                     "  • Mutational - biophysical impact of the amino acid changes, from the "
                     "optimized MIYATA_EVO substitution matrix.\n"
                     "  • Protein - depth of UniProt characterization, a weighted sum of ten "
                     "annotation-evidence categories.",
                     fontsize=8.5, va="top", family="monospace")
            fig.text(0.5, 0.03, CITATION, ha="center", va="bottom", fontsize=7.5,
                     color="#555555")
            rect, title = [0.015, 0.04, 0.97, 0.64], f"Top {len(top)} positions by KSS (all variants)"
        else:
            fig.suptitle(f"KSS report - {dataset_name}  (top positions, continued)",
                         fontsize=13, fontweight="bold", y=0.97)
            rect, title = [0.015, 0.04, 0.97, 0.88], "Top positions by KSS (continued)"
        ax = fig.add_axes(rect); ax.axis("off")
        ax.set_title(title, fontsize=11, loc="left", pad=6)
        _render_table(ax, chunk, rank)
        rank += len(chunk)
        pdf.savefig(fig); plt.close(fig)


def _assign_levels(xs: List[float], min_dx: float) -> List[int]:
    """Give each label the lowest stacking level whose previous x is at least min_dx away, so
    labels closer than min_dx are separated vertically rather than horizontally."""
    last_x: List[float] = []
    levels = []
    for x in xs:
        level = next((k for k, lx in enumerate(last_x) if x - lx >= min_dx), len(last_x))
        if level == len(last_x):
            last_x.append(x)
        else:
            last_x[level] = x
        levels.append(level)
    return levels


def _gene_page(pdf, gene, positions, protein, labels):
    """KSS with the top positions labelled, then discriminative (up) and mutational (down)
    mirrored about a shared zero, laid out like Fig. 3 of the manuscript."""
    items = sorted((int(p), rec) for p, rec in positions.items())
    xs = [p for p, _ in items]
    span = (xs[-1] - xs[0]) if len(xs) > 1 else 1
    # Bar width follows the smallest gap between positions, so bars never overlap.
    gaps = [b - a for a, b in zip(xs, xs[1:])]
    width = max((min(gaps) if gaps else 1) * 0.9, span * 1e-4)

    # Labels stay on their bar's x; when two are too close, one moves up a level (never sideways).
    labelled = sorted((pos, positions[str(pos)].get("kss", 0.0), text)
                      for pos, text in labels.items())
    levels = _assign_levels([x for x, _, _ in labelled], span * 0.013) if labelled else []
    max_level = max(levels) if levels else 0
    step = 0.5
    top_y = 1.0 + (max_level + 1) * step + 0.1

    fig, (ax_k, ax_m) = plt.subplots(2, 1, figsize=(11, 5.2 + 0.7 * max_level), sharex=True,
                                     gridspec_kw={"height_ratios": [top_y, 1.3]})
    fig.suptitle(f"{gene}  -  scores per position", fontsize=13, fontweight="bold", y=0.99)

    ax_k.bar(xs, [rec.get("kss", 0.0) for _, rec in items], width=width,
             color=COLORS["kss"], edgecolor="white", linewidth=0.2)
    for (x, y, text), level in zip(labelled, levels):
        base = 1.0 + level * step
        ax_k.plot([x, x], [y, base], color="#999999", linewidth=0.5, zorder=1)
        ax_k.annotate(text, (x, base + 0.03), ha="center", va="bottom", fontsize=5.5,
                      rotation=90, zorder=5,
                      bbox=dict(boxstyle="square,pad=0.1", facecolor="white", edgecolor="none"))
    ax_k.set_ylim(0, top_y); ax_k.set_yticks([0, 0.5, 1.0]); ax_k.set_ylabel("KSS")
    ax_k.spines[["top", "right"]].set_visible(False)
    ax_k.grid(axis="y", color="#EEEEEE", linewidth=0.6); ax_k.set_axisbelow(True)

    ax_m.bar(xs, [rec.get("discriminative_score", 0.0) for _, rec in items], width=width,
             color=COLORS["discriminative"], edgecolor="white", linewidth=0.2,
             label="Discriminative")
    ax_m.bar(xs, [-rec.get("mutational_score", 0.0) for _, rec in items], width=width,
             color=COLORS["mutational"], edgecolor="white", linewidth=0.2, label="Mutational")
    ax_m.axhline(0, color="#888888", linewidth=0.6)
    ax_m.set_ylim(-1, 1); ax_m.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax_m.set_yticklabels(["1", "0.5", "0", "0.5", "1"])
    ax_m.set_ylabel("Mutational      Discriminative")
    ax_m.spines[["top", "right"]].set_visible(False)
    ax_m.grid(axis="y", color="#EEEEEE", linewidth=0.6); ax_m.set_axisbelow(True)
    ax_m.set_xlabel("genomic position (codon)")
    ax_m.set_xlim(xs[0] - width, xs[-1] + width)   # tight to the data, no wasted margins

    # Footer: protein name and protein characterization score, then the UniProt and GenBank links.
    p_score = next(iter(positions.values())).get("protein_score", 0.0)
    links = []
    if protein["accession"]:
        links.append(f"UniProt {protein['accession']}: {protein['link']}")
    if protein["genbank"]:
        links.append(f"GenBank reference {protein['genbank']}: {protein['genbank_link']}")
    fig.text(0.5, 0.035, f"{protein['name']}   |   protein characterization score {p_score:.3f}",
             ha="center", fontsize=9, color="#333333")
    fig.text(0.5, 0.012, "      ".join(links), ha="center", fontsize=8, color="#555555")
    fig.tight_layout(rect=[0, 0.055, 1, 0.96])
    pdf.savefig(fig); plt.close(fig)


def generate_report(results_dir: str, input_folder: str, dataset_name: str,
                    genes: List[str], weights: Dict[str, float],
                    output_path: str, top_n: int = 12) -> Optional[str]:
    """Write a PDF report of a KSS run. Returns the path, or None if no scores are found."""
    scores_by_gene = {g: s for g in genes if (s := _load_gene_scores(results_dir, g))}
    if not scores_by_gene:
        return None

    top = _top_positions(scores_by_gene, top_n)
    with PdfPages(output_path) as pdf:
        _table_pages(pdf, dataset_name, weights, scores_by_gene, top)
        for gene, positions in scores_by_gene.items():
            # Label the gene's top positions (those in the overall top-N) by their change.
            labels = {r["Position"]: r["variants"][0]["label_plain"] for r in top
                      if r["Gene"] == gene and r["variants"]}
            _gene_page(pdf, gene, positions, _protein_info(input_folder, gene), labels)
    return output_path
