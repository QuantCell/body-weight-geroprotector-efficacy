"""Regenerate manuscript main and extended-data figures from packaged data."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd


BLUE = "#2F7ED8"
ORANGE = "#F26B38"
GREY = "#8F8C85"
INK = "#1A1A1A"
GRID = "#E3E3E3"
SLATE = "#55534E"
PRIMARY_BLUE = "#0072B2"
SENSITIVITY_GREY = "#A7A9AC"
DEVELOPMENTAL_BLUE = "#1F5FA9"
DEVELOPMENTAL_ORANGE = "#D4692A"
FOUNDERS = ["BALB/cByJ", "C57BL/6J", "C3H/HeJ", "DBA/2J"]


def format_p(value: float) -> str:
    if value >= 0.0001:
        return f"P = {value:.4f}"
    coefficient, exponent = f"{value:.1e}".split("e")
    superscript = str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"P = {coefficient} × 10{superscript}"


def save_figure(fig, output: Path, stem: str, dpi: int = 300,
                target_height_pixels: int | None = None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    bbox = "tight"
    pad = None
    if target_height_pixels is not None:
        fig.canvas.draw()
        tight = fig.get_tightbbox(fig.canvas.get_renderer()).padded(0.1)
        bbox = Bbox.from_extents(tight.x0, tight.y0, tight.x1,
                                 tight.y0 + target_height_pixels / dpi)
        pad = 0
    fig.savefig(output / f"{stem}.png", dpi=dpi, bbox_inches=bbox,
                pad_inches=pad, facecolor="white")
    fig.savefig(output / f"{stem}.tif", dpi=600, bbox_inches=bbox,
                pad_inches=pad, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(output / f"{stem}.pdf", dpi=dpi, bbox_inches=bbox,
                pad_inches=pad, facecolor="white",
                metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def base_style(font_size: float = 8.7) -> None:
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.family": "DejaVu Sans", "font.size": font_size,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": True, "axes.spines.bottom": True,
        "axes.edgecolor": "#AAA79E", "axes.labelcolor": SLATE,
        "xtick.color": "#77746E", "ytick.color": "#77746E",
        "xtick.major.size": 0, "ytick.major.size": 0,
        "xtick.minor.size": 0, "ytick.minor.size": 0,
    })


QUARTER_NAMES = ["Lightest", "Second", "Third", "Heaviest"]
MEDIAN_COLUMN = {0: "standardized_control_median", 1: "standardized_treated_median"}


def panel_header(ax, letter: str, text: str) -> None:
    ax.annotate(letter, (0, 1), xycoords="axes fraction", xytext=(0, 7), textcoords="offset points",
                fontsize=11, fontweight="bold", va="bottom")
    ax.annotate(text, (0, 1), xycoords="axes fraction", xytext=(14, 8), textcoords="offset points",
                fontsize=8.6, color="#333333", va="bottom")


def figure1(package: Path, output: Path) -> None:
    """Standardized survival by matched pre-treatment quartile in male weight-lowering experiments."""
    frame = pd.read_csv(package / "data/derived/matched_quartile_standardized_survival.csv").set_index("quartile")
    stored = pd.read_csv(package / "data/derived/matched_quartile_survival_curves.csv")
    curves = {}
    for q in range(1, 5):
        for tx in (0, 1):
            z = stored[(stored.quartile == q) & (stored.Tx == tx)].sort_values("time")
            curves[q, tx] = np.r_[float(z.entry.iloc[0]), z.time], np.r_[1., z.survival]
            observed = float(z.loc[z.survival.le(.5), "time"].iloc[0])
            if observed != frame.loc[q, MEDIAN_COLUMN[tx]]:
                raise AssertionError(f"Curve median differs from summary: quartile {q}, treated {tx}")
    blues = ["#9DC3EA", "#5A95D6", DEVELOPMENTAL_BLUE, "#0B2E5E"]
    oranges = ["#F0A57A", "#DA7337", "#A94A14", "#5E2406"]

    def survival_axes(ax):
        ax.set_xlim(150, 1450)
        ax.set_xticks(range(200, 1401, 200))
        ax.set_ylim(0, 1.03)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "25", "50", "75", "100"])
        ax.set_xlabel("Age (days)")
        ax.grid(color="#ECEAE4", lw=0.7, zorder=0)
        ax.axhline(0.5, color="#B9B6AE", lw=0.8, ls=(0, (2, 2)), zorder=1)

    def headline(ax, label, value):
        ax.text(0.97, 0.90, label, transform=ax.transAxes, ha="right", fontsize=8.2, color=SLATE)
        ax.text(0.97, 0.82, value, transform=ax.transAxes, ha="right", fontsize=9.2, fontweight="bold", color=INK)

    base_style(8.7)
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 6.7), sharey=True)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.94, bottom=0.075, hspace=0.42, wspace=0.1)

    for ax, letter, q in ((axes[0, 0], "a", 1), (axes[0, 1], "b", 4)):
        for tx, colour, label in ((0, DEVELOPMENTAL_BLUE, "Control"), (1, DEVELOPMENTAL_ORANGE, "Treated")):
            ax.step(*curves[q, tx], where="post", color=colour, lw=1.8, label=label, zorder=3)
        mc, mt = frame.loc[q, MEDIAN_COLUMN[0]], frame.loc[q, MEDIAN_COLUMN[1]]
        for m, colour, side in ((mc, DEVELOPMENTAL_BLUE, -1), (mt, DEVELOPMENTAL_ORANGE, 1)):
            ax.plot([m, m], [0, 0.5], color=colour, lw=0.9, ls=(0, (2, 2)), zorder=2)
            ax.annotate(f"{m:.0f}", (m, 0.015), xytext=(3 * side, 0), textcoords="offset points",
                        ha="left" if side > 0 else "right", fontsize=8, color=SLATE)
        ax.annotate("", (mt, 0.5), (mc, 0.5), zorder=4,
                    arrowprops=dict(arrowstyle="->", color=INK, lw=1.2, shrinkA=0, shrinkB=0))
        headline(ax, "Median increase", f"+{mt - mc:.0f} days (+{frame.loc[q, 'gain_pct']:.1f}%)")
        survival_axes(ax)
        panel_header(ax, letter, f"{QUARTER_NAMES[q - 1]} quartile, control and treated")
    axes[0, 0].legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 0.04), fontsize=8.4)

    for ax, letter, tx, ramp, title in ((axes[1, 0], "c", 0, blues, "Controls, all four quartiles"),
                                        (axes[1, 1], "d", 1, oranges, "Treated mice, all four quartiles")):
        medians = frame[MEDIAN_COLUMN[tx]].to_numpy()
        for q in range(1, 5):
            ax.step(*curves[q, tx], where="post", color=ramp[q - 1], lw=1.7, label=QUARTER_NAMES[q - 1], zorder=3)
            ax.scatter([medians[q - 1]], [0.5], s=22, color=ramp[q - 1], edgecolors="white", linewidths=0.8, zorder=4)
        headline(ax, "Lightest to heaviest median", f"−{medians[0] - medians[-1]:.0f} days")
        survival_axes(ax)
        ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 0.02), fontsize=7.8, handlelength=1.6,
                  title="Pre-treatment quartile", title_fontsize=7.8, alignment="left")
        panel_header(ax, letter, title)
    for row in axes:
        row[0].set_ylabel("Surviving from six months (%)")
    save_figure(fig, output, "Figure_1_weight_stratified_survival")


def figure2(package: Path, output: Path) -> None:
    """Matched pre-treatment quartiles; weight labels show mean control weight."""
    plt.rcdefaults()
    data = package / "data" / "derived"
    frame = pd.read_csv(data / "matched_quartile_standardized_survival.csv")
    trend = pd.read_csv(data / "matched_quartile_trend.csv").iloc[0]
    sensitivity = pd.read_csv(data / "descriptive_sensitivities.csv")
    conventional = (sensitivity[
        sensitivity.analysis.eq("weight-lowering conventional quartiles")
    ].sort_values("bin").gain_pct.to_numpy())
    if conventional.shape != (4,):
        raise AssertionError("Expected four conventional-quartile sensitivity estimates")
    x = np.arange(4)
    labels = ["Lightest\nquartile", "Second\nquartile", "Third\nquartile", "Heaviest\nquartile"]

    base_style(8.7)
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 4.10), dpi=300,
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    ax.plot(x, frame.standardized_control_median, "o-", lw=2.8, ms=8,
            color=DEVELOPMENTAL_BLUE, label="Control")
    ax.plot(x, frame.standardized_treated_median, "o-", lw=2.8, ms=8,
            color=DEVELOPMENTAL_ORANGE, label="Treated")
    for i, row in frame.iterrows():
        ax.vlines(i, row.standardized_control_median, row.standardized_treated_median,
                  color="#B0B0B0", linewidth=1.0, zorder=1)
        first = i == 0
        ax.annotate(f"{row.mean_control_weight:.1f} g", (i, row.standardized_control_median),
                    xytext=(10 if first else 0, -25 if first else -19),
                    textcoords="offset points",
                    ha="left" if first else "center", color=SLATE, fontsize=8.2)
        ax.annotate(f"+{row.gain_pct:.1f}%",
                    (i, (row.standardized_control_median + row.standardized_treated_median) / 2),
                    xytext=(7, 0), textcoords="offset points", va="center",
                    color="#6A401E", fontsize=8.5, fontweight="bold")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Standardized median lifespan (days)")
    ax.set_ylim(620, 990)
    ax.set_yticks(range(650, 951, 50))
    ax.grid(axis="y", color="#DEDBD3", lw=0.8)
    ax.tick_params(axis="both", labelsize=8.6)
    ax.legend(frameon=False, fontsize=8.2)
    ax.set_title("a", loc="left", fontsize=11, fontweight="bold", pad=12)
    ax.text(0, 1.025, "Median lifespan across the four quartiles", transform=ax.transAxes,
            fontsize=8.5, color="#333333")

    ax = axes[1]
    width = 0.36
    ax.bar(x - width / 2, frame.gain_pct, width=width, color=PRIMARY_BLUE, label="Standardized")
    ax.bar(x + width / 2, conventional, width=width, color=SENSITIVITY_GREY, label="Pooled sensitivity")
    ax.set_xticks(x, ["Lowest", "Second", "Third", "Highest"])
    ax.tick_params(axis="both", labelsize=8.6)
    ax.set_ylabel("Median lifespan increase (%)")
    ax.set_ylim(0, 24)
    ax.grid(axis="y", color="#DEDBD3", lw=0.8)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2,
              fontsize=8.2, frameon=False, handlelength=1.5, columnspacing=1.3)
    ax.set_title("b", loc="left", fontsize=11, fontweight="bold", pad=12)
    ax.text(0, 1.025, "Standardized and site-pooled increases", transform=ax.transAxes,
            fontsize=8.5, color="#333333")
    ax.text(0.02, 0.96,
            f"Trend HR {trend.HR_per_quartile:.2f} ({trend.lo:.2f}–{trend.hi:.2f})\n{format_p(trend.p)}",
            transform=ax.transAxes, va="top", fontsize=8.2, fontweight="bold")
    fig.tight_layout(rect=[0, 0.12, 1, 0.93])
    save_figure(fig, output, "Figure_2_matched_quartiles")


def load_strain_means(path: Path, sex: str) -> dict[str, float]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    rows = raw.get("strainmeans", raw) if isinstance(raw, dict) else raw
    return {row["strain"]: float(row["mean"]) for row in rows
            if str(row.get("sex", "")).lower().startswith(sex) and row.get("mean") is not None}


def load_founder_lifespans(package: Path) -> pd.DataFrame:
    """Read sex-specific founder medians from Yuan et al. Table S1."""
    from zipfile import ZipFile
    import xml.etree.ElementTree as ET
    ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    with ZipFile(package / 'data/inputs/yuan2020_supplementary_tables.docx') as archive:
        document = ET.fromstring(archive.read('word/document.xml'))
    table = document.findall('.//w:tbl', ns)[0]
    rows = []
    for row in table.findall('w:tr', ns)[2:]:
        cells = [''.join(node.text or '' for node in cell.findall('.//w:t', ns)) for cell in row.findall('w:tc', ns)]
        if len(cells) < 6 or cells[0] not in FOUNDERS:
            continue
        for sex, column in [('F', 4), ('M', 5)]:
            rows.append({'Strain': cells[0], 'Sex': sex, 'median': float(cells[column])})
    frame = pd.DataFrame(rows)
    if len(frame) != 8 or frame.duplicated(['Strain', 'Sex']).any():
        raise ValueError('Expected one median for each founder strain and sex')
    return frame


from matplotlib.colors import to_rgb
from matplotlib.patches import Ellipse, PathPatch
from matplotlib.path import Path as MPath

def cubic(points, n: int = 40) -> np.ndarray:
    """Points along a cubic Bézier curve."""
    t = np.linspace(0, 1, n)[:, None]
    p = np.asarray(points, float)
    return (1 - t) ** 3 * p[0] + 3 * (1 - t) ** 2 * t * p[1] + 3 * (1 - t) * t ** 2 * p[2] + t ** 3 * p[3]

def closed_path(start, segments) -> MPath:
    """Closed path of cubic Bézier segments."""
    verts, codes = [tuple(start)], [MPath.MOVETO]
    for segment in segments:
        verts += [tuple(point) for point in segment]
        codes += [MPath.CURVE4] * 3
    return MPath(verts + [tuple(start)], codes + [MPath.CLOSEPOLY])

def draw_mouse(ax, width_in: float) -> None:
    """Side view of a laboratory mouse drawn with vector shapes."""
    lw = lambda points: max(points * width_in / 9.6, 0.1)
    pink, pink_edge = "#E0A497", "#C2867A"
    fur_edge, fur_line = "#4A423C", "#4E4640"

    tail = np.vstack([cubic([(3.9, 3.2), (2.2, 3.35), (0.1, 3.15), (0.55, 2.15)], 60),
                      cubic([(0.55, 2.15), (0.95, 1.35), (3.1, 1.45), (4.95, 1.05)], 60)[1:]])
    step = np.gradient(tail, axis=0)
    normal = np.c_[-step[:, 1], step[:, 0]] / np.hypot(step[:, 0], step[:, 1])[:, None]
    half = np.linspace(0.15, 0.018, len(tail))[:, None]
    ax.fill(*np.vstack([tail + normal * half, (tail - normal * half)[::-1]]).T, fc=pink, ec=pink_edge,
            lw=lw(0.5), zorder=1)
    ax.plot(*(tail + normal * half * 0.35)[:90].T, color="#EFC3B8", lw=lw(0.8), alpha=0.8, zorder=1,
            solid_capstyle="round")
    ax.add_patch(Ellipse((9.25, 5.98), 0.56, 0.66, angle=-12, fc="#CF9A90", ec=pink_edge, lw=lw(0.5), zorder=1.5))

    body = closed_path((3.85, 3.15), [
        ((3.55, 4.1), (4.0, 4.95), (4.8, 5.15)), ((5.5, 5.5), (6.4, 5.85), (7.3, 5.8)),
        ((8.2, 5.78), (8.9, 6.0), (9.5, 5.85)), ((10.2, 5.65), (10.85, 4.85), (11.05, 4.47)),
        ((11.1, 4.2), (10.8, 3.98), (10.4, 3.95)), ((9.85, 3.9), (9.25, 3.65), (8.95, 3.3)),
        ((8.72, 3.0), (8.62, 2.72), (8.4, 2.6)), ((7.6, 2.38), (6.6, 2.3), (5.9, 2.35)),
        ((5.0, 2.38), (4.45, 2.5), (4.12, 2.75)), ((3.9, 2.92), (3.8, 3.03), (3.85, 3.15))])
    outline = PathPatch(body, fc="none", ec="none", zorder=2)
    ax.add_patch(outline)
    stops = [(0, "#A39A90"), (0.14, "#877D73"), (0.35, "#6E645B"), (0.7, "#625850"), (1, "#554C45")]
    level = np.linspace(0, 1, 400)
    ramp = np.stack([np.interp(level, [s for s, _ in stops], [to_rgb(c)[k] for _, c in stops]) for k in range(3)], -1)
    image = ax.imshow(np.repeat(ramp[:, None, :], 400, axis=1), extent=(3.5, 11.2, 2.25, 6.1), origin="lower",
                      interpolation="bicubic", aspect="auto", zorder=2)
    image.set_clip_path(outline)
    for shape in (Ellipse((5.25, 3.75), 2.6, 2.45, angle=10, fc="#7A7067", ec="none", alpha=0.55, zorder=3),
                  Ellipse((5.25, 3.75), 2.6, 2.45, angle=10, fc="none", ec=fur_line, lw=lw(0.6), alpha=0.6, zorder=3),
                  Ellipse((10.0, 4.45), 1.7, 1.0, angle=-20, fc="#80766D", ec="none", alpha=0.4, zorder=3)):
        ax.add_patch(shape)
        shape.set_clip_path(outline)
    ax.plot(*cubic([(8.15, 4.6), (8.0, 4.0), (8.15, 3.35), (8.55, 2.75)]).T, color=fur_line, lw=lw(0.6),
            alpha=0.55, zorder=3, solid_capstyle="round")
    ax.add_patch(PathPatch(body, fc="none", ec=fur_edge, lw=lw(0.6), zorder=4))

    for heel, toe, t in (((5.35, 2.42), (6.95, 2.18), 0.17), ((8.45, 2.75), (9.3, 2.42), 0.12),
                         ((8.9, 2.95), (9.75, 2.62), 0.11)):
        heel, toe = np.asarray(heel), np.asarray(toe)
        a = toe - heel
        n = np.array([-a[1], a[0]]) / np.hypot(*a)
        foot = closed_path(heel + n * t, [(heel + a * 0.5 + n * t * 1.1, toe + n * t * 0.6, toe),
                                          (toe - n * t * 0.9, heel + a * 0.5 - n * t * 0.9, heel - n * t * 0.6),
                                          (heel - a * 0.12 - n * t * 0.3, heel - a * 0.12 + n * t * 0.7, heel + n * t)])
        ax.add_patch(PathPatch(foot, fc=pink, ec=pink_edge, lw=lw(0.5), zorder=5))
        for k in (0.72, 0.86):
            p = heel + a * k
            ax.plot([p[0], p[0] + n[0] * t * 0.5], [p[1] - t * 0.75, p[1] - t * 0.25], color=pink_edge, lw=lw(0.4),
                    zorder=6)

    ax.add_patch(Ellipse((8.72, 5.68), 0.78, 0.98, angle=-18, fc="#DCA79C", ec=pink_edge, lw=lw(0.6), zorder=6))
    ax.add_patch(Ellipse((8.76, 5.64), 0.52, 0.7, angle=-18, fc="#C98A80", ec="none", zorder=6))
    ax.add_patch(Ellipse((8.66, 5.7), 0.32, 0.52, angle=-18, fc="#DCA79C", ec="none", alpha=0.55, zorder=6))
    ax.add_patch(Ellipse((9.92, 4.95), 0.36, 0.33, fc="#120C09", zorder=7))
    ax.add_patch(Ellipse((9.98, 5.02), 0.09, 0.08, fc="white", alpha=0.85, zorder=8))
    ax.add_patch(Ellipse((11.0, 4.42), 0.26, 0.22, angle=-20, fc="#D98C80", ec="#B9736A", lw=lw(0.3), zorder=7))
    ax.plot(*cubic([(10.55, 4.02), (10.65, 3.97), (10.75, 3.98), (10.85, 4.05)], 20).T, color=fur_edge, lw=lw(0.4),
            zorder=7)
    for end, bend in (((12.1, 5.05), 0.25), ((12.25, 4.6), 0.12), ((12.15, 4.1), -0.05), ((11.85, 3.62), -0.2)):
        start, end = np.array((10.75, 4.3)), np.array(end)
        middle = (start + end) / 2 + np.array([0, bend])
        ax.plot(*cubic([start, start + (middle - start) * 0.7, middle + (end - middle) * 0.4, end], 30).T,
                color="#C9C2BA", lw=lw(0.35), zorder=7, solid_capstyle="round")
    ax.set_xlim(0, 12.4)
    ax.set_ylim(0.8, 6.8)
    ax.set_aspect("equal")
    ax.axis("off")


def figure3(package: Path, output: Path) -> None:
    """Body-weight, lifespan and adiposity strips; vertical positions are jittered for visibility."""
    inputs = package / "data" / "inputs"
    weight = load_strain_means(inputs / "mpd_32603_strainmeans.json", "m")
    age = load_strain_means(inputs / "mpd_32680_strainmeans.json", "m")
    matched = {strain: value for strain, value in weight.items()
               if strain in age and age[strain] <= 36.0}
    if len(matched) != 27 or not all(founder in matched for founder in FOUNDERS):
        raise AssertionError(f"Expected 27 age-matched strains including all founders, found {len(matched)}")
    survey_mean = float(np.mean(list(matched.values())))
    if round(survey_mean, 1) != 30.3:
        raise AssertionError(f"Unexpected survey mean: {survey_mean}")

    quartiles = pd.read_csv(package / "data" / "derived" /
                            "matched_quartile_standardized_survival.csv")
    q = quartiles.sort_values('quartile')
    summary = pd.read_csv(package / "data/derived/founder_strain_weight_comparison.csv").set_index("quantity").value
    control_mean = float(summary.loc["itp_male_control_6m_mean_weight_g"])
    quartile_marks = [("ITP lightest quartile", float(q.iloc[0].mean_control_weight)),
                      ("ITP overall mean", control_mean),
                      ("ITP heaviest quartile", float(q.iloc[-1].mean_control_weight))]
    body = pd.read_excel(package / "supplementary" / "Supplementary_Data_1-10.xlsx",
                         sheet_name="SD5_body_composition")
    rapa = body[(body["Intervention"] == "Rapamycin") &
                (body["strain"] == "UM-HET3") &
                body["dosage"].astype(str).str.contains("42 ppm", regex=False)]
    adiposity = rapa.set_index("sex")["bodyfat_pct_control"].astype(float)
    if not {"Male", "Female"}.issubset(adiposity.index):
        raise AssertionError("Missing UM-HET3 control adiposity values")

    panels = [
        ("a", "Body weight, ♂", matched, "Body weight (g)",
         quartile_marks,
         survey_mean),
        ("d", "Body fat, ♂", load_strain_means(inputs / "mpd_10331_strainmeans.json", "m"),
         "Body fat (%)", [("ITP UM-HET3", adiposity.loc["Male"])], None),
        ("e", "Body fat, ♀", load_strain_means(inputs / "mpd_10331_strainmeans.json", "f"),
         "Body fat (%)", [("ITP UM-HET3", adiposity.loc["Female"])], None),
    ]

    base_style(10)
    plt.rcParams["axes.spines.left"] = False
    fig = plt.figure(figsize=(8.8, 7.7))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.0], hspace=.54, wspace=.10)
    axes = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])]
    survival_axes = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    for ax, (letter, title, values_by_strain, x_label, marks, reference_mean) in zip(axes, panels):
        values = np.array(list(values_by_strain.values()))
        names = list(values_by_strain.keys())
        y = np.random.default_rng(0).uniform(-0.17, 0.17, len(values))
        founder_mask = np.array([name in FOUNDERS for name in names])
        ax.scatter(values[~founder_mask], y[~founder_mask], s=26, c=GREY, alpha=0.55,
                   edgecolors="none", label=f"other strains ({(~founder_mask).sum()})", zorder=2)
        ax.scatter(values[founder_mask], y[founder_mask], s=64, c=BLUE,
                   edgecolors="white", linewidths=0.8, label="UM-HET3 founders (4)", zorder=4)
        for name, value, yy in sorted(zip(names, values, y), key=lambda item: item[1]):
            if name in FOUNDERS:
                below = FOUNDERS.index(name) % 2 == 0
                ax.annotate(name, (value, yy), xytext=(0, -17 if below else 13),
                            textcoords="offset points", ha="center", fontsize=7.2,
                            color="#3A5F8A", zorder=6)
        if reference_mean is not None:
            ax.plot([reference_mean, reference_mean], [-0.253, 0.690], ls=(0, (4, 3)),
                    color=SLATE, lw=1.6, zorder=5)
            ax.annotate(f"Survey mean\n{reference_mean:.1f} g", (reference_mean, 0.80 * 1.55 - 0.55), xytext=(0, 2), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7.6, color=SLATE,
                        fontweight="bold", zorder=7)
        for index, (label, value) in enumerate(marks):
            high = 0.80
            ax.axvline(value, color=ORANGE, lw=2.0, ymin=0.34, ymax=high, zorder=5)
            ax.annotate(f"{label}\n{value:.1f}", (value, high * 1.55 - 0.55), xytext=(0, 2),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=7.6, color="#A8471C", fontweight="bold", zorder=7)
        ax.set_ylim(-0.95 if reference_mean is not None else -0.72, 1.10)
        ax.set_yticks([])
        ax.set_xlabel(x_label)
        ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", fontsize=9.3,
                     color="#1A1A1A", pad=11)
        ax.grid(axis="x", color="#E4E1DA", lw=0.7, zorder=0)
        if len(marks) == 1:
            note = f"{int((values >= marks[0][1]).sum())} of {len(values)} strains at or above the ITP value"
        else:
            heaviest_count = int((values >= marks[-1][1]).sum())
            heaviest_text = "none" if heaviest_count == 0 else str(heaviest_count)
            note = (f"{int((values >= marks[1][1]).sum())} of {len(values)} strains at or above the ITP value; "
                    f"{int((values >= marks[0][1]).sum())} reach the lightest ITP quartile, "
                    f"{heaviest_text} the heaviest.")
        ax.annotate(note, (0.995, 0.02), xycoords="axes fraction", ha="right",
                    fontsize=8.4, color=SLATE)
    controls = pd.read_excel(package / "supplementary/Supplementary_Data_1-10.xlsx", sheet_name="SD10_control_lifespan")
    induced = controls.lifespan_evidence_quote.fillna("").str.contains("doxorubicin", case=False)
    controls = controls[~induced & controls.sex.isin(["Female", "Male"])].copy()
    controls["origin"] = np.where(controls.itp_status.astype(str).str.upper().eq("ITP"), "ITP", "Non-ITP")
    controls = controls.drop_duplicates(["sex", "origin", "pmid", "lifespan_control_days"])
    reference = pd.read_csv(package / "data/derived/control_lifespan_ceiling.csv")
    reference = reference[reference.specification.eq("unique_control_groups")].set_index("sex")
    for ax, sex, letter in zip(survival_axes, ["Male", "Female"], ["b", "c"]):
        selected=controls[controls.sex.eq(sex)&controls.origin.eq('Non-ITP')].sort_values(['lifespan_control_days','pmid'])
        values=selected.lifespan_control_days.to_numpy()
        median=float(reference.loc[sex,'itp_median_days'])
        assert len(values)==int(reference.loc[sex,'non_itp_n'])
        y=np.random.default_rng(0).uniform(-.17,.17,len(values))
        ax.scatter(values,y,s=26,c=GREY,alpha=.55,edgecolors='none',zorder=2)
        high=.80
        ax.axvline(median,color=ORANGE,lw=2.0,ymin=.34,ymax=high,zorder=5)
        ax.annotate(f'ITP median\n{int(np.floor(median + 0.5))} days',(median,high*1.55-.55),xytext=(0,2),textcoords='offset points',ha='center',va='bottom',fontsize=7.6,color='#A8471C',fontweight='bold',zorder=7)
        founder_data=load_founder_lifespans(package)
        founder_data=founder_data[founder_data.Sex.eq('M' if sex=='Male' else 'F')].set_index('Strain')
        config={'Male':[('BALB/cByJ',.37,0,14),('C57BL/6J',.37,0,14),('DBA/2J',-.36,-25,-14),('C3H/HeJ',-.36,30,-14)],'Female':[('DBA/2J',.37,0,14),('BALB/cByJ',.37,0,14),('C3H/HeJ',-.36,-4,-14),('C57BL/6J',-.36,0,-14)]}
        for strain,yy,dx,dy in config[sex]:
            value=float(founder_data.loc[strain,'median'])
            ax.scatter(value,yy,s=64,c=BLUE,edgecolors='white',linewidths=.8,zorder=6)
            ax.annotate(strain,(value,yy),xytext=(dx,dy),textcoords='offset points',ha='center',va='bottom' if dy>0 else 'top',fontsize=7.2,color='#3A5F8A',zorder=7)
        ax.set_ylim(-.95,1.10);ax.set_yticks([])
        ax.set_xlim(400,1050);ax.set_xticks([400,600,800,1000])
        ax.set_xlabel('Control-group lifespan (days)')
        ax.set_title(f"{letter}  Control survival, {'♂' if sex=='Male' else '♀'}",loc='left',fontweight='bold',fontsize=9.3,color='#1A1A1A',pad=11)
        ax.grid(axis='x',color='#E4E1DA',lw=.7,zorder=0)
        note=f'{int((values>=median).sum())} of {len(values)} non-ITP groups reach the ITP median'
        ax.annotate(note,(.995,.02),xycoords='axes fraction',ha='right',fontsize=8.4,color=SLATE)
        ax.spines['left'].set_visible(False)
    from matplotlib.lines import Line2D
    handles=[Line2D([],[],marker='o',linestyle='none',markerfacecolor=GREY,markeredgecolor='none',markersize=5,label='Laboratory strains (a, d, e) or non-ITP control groups (b, c)'),Line2D([],[],marker='o',linestyle='none',markerfacecolor=BLUE,markeredgecolor='none',markersize=5,label='UM-HET3 founder strains'),Line2D([],[],color=ORANGE,lw=2,label='ITP reference values')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.58,1.005),ncol=3,frameon=False,fontsize=6.7,handletextpad=.5,columnspacing=1.0)
    fig.subplots_adjust(left=.10, right=.97, bottom=.07, top=.93)
    mouse = fig.add_axes([.095, .965, .090, .060])
    draw_mouse(mouse, .09 * fig.get_figwidth())
    output.mkdir(parents=True, exist_ok=True)
    save_figure(fig, output, 'Figure_3_strain_context', dpi=300)

def figure5(package: Path, output: Path) -> None:
    """Early-life and adult-weight mortality estimates."""
    data = package / "data" / "derived"
    age = pd.read_csv(data / "developmental_figure_estimates.csv")
    control = pd.read_csv(data / "control_weight_gradient_by_age.csv")

    early_age = age[age.panel.eq("age")].set_index("sex")
    control_6_12 = control[control.weight_age.isin([6, 12])]
    gradients = {
        ("m", 1.4): early_age.loc["m", ["HR", "lo", "hi"]].astype(float).to_numpy(),
        ("f", 1.4): early_age.loc["f", ["HR", "lo", "hi"]].astype(float).to_numpy(),
    }
    for row in control_6_12.itertuples(index=False):
        gradients[(row.sex, float(row.weight_age))] = np.array([row.HR_per_SD, row.lo, row.hi], float)

    decomp = age[age.panel.eq("decomposition")].copy()
    expected_labels = ["1.4-month alone", "6-month alone", "1.4-month adjusted", "6-month adjusted"]
    if decomp.label.tolist() != expected_labels:
        raise AssertionError(f"Unexpected decomposition rows: {decomp.label.tolist()}")

    base_style(8.5)
    fig = plt.figure(figsize=(7.09, 3.10))
    grid = fig.add_gridspec(1, 2, wspace=0.44, left=0.075, right=0.985, top=0.85, bottom=0.235)

    ax = fig.add_subplot(grid[0, 0])
    ages = [1.4, 6, 12]
    male = np.array([gradients[("m", value)] for value in ages])
    female = np.array([gradients[("f", value)] for value in ages])
    ax.fill_between(ages, male[:, 1], male[:, 2], color=DEVELOPMENTAL_BLUE, alpha=0.15, lw=0)
    ax.fill_between(ages, female[:, 1], female[:, 2], color=DEVELOPMENTAL_ORANGE, alpha=0.15, lw=0)
    ax.plot(ages, male[:, 0], "-o", color=DEVELOPMENTAL_BLUE, lw=1.7, ms=4.4,
            mfc=DEVELOPMENTAL_BLUE, mec="white", mew=1)
    ax.plot(ages, female[:, 0], "-^", color=DEVELOPMENTAL_ORANGE, lw=1.7, ms=4.4,
            mfc=DEVELOPMENTAL_ORANGE, mec="white", mew=1)
    ax.axhline(1, color="#999999", lw=0.9)
    ax.set_xticks(ages, ["1.4", "6", "12"])
    ax.set_xlim(0.2, 13.2)
    ax.set_ylim(0.86, 1.50)
    ax.set_xlabel("Age at weighing (months)", fontsize=8.8, color=INK)
    ax.set_ylabel("Mortality hazard per +1 s.d.", fontsize=8.8, color=INK)
    ax.annotate("males", (12, male[-1, 0]), xytext=(-4, 10), textcoords="offset points",
                color=DEVELOPMENTAL_BLUE, fontsize=7.8, fontweight="bold", ha="right")
    ax.annotate("females", (12, female[-1, 0]), xytext=(-4, -15), textcoords="offset points",
                color=DEVELOPMENTAL_ORANGE, fontsize=7.8, fontweight="bold", ha="right")
    ax.set_title("a", loc="left", fontsize=10.5, fontweight="bold", pad=11)
    ax.text(0, 1.035, "Weight predicts mortality", transform=ax.transAxes, fontsize=8.2, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.6)

    ax = fig.add_subplot(grid[0, 1])
    labels = ["1.4-mo\nalone", "6-mo\nalone", "1.4-mo\nadj.", "6-mo\nadj."]
    colours = [GREY, DEVELOPMENTAL_BLUE, GREY, DEVELOPMENTAL_BLUE]
    x = np.arange(4)
    for index, row in enumerate(decomp.itertuples(index=False)):
        ax.plot([x[index], x[index]], [row.lo, row.hi], color=colours[index], lw=1.6)
        ax.plot(x[index], row.HR, "o", color=colours[index], ms=6.4,
                mec="white", mew=1.2, zorder=3)
    ax.axhline(1, color="#999999", lw=0.9)
    ax.axvline(1.5, color="#CCCCCC", lw=0.9, ls=":")
    ax.set_xticks(x, labels, fontsize=6.8)
    ax.set_ylabel("Mortality hazard per +1 s.d.", fontsize=8.8, color=INK)
    ax.set_ylim(0.82, 1.78)
    ax.text(0.5, 1.70, "separate", ha="center", fontsize=7.6, color=INK)
    ax.text(2.5, 1.70, "adjusted", ha="center", fontsize=7.6, color=INK)
    ax.set_title("b", loc="left", fontsize=10.5, fontweight="bold", pad=11)
    ax.text(0, 1.035, "Adult, not early-life, weight", transform=ax.transAxes,
            fontsize=8.2, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.6)

    for axis in fig.axes:
        axis.tick_params(axis="both", which="major", labelsize=7.8,
                         colors=INK, labelcolor=INK)
        axis.xaxis.label.set_color(INK)
        axis.yaxis.label.set_color(INK)
        axis.set_axisbelow(True)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
    save_figure(fig, output, "Figure_5_developmental", dpi=600)


def figure4(package: Path, output: Path) -> dict:
    """Median-lifespan increase in the lightest and heaviest pre-treatment quartile, by compound."""
    def signed_percent(value):
        return f"{value:+.1f}%".replace("-", "−")

    def p_text(value):
        if value < 0.001:
            mantissa, exponent = f"{value:.1e}".split("e")
            superscript = str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
            return f"P = {mantissa} × 10{superscript}"
        return f"P = {value:.2g}"

    compound_names = {"17alpha-estradiol": "17α-estradiol", "2-(2-hydroxyphenyl)benzoxazole": "HBX"}
    base_style()
    table = pd.read_csv(package / "data/derived/compound_descriptive.csv")
    table["compound"] = table.compound.str.replace("2-(2-Hydroxyphenyl)", "2-(2-hydroxyphenyl)", regex=False)
    common = pd.read_csv(package / "data/derived/compound_common_interaction.csv").iloc[0]
    is_pooled = table.compound.str.contains("pool", case=False)
    pooled = table[is_pooled].iloc[0]
    compounds = table[~is_pooled].assign(midpoint=lambda d: (d.lowest_gain_pct + d.highest_gain_pct) / 2)
    compounds = compounds.sort_values("midpoint")

    def dumbbell(ax, y, low, high, line, low_edge, high_fill, size, width):
        ax.plot([low, high], [y, y], color=line, lw=width, zorder=2)
        ax.scatter([low], [y], color="white", edgecolors=low_edge, s=size, lw=1.3, zorder=3)
        ax.scatter([high], [y], color=high_fill, s=size, zorder=3)
        (left, left_colour), (right, right_colour) = sorted([(low, low_edge), (high, high_fill)])
        for x, colour, align, shift in ((left, left_colour, "right", -1.1), (right, right_colour, "left", 1.1)):
            ax.text(x + shift, y, signed_percent(x), ha=align, va="center", fontsize=7.2, color=colour,
                    zorder=4, bbox=dict(facecolor="white", edgecolor="none", pad=0.3))

    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    rows = list(range(len(compounds)))
    for y, (_, row) in zip(rows, compounds.iterrows()):
        dumbbell(ax, y, row.lowest_gain_pct, row.highest_gain_pct, BLUE + "73", GREY, BLUE, 34, 1.6)
    pooled_row = len(compounds) + 1.1
    dumbbell(ax, pooled_row, pooled.lowest_gain_pct, pooled.highest_gain_pct, INK, INK, INK, 44, 2.0)
    ax.text((pooled.lowest_gain_pct + pooled.highest_gain_pct) / 2, pooled_row - 0.62,
            f"Common treatment × weight HR {common.estimate:.2f} ({common.lower_95:.2f}–{common.upper_95:.2f}); "
            f"{p_text(common.p)}", ha="center", va="center", fontsize=7.2, color=SLATE)
    labels = [f"{compound_names.get(name, name)} ({int(n)}; {int(mice):,})"
              for name, n, mice in zip(compounds.compound, compounds.experiments, compounds.unique_mice)]
    labels.append(f"All 14 compounds ({int(pooled.experiments)}; {int(pooled.unique_mice):,})")
    ax.set_yticks(rows + [pooled_row], labels)
    ax.axvline(0, color="#AAA79E", lw=1, zorder=1)
    ax.set_xlim(-17, 39)
    ax.set_xticks([-10, 0, 10, 20, 30])
    ax.set_xlabel("Median lifespan increase versus controls (%)")
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color="#E4E1DA", lw=0.7, zorder=0)
    ax.scatter([], [], color="white", edgecolors=GREY, s=34, lw=1.3, label="Lightest quartile")
    ax.scatter([], [], color=BLUE, s=34, label="Heaviest quartile")
    ax.legend(frameon=False, loc="lower right", fontsize=7.8)
    fig.tight_layout()
    save_figure(fig, output, "Figure_4_compound_weight_strata")
    names = {"lowest_control_median": "lightest_control_median_days", "lowest_treated_median": "lightest_treated_median_days",
             "lowest_gain_pct": "lightest_increase_percent", "highest_control_median": "heaviest_control_median_days",
             "highest_treated_median": "heaviest_treated_median_days", "highest_gain_pct": "heaviest_increase_percent"}
    panel = table[["compound", "experiments", "unique_mice", *names]].rename(columns=names)
    return {"Fig4": panel.round({"lightest_increase_percent": 2, "heaviest_increase_percent": 2}),
            "Fig4_interaction": pd.DataFrame({"quantity": ["common treatment × weight HR", "lower 95% CI", "upper 95% CI", "P"],
                                              "value": [common.estimate, common.lower_95, common.upper_95, common.p]})}


def extended_figure1(package: Path, output: Path) -> None:
    """All-experiment treatment-by-weight estimates."""
    plt.rcdefaults()
    data = package / "data"
    interactions = pd.read_csv(data / "derived" / "experiment_interactions.csv")
    weight = pd.read_csv(data / "derived" / "experiment_weight_effects.csv")
    post = weight[weight.post_dosing.astype(str).str.lower().eq("true")].copy()
    post["weighted_diff"] = post.diff_g * post.n_tx
    lowering = post.groupby(["cohort", "group"], as_index=False).agg(
        weighted_diff=("weighted_diff", "sum"), n_tx=("n_tx", "sum")
    )
    lowering["weight_lowering"] = lowering.weighted_diff / lowering.n_tx < -1
    interactions = interactions.merge(
        lowering[["cohort", "group", "weight_lowering"]],
        on=["cohort", "group"], how="left", validate="many_to_one"
    )
    interactions["weight_lowering"] = interactions.weight_lowering.eq(True)
    observed = interactions.groupby("sex").weight_lowering.sum().to_dict()
    if observed != {"m": 42, "f": 40}:
        raise AssertionError(f"Unexpected weight-lowering counts: {observed}")
    pooled = pd.read_csv(data / "derived" / "primary_moderation.csv").set_index("sex")

    base_style(10)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 6.2), sharex=True)
    panels = zip(axes, ["m", "f"], ["a", "b"],
                 [f"Males ({int(pooled.loc['m', 'experiments'])} experiments)",
                  f"Females ({int(pooled.loc['f', 'experiments'])} experiments)"])
    for ax, sex, letter, title in panels:
        frame = interactions[interactions.sex == sex].sort_values("int_hr").reset_index(drop=True)
        y = np.arange(len(frame))
        colours = np.where(frame.weight_lowering, BLUE, GREY)
        for yi, lo, hi, colour in zip(y, frame.int_lo, frame.int_hi, colours):
            ax.plot([max(lo, 0.35), min(hi, 2.1)], [yi, yi], color=colour, alpha=0.25, lw=0.8)
        ax.scatter(frame.int_hr.clip(0.35, 2.1), y, c=colours, s=12, alpha=0.8, edgecolors="none")
        ax.axvline(1, color="#AAA79E", lw=1.2)
        estimate = pooled.loc[sex]
        ax.errorbar(estimate.HR, -7,
                    xerr=[[estimate.HR - estimate.lo], [estimate.hi - estimate.HR]],
                    fmt="D", ms=8, color="black", capsize=3, lw=2)
        ax.text(0.37, -11, f"Pooled HR {estimate.HR:.2f} ({estimate.lo:.2f}–{estimate.hi:.2f})",
                fontsize=9, fontweight="bold")
        ax.set_ylim(len(frame) + 2, -15)
        ax.set_yticks([])
        ax.set_title(letter, loc="left", fontweight="bold", fontsize=13, pad=22)
        ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 5),
                    textcoords="offset points", ha="left", va="bottom",
                    fontsize=11, color="#333333")
        ax.set_xlabel("Treatment × weight interaction HR per +1 s.d.")
        ax.grid(axis="x", color="#E4E1DA", lw=0.7)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker='o', linestyle='none', color=BLUE, label='Weight-lowering experiments'),
                        Line2D([], [], marker='o', linestyle='none', color=GREY, label='Non-weight-lowering experiments'),
                        Line2D([], [], marker='D', linestyle='none', color='black', label='Pooled estimate')],
               loc='lower center', bbox_to_anchor=(.5, -.015), ncol=3, frameon=False, fontsize=9)
    fig.tight_layout(rect=[0, 0.07, 1, 0.98])
    save_figure(fig, output, "Extended_Data_Figure_1_all_experiment_interactions")



def figure_source_tables(package: Path) -> dict:
    """Panel-labelled numerical data used by the manuscript figures."""
    data = package / "data/derived"
    read = lambda name: pd.read_csv(data / (name + ".csv"))
    matched = read("matched_quartile_standardized_survival")
    curves = read("matched_quartile_survival_curves").rename(columns={"Tx": "treated", "entry": "entry_days", "time": "age_days", "survival": "survival_fraction"})
    curves["group"] = curves.treated.map({0: "Control", 1: "Treated"})
    sensitivities = read("descriptive_sensitivities")
    sensitivities = sensitivities[sensitivities.analysis.eq("weight-lowering conventional quartiles")].sort_values("bin")
    inputs = package / "data/inputs"
    weights = load_strain_means(inputs / "mpd_32603_strainmeans.json", "m")
    ages = load_strain_means(inputs / "mpd_32680_strainmeans.json", "m")
    strains = pd.DataFrame([{"strain": s, "sex": "Male", "mean_weight_g": v, "age_weeks": ages[s], "founder": s in FOUNDERS}
                           for s, v in weights.items() if s in ages and ages[s] <= 36.0])
    fat = pd.DataFrame([{"strain": s, "sex": sex, "body_fat_pct": v, "founder": s in FOUNDERS}
                        for sex, key in [("Male", "m"), ("Female", "f")]
                        for s, v in load_strain_means(inputs / "mpd_10331_strainmeans.json", key).items()])
    workbook = package / "supplementary/Supplementary_Data_1-10.xlsx"
    controls = pd.read_excel(workbook, sheet_name="SD10_control_lifespan")
    controls = controls[~controls.lifespan_evidence_quote.fillna("").str.contains("doxorubicin", case=False) & controls.sex.isin(["Female", "Male"])].copy()
    controls["origin"] = np.where(controls.itp_status.astype(str).str.upper().eq("ITP"), "ITP", "Non-ITP")
    controls = controls.drop_duplicates(["sex", "origin", "pmid", "lifespan_control_days"])
    controls = controls[["sex", "origin", "first_author", "year", "pmid", "strain", "lifespan_measure", "lifespan_control_days"]].sort_values(["sex", "origin", "lifespan_control_days", "pmid"]).reset_index(drop=True)
    body = pd.read_excel(workbook, sheet_name="SD5_body_composition")
    body = body[body.Intervention.eq("Rapamycin") & body.strain.eq("UM-HET3") & body.dosage.astype(str).str.contains("42 ppm", regex=False)]
    ref = read("control_lifespan_ceiling");ref = ref[ref.specification.eq("unique_control_groups")].reset_index(drop=True)
    references = [{"panel": "a", "reference": "Laboratory-strain mean", "sex": "Male", "value": strains.mean_weight_g.mean(), "unit": "g"},
                  {"panel": "a", "reference": "ITP overall mean", "sex": "Male", "value": float(read("founder_strain_weight_comparison").set_index("quantity").value.loc["itp_male_control_6m_mean_weight_g"]), "unit": "g"}]
    references += [{"panel": "a", "reference": f"ITP {name} quartile", "sex": "Male", "value": float(matched.set_index("quartile").loc[q, "mean_control_weight"]), "unit": "g"} for q, name in [(1, "lightest"), (4, "heaviest")]]
    references += [{"panel": "b" if row.sex == "Male" else "c", "reference": "ITP median across control groups", "sex": row.sex, "value": row.itp_median_days, "unit": "days"} for row in ref.itertuples()]
    references += [{"panel": "d" if row.sex == "Male" else "e", "reference": "ITP body fat", "sex": row.sex, "value": float(row.bodyfat_pct_control), "unit": "%"} for row in body.itertuples()]
    compounds = read("compound_descriptive")
    compounds["sort_mean_pct"] = (compounds.lowest_gain_pct + compounds.highest_gain_pct) / 2
    compounds = pd.concat([compounds[compounds.compound.eq("Pooled single-agent set")], compounds[~compounds.compound.eq("Pooled single-agent set")].sort_values("sort_mean_pct", ascending=False)], ignore_index=True)
    names = {"17alpha-estradiol": "17α-estradiol", "2-(2-Hydroxyphenyl)benzoxazole": "HBX", "Pooled single-agent set": "All 14 compounds"}
    compounds["display_label"] = [f"{names.get(r.compound, r.compound)} ({int(r.experiments)}; {int(r.unique_mice):,})" for r in compounds.itertuples()]
    developmental = read("developmental_figure_estimates")
    gradients = read("control_weight_gradient_by_age");gradients = gradients[gradients.weight_age.isin([6, 12])].rename(columns={"weight_age": "weight_age_months", "HR_per_SD": "HR"})
    early = developmental[developmental.panel.eq("age")].copy();early["weight_age_months"] = 1.4
    gradients = pd.concat([early[["sex", "weight_age_months", "n", "HR", "lo", "hi", "p"]], gradients[["sex", "weight_age_months", "n", "HR", "lo", "hi", "p"]]], ignore_index=True)
    interactions = read("experiment_interactions")
    effects = read("experiment_weight_effects");effects = effects[effects.post_dosing.astype(str).str.lower().eq("true")].copy()
    effects["weighted_diff"] = effects.diff_g * effects.n_tx
    effects = effects.groupby(["cohort", "group"], as_index=False).agg(weighted_diff=("weighted_diff", "sum"), n_tx=("n_tx", "sum"))
    effects["weight_lowering"] = effects.weighted_diff / effects.n_tx < -1
    interactions = interactions.merge(effects[["cohort", "group", "weight_lowering"]], on=["cohort", "group"], how="left", validate="many_to_one")
    interactions["weight_lowering"] = interactions.weight_lowering.eq(True)
    interactions = interactions.sort_values(["sex", "int_hr"]).reset_index(drop=True)
    matched = matched.rename(columns={"mean_control_weight": "mean_control_weight_g", "standardized_control_median": "control_median_days", "standardized_treated_median": "treated_median_days", "gain_days": "median_increase_days", "gain_pct": "median_increase_percent"})
    curves["surviving_percent"] = curves.survival_fraction * 100
    common_names = {"lo": "lower_95", "hi": "upper_95", "p": "P"}
    curves = curves.drop(columns=["treated"])
    return {
        "Figure_1": {"Fig1a": curves[curves.quartile.eq(1)].reset_index(drop=True), "Fig1b": curves[curves.quartile.eq(4)].reset_index(drop=True), "Fig1c": curves[curves.group.eq("Control")].reset_index(drop=True), "Fig1d": curves[curves.group.eq("Treated")].reset_index(drop=True), "Fig1_medians": matched},
        "Figure_2": {"Fig2a": matched, "Fig2b": sensitivities.rename(columns={"mean_control_weight": "mean_control_weight_g", "control_median": "control_median_days", "treated_median": "treated_median_days", "gain_pct": "median_increase_percent"}).reset_index(drop=True), "Fig2b_trend": read("matched_quartile_trend").rename(columns=common_names)},
        "Figure_3": {"Fig3a": strains, "Fig3bc_controls": controls, "Fig3bc_founders": load_founder_lifespans(package).rename(columns={"median": "median_lifespan_days"}), "Fig3de_body_fat": fat, "Fig3_reference_lines": pd.DataFrame(references)},
        "Figure_4": {"Fig4": compounds.rename(columns={"lowest_control_median": "lightest_control_median_days", "lowest_treated_median": "lightest_treated_median_days", "lowest_gain_pct": "lightest_increase_percent", "highest_control_median": "heaviest_control_median_days", "highest_treated_median": "heaviest_treated_median_days", "highest_gain_pct": "heaviest_increase_percent", "overall_control_median": "overall_control_median_days", "overall_treated_median": "overall_treated_median_days", "overall_gain_pct": "overall_increase_percent"}), "Fig4_interaction": read("compound_common_interaction").rename(columns=common_names)},
        "Figure_5": {"Fig5a": gradients.rename(columns=common_names), "Fig5b": developmental[developmental.panel.eq("decomposition")].reset_index(drop=True).rename(columns=common_names)},
        "Extended_Data_Figure_1": {"EDFig1ab": interactions.rename(columns={"group": "intervention_group", "weight_age": "weight_age_months", "int_hr": "interaction_HR_per_SD", "int_lo": "lower_95", "int_hi": "upper_95", "int_p": "P"}), "EDFig1ab_pooled": read("primary_moderation").rename(columns=common_names)},
    }


def verify_figure_source_data(package: Path, source_directory: Path | None = None) -> None:
    """Compare the released figure source data with the current analysis results."""
    for figure, sheets in figure_source_tables(package).items():
        path = (source_directory or package / "figures/source_data") / ("Source_Data_" + figure + ".xlsx")
        for sheet, expected in sheets.items():
            actual = pd.read_excel(path, sheet_name=sheet)
            expected = expected.reset_index(drop=True)
            if actual.shape != expected.shape or list(actual.columns) != list(expected.columns):
                raise AssertionError(f"Figure source-data structure differs: {figure}, {sheet}")
            for column in expected:
                a, b = actual[column], expected[column]
                if not np.array_equal(a.isna(), b.isna()):
                    raise AssertionError(f"Missing values differ: {figure}, {sheet}, {column}")
                use = b.notna()
                if pd.api.types.is_numeric_dtype(b) and not pd.api.types.is_bool_dtype(b):
                    equal = np.allclose(pd.to_numeric(a[use]), b[use].astype(float), rtol=1e-9, atol=0)
                else:
                    equal = a[use].astype(str).tolist() == b[use].astype(str).tolist()
                if not equal:
                    raise AssertionError(f"Figure source-data values differ: {figure}, {sheet}, {column}")
    print("Verified source data for all six figures.")


def export_figure_source_data(package: Path, destination: Path) -> None:
    """Populate the formatted source-data workbooks from current analysis results."""
    from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
    import xml.etree.ElementTree as ET
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    relns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    def tag(name):
        return "{" + ns + "}" + name
    def column_name(index):
        value = ""
        while index:
            index, digit = divmod(index - 1, 26)
            value = chr(65 + digit) + value
        return value
    destination.mkdir(parents=True, exist_ok=True)
    for figure, sheets in figure_source_tables(package).items():
        name = "Source_Data_" + figure + ".xlsx"
        with ZipFile(package / "figures/source_data" / name) as source:
            contents = {n: source.read(n) for n in source.namelist()}
        workbook = ET.fromstring(contents["xl/workbook.xml"])
        relationships = {r.get("Id"): r.get("Target") for r in ET.fromstring(contents["xl/_rels/workbook.xml.rels"])}
        paths = {sheet.get("name"): relationships[sheet.get("{" + relns + "}id")] for sheet in workbook.find(tag("sheets"))}
        for sheet, frame in sheets.items():
            target = paths[sheet]
            target = target.lstrip("/") if target.startswith("/") else "xl/" + target
            root = ET.fromstring(contents[target])
            old = root.find(tag("sheetData")); position = list(root).index(old); root.remove(old)
            data = ET.Element(tag("sheetData"))
            for row_index, values in enumerate([list(frame.columns)] + list(frame.itertuples(index=False, name=None)), 1):
                row = ET.SubElement(data, tag("row"), r=str(row_index))
                for column_index, value in enumerate(values, 1):
                    if pd.isna(value):
                        continue
                    cell = ET.SubElement(row, tag("c"), r=column_name(column_index) + str(row_index))
                    if isinstance(value, (bool, np.bool_)):
                        cell.set("t", "b"); ET.SubElement(cell, tag("v")).text = "1" if value else "0"
                    elif isinstance(value, (int, float, np.number)):
                        ET.SubElement(cell, tag("v")).text = repr(value.item() if hasattr(value, "item") else value)
                    else:
                        cell.set("t", "inlineStr")
                        ET.SubElement(ET.SubElement(cell, tag("is")), tag("t")).text = str(value)
            root.insert(position, data)
            dimension = root.find(tag("dimension"))
            if dimension is not None:
                dimension.set("ref", "A1:" + column_name(len(frame.columns)) + str(len(frame) + 1))
            contents[target] = ET.tostring(root, encoding="UTF-8", xml_declaration=True)
        with ZipFile(destination / name, "w", ZIP_DEFLATED) as target:
            for path, data in contents.items():
                target.writestr(ZipInfo(path, (1980, 1, 1, 0, 0, 0)), data, compress_type=ZIP_DEFLATED)



def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, default=Path(__file__).absolute().parents[2])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    package = args.package.absolute()
    output = (args.output or package / "figures").absolute()
    figure1(package, output)
    figure2(package, output)
    figure3(package, output)
    figure4(package, output)
    figure5(package, output)
    extended_figure1(package, output)
    verify_figure_source_data(package)
    export_figure_source_data(package, output / "source_data")
    verify_figure_source_data(package, output / "source_data")
    print(f"Regenerated five main figures and one extended-data figure in {output}")


if __name__ == "__main__":
    main()
