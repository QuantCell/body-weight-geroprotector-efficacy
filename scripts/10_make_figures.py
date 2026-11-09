"""Regenerate the five manuscript figures and the Extended Data figure, with their source data, from the package inputs and outputs."""
from __future__ import annotations

import argparse
import io
import json
import os
import re
import zipfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.patches import Ellipse, PathPatch
from matplotlib.path import Path as MPath
from PIL import Image

from common import DERIVED, read_result


BLUE = "#2F7ED8"
ORANGE = "#F26B38"
AMBER = "#E69F00"
DARK_ORANGE = "#A8471C"
INK = "#1A1A1A"
SLATE = "#55534E"
DEVELOPMENTAL_BLUE = "#1F5FA9"
DEVELOPMENTAL_ORANGE = "#D4692A"
CONTROL, TREATED = DEVELOPMENTAL_BLUE, DEVELOPMENTAL_ORANGE
MALE = "#00798C"
STRAIN = "#2E9E6E"
FOUNDERS = ["BALB/cByJ", "C57BL/6J", "C3H/HeJ", "DBA/2J"]
SEX_SYMBOL = {"Male": "♂", "Female": "♀"}
# one type scale for every figure (all are shown at the same reduction in the manuscript)
LETTER, TITLE, LABEL, TICK, TEXT, NOTE, HEAD = 11, 8.8, 8.8, 8.4, 8.4, 7.8, 9.6


def save_figure(fig, output: Path, stem: str, dpi: int = 300, tight: bool = True) -> None:
    output.mkdir(parents=True, exist_ok=True)
    bbox = "tight" if tight else None
    fig.savefig(output / f"{stem}.png", dpi=dpi, bbox_inches=bbox, facecolor="white")
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=600, bbox_inches=bbox, facecolor="white")
    image = Image.open(buffer).convert("RGBA")
    flat = Image.alpha_composite(Image.new("RGBA", image.size, "white"), image).convert("RGB")
    flat.save(output / f"{stem}.tif", compression="tiff_lzw", dpi=(600, 600))
    fig.savefig(output / f"{stem}.pdf", bbox_inches=bbox, facecolor="white", metadata={"CreationDate": None})
    plt.close(fig)


def save_source_data(path: Path, sheets: dict) -> None:
    """Write one source-data workbook, one sheet per panel."""
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as book:
        book.book.properties.creator = "Authors"
        for name, frame in sheets.items():
            frame.to_excel(book, sheet_name=name, index=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(buffer) as source, zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as target:
        for item in source.infolist():
            data = source.read(item.filename)
            if item.filename == "docProps/core.xml":
                data = re.sub(rb"<dcterms:(created|modified)\b.*?</dcterms:\1>", b"", data)
            target.writestr(zipfile.ZipInfo(item.filename, (1980, 1, 1, 0, 0, 0)), data, zipfile.ZIP_DEFLATED)


SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def p_text(p: float) -> str:
    """P as the manuscript prints it."""
    if p < 0.001:
        mantissa, exponent = f"{p:.1e}".split("e")
        return f"P = {mantissa} × 10{str(int(exponent)).translate(SUPERSCRIPT)}"
    return f"P = {p:.2g}"


def base_style(font_size: float = TICK) -> None:
    """Reset every property the figures depend on."""
    shared = font_size == TICK
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.family": ["Arial", "DejaVu Sans"], "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.size": font_size, "axes.labelsize": LABEL if shared else "medium",
        "xtick.labelsize": TICK if shared else "medium", "ytick.labelsize": TICK if shared else "medium",
        "legend.fontsize": TEXT if shared else "medium",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": True, "axes.spines.bottom": True,
        "axes.edgecolor": "#6E6B64", "axes.labelcolor": INK,
        "xtick.color": "#333333", "ytick.color": "#333333",
        "xtick.major.size": 0, "ytick.major.size": 0,
        "xtick.minor.size": 0, "ytick.minor.size": 0,
    })


def panel_title(ax, letter: str, title: str, pad: float = 6) -> None:
    """Bold panel letter followed by a plain title, the same size in every figure."""
    if letter:
        ax.annotate(letter, (0, 1), xycoords="axes fraction", xytext=(0, pad), textcoords="offset points",
                    ha="left", va="bottom", fontsize=LETTER, fontweight="bold", color=INK)
    ax.annotate(title, (0, 1), xycoords="axes fraction", xytext=(13, pad), textcoords="offset points",
                ha="left", va="bottom", fontsize=TITLE, color=INK)


def standardized_km(rows: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Kaplan-Meier survival from six months, with each experiment and site weighted equally."""
    entry, time = rows.entry.to_numpy(), rows.lifespan_days.to_numpy()
    dead, weight = rows.dead.to_numpy().astype(bool), rows.standardization_weight.to_numpy()
    times, survival, s = [entry.min()], [1.0], 1.0
    for t in np.unique(time[dead]):
        s *= 1 - weight[(time == t) & dead].sum() / weight[(entry < t) & (time >= t)].sum()
        times.append(t)
        survival.append(s)
    return np.array(times), np.array(survival)


def km_median(times: np.ndarray, survival: np.ndarray) -> float:
    """Median as survival::survfit reports it, including the midpoint of a flat stretch at one half."""
    tolerance = np.finfo(float).eps ** 0.5
    i = int(np.argmax(survival < 0.5 + tolerance))
    below = np.flatnonzero(survival < survival[i])
    if abs(survival[i] - 0.5) < tolerance and below.size:
        return float((times[i] + times[below[0]]) / 2)
    return float(times[i])


QUARTER_NAMES = ["Lightest", "Second", "Third", "Heaviest"]
MEDIAN_COLUMN = {0: "standardized_control_median", 1: "standardized_treated_median"}


def mouse_icons(ax, xs, y, height, fats):
    """Same-length mouse drawings of increasing fatness centred on data x positions, at axes-fraction height y."""
    from matplotlib.transforms import blended_transform_factory

    blend = blended_transform_factory(ax.transData, ax.transAxes)
    span = np.diff(ax.get_xlim())[0]
    width = 0.19 * span
    for x, fat in zip(xs, fats):
        draw_mouse(ax.inset_axes([x - width / 2, y - height / 2, width, height], transform=blend), 0.35, fat)


def figure1(package: Path, output: Path) -> dict:
    """Two scenarios for how treatment benefit could depend on weight, and what the matched quartiles show."""
    rows = pd.read_csv(DERIVED / "quartile_comparisons.csv")
    frame = read_result("matched_quartile_survival").set_index("quartile").sort_index()
    curves = {}
    for q in range(1, 5):
        for tx in (0, 1):
            curves[q, tx] = standardized_km(rows[(rows.matched_quartile == q) & (rows.Tx == tx)])
            if km_median(*curves[q, tx]) != frame.loc[q, MEDIAN_COLUMN[tx]]:
                raise AssertionError(f"Survival-curve median differs from the quartile table: quartile {q}, treated {tx}")
    pd.concat([pd.DataFrame({"quartile": q, "group": "treated" if tx else "control", "age_days": times,
                             "surviving": survival}) for (q, tx), (times, survival) in curves.items()]
              ).to_csv(package / "data" / "outputs" / "figure_1_survival_curves.csv", index=False, lineterminator="\n")
    control, treated = frame[MEDIAN_COLUMN[0]].to_numpy(), frame[MEDIAN_COLUMN[1]].to_numpy()
    fats = [0, 1 / 3, 2 / 3, 1]
    x = np.arange(4)

    def two_lines(ax, ctl, trt, gains, avoid=None):
        ax.plot(x, ctl, "o-", color=CONTROL, lw=2.4, ms=6.5, zorder=3)
        ax.plot(x, trt, "o-", color=TREATED, lw=2.4, ms=6.5, zorder=3)
        for i, gain in zip(x, gains):
            ax.annotate("", (i, trt[i]), (i, ctl[i]), zorder=4,
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2, shrinkA=5, shrinkB=5))
            if gain:
                y = trt[i] + 16
                if avoid is not None and y - 4 < avoid < y + 22:
                    y = avoid + 6
                ax.text(i + 0.05, y, gain, ha="left", va="bottom", fontsize=NOTE, color=INK, fontweight="bold",
                        bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"), zorder=5)
        ax.text(3.25, trt[3], "Treated", ha="left", va="center", fontsize=TEXT, fontweight="bold", color=TREATED)
        ax.text(3.25, ctl[3], "Control", ha="left", va="center", fontsize=TEXT, fontweight="bold", color=CONTROL)
        ax.set_xlim(-0.45, 4.45)
        ax.set_xticks([])
        mouse_icons(ax, x, -0.09, 0.15, fats)
        ax.text(1.5, -0.25, "Heavier mice →", transform=ax.get_xaxis_transform(), ha="center", va="center",
                fontsize=TEXT, color=INK)

    def lightest_line(ax, level):
        ax.plot([-0.2, 3.2], [level, level], color=CONTROL, lw=0.9, ls=(0, (3, 2)), zorder=1)
        ax.text(3.3, level, "lightest\ncontrol", ha="left", va="center", fontsize=NOTE - 0.6, color=CONTROL)

    def scenario(ax, letter, title, ctl, trt, caption, below=False):
        two_lines(ax, np.array(ctl), np.array(trt), [""] * 4)
        if below:
            lightest_line(ax, ctl[0])
        ax.set_ylim(600, 1100)
        ax.set_yticks([])
        ax.set_ylabel("Lifespan →")
        ax.text(0.03, 0.97, caption, transform=ax.transAxes, ha="left", va="top", fontsize=TEXT, color=INK,
                fontweight="bold")
        panel_title(ax, letter, title)

    base_style()
    # same canvas as the survival-only figure this one replaced, so it fills the same space in the manuscript
    fig = plt.figure(figsize=(2109 / 300, 1972 / 300))
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 1.1], hspace=0.5, left=0.075, right=0.952, top=0.955, bottom=0.075)
    # equal panel sizes; the wider second gap holds panel b's axis labels
    top = grid[0].subgridspec(1, 5, width_ratios=[1, 0.39, 1, 0.61, 1], wspace=0)
    scenario(fig.add_subplot(top[0, 0]), "a", "Scenario 1: ceiling raised",
             [880, 850, 820, 790], [990, 960, 930, 900], "Same gain at every weight")
    scenario(fig.add_subplot(top[0, 2]), "", "Scenario 2: risk rescued",
             [890, 840, 765, 675], [945, 905, 860, 810], "Gain grows with weight;\nheavy treated mice stay\nbelow light control mice",
             below=True)
    ax_b = fig.add_subplot(top[0, 4])
    two_lines(ax_b, control, treated, [f"+{t - c:.0f} d" for c, t in zip(control, treated)], avoid=control[0])
    lightest_line(ax_b, control[0])
    ax_b.set_ylim(600, 1100)
    ax_b.set_yticks(range(600, 1101, 100))
    ax_b.set_ylabel("Median lifespan (days)")
    ax_b.grid(axis="y", color="#ECEAE4", lw=0.7, zorder=0)
    panel_title(ax_b, "b", "Observed: ITP males")
    for i, w in zip(x, frame.mean_control_weight):
        ax_b.text(i, -0.17, f"{w:.0f} g", transform=ax_b.get_xaxis_transform(), ha="center", va="center", fontsize=NOTE,
                  color=INK)

    surv = grid[1].subgridspec(1, 2, wspace=0.12)
    first = None
    sheets = {}
    for col, (letter, q) in enumerate((("c", 1), ("d", 4))):
        ax = fig.add_subplot(surv[0, col], sharey=first)
        first = first or ax
        for tx, colour, name in ((0, CONTROL, "Control"), (1, TREATED, "Treated")):
            times, survival = curves[q, tx]
            ax.step(times, survival, where="post", color=colour, lw=2, label=name, zorder=3)
        mc, mt = frame.loc[q, MEDIAN_COLUMN[0]], frame.loc[q, MEDIAN_COLUMN[1]]
        ax.annotate("", (mt, 0.5), (mc, 0.5), zorder=4,
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.4, shrinkA=0, shrinkB=0))
        ax.text(0.97, 0.90, "Median increase", transform=ax.transAxes, ha="right", fontsize=TEXT, color=SLATE)
        ax.text(0.97, 0.82, f"+{mt - mc:.0f} days (+{frame.loc[q, 'gain_pct']:.1f}%)", transform=ax.transAxes,
                ha="right", fontsize=HEAD, fontweight="bold", color=INK)
        ax.axvline(500, color=INK, lw=0.8, ls=(0, (1, 2)), zorder=1)
        ax.text(510, 0.04, "day 500", fontsize=NOTE, color=INK)
        ax.set_xlim(150, 1450)
        ax.set_xticks(range(200, 1401, 200))
        ax.set_ylim(0, 1.03)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "25", "50", "75", "100"])
        ax.set_xlabel("Age (days)")
        ax.grid(color="#ECEAE4", lw=0.7, zorder=0)
        if col == 0:
            ax.set_ylabel("Alive (%)")
            ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.29, 0.11))
        else:
            ax.tick_params(labelleft=False)
        panel_title(ax, letter, f"{QUARTER_NAMES[q - 1]} quartile")
        draw_mouse(ax.inset_axes([0.77, 0.6, 0.2, 0.15]), 0.4, fats[q - 1])
        sheets[f"Fig1{letter}"] = pd.concat([pd.DataFrame({"group": "treated" if tx else "control", "age_days": t,
                                                           "alive_percent": 100 * s})
                                            for tx in (0, 1) for t, s in [curves[q, tx]]], ignore_index=True)
    save_figure(fig, output, "Figure_1_quartile_survival", tight=False)
    medians = pd.DataFrame({"quartile": QUARTER_NAMES, "mean_control_weight_g": frame.mean_control_weight.round(2).to_numpy(),
                            "control_median_days": control, "treated_median_days": treated,
                            "increase_days": treated - control, "increase_percent": frame.gain_pct.round(2).to_numpy()})
    return {"Fig1b": medians, **sheets}


def figure2(package: Path, output: Path) -> dict:
    """Weight dependence of the treatment effect by sex and intervention class."""
    from matplotlib.transforms import blended_transform_factory

    primary = read_result("primary_interaction").set_index("sex")
    split = read_result("weight_class").set_index("analysis")
    trend = read_result("quartile_trend").iloc[0]
    estimates = [("Males: all experiments", primary.loc["m"]),
                 ("Males: weight-lowering", split.loc["male_weight_lowering"]),
                 ("Males: non-weight-lowering", split.loc["male_nonlowering"]),
                 ("Females: all experiments", primary.loc["f"])]

    base_style()
    # same canvas as the figure this one replaced
    fig = plt.figure(figsize=(2165 / 300, 1039 / 300))
    ax_a = fig.add_axes([0.015, 0.33, 0.415, 0.57])
    ax_b = fig.add_axes([0.54, 0.16, 0.28, 0.74])

    for y, (label, est) in zip([3, 2, 1, 0], estimates):
        ax_a.plot([est.lo, est.hi], [y, y], color=MALE, lw=2, zorder=2)
        ax_a.scatter(est.HR, y, s=58, facecolors=MALE if est.p < 0.05 else "white", edgecolors=MALE, linewidths=1.7,
                     zorder=3)
        ax_a.text(0.835, y + 0.31, label, ha="left", va="center", fontsize=TEXT, color=INK)
        ax_a.text(1.105, y + 0.31, f"{est.HR:.2f} ({est.lo:.2f}–{est.hi:.2f})", ha="right", va="center",
                  fontsize=NOTE, color=MALE)
        if est.p < 0.05:
            ax_a.text(1.105, y - 0.02, f"{100 * (1 - est.HR):.0f}% more favourable", ha="right", va="center",
                      fontsize=NOTE, color=MALE, bbox=dict(facecolor="white", edgecolor="none", pad=1))
    ax_a.axvline(1, color=INK, ls="--", lw=0.9, zorder=1)
    ax_a.text(1, 3.66, "1 = no change", ha="center", va="center", fontsize=NOTE, color=INK,
              bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    ax_a.text(0.97, -0.5, "Filled: P < 0.05; open: P ≥ 0.05", ha="center", va="center", fontsize=NOTE, color=INK,
              bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    ax_a.set_xlim(0.83, 1.11)
    ax_a.set_ylim(-0.7, 3.95)
    ax_a.set_yticks([])
    ax_a.set_xticks([0.85, 0.9, 1, 1.1], ["0.85", "0.90", "1.00", "1.10"])
    ax_a.set_xlabel("Treatment HR change per 1 s.d. higher weight", labelpad=6)
    ax_a.spines["left"].set_visible(False)
    blend = blended_transform_factory(ax_a.transData, ax_a.transAxes)
    for start, end, x, label in ((0.985, 0.85, 0.918, "Greater benefit\nin heavier mice"),
                                 (1.015, 1.1, 1.058, "Greater benefit\nin lighter mice")):
        ax_a.annotate("", xy=(end, -0.26), xytext=(start, -0.26), xycoords=blend, textcoords=blend,
                      arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.4), annotation_clip=False)
        ax_a.text(x, -0.3, label, ha="center", va="top", fontsize=TEXT, color=INK, transform=blend, clip_on=False)
    panel_title(ax_a, "a", "Heavier males benefit more")

    sensitivity = read_result("quartile_sensitivities")
    lines = [("weight-lowering standardized quarters", DEVELOPMENTAL_BLUE, 2.8, 6.5),
             ("non-lowering conventional quartiles", DEVELOPMENTAL_ORANGE, 2.6, 6)]
    for analysis, colour, width, size in lines:
        block = sensitivity[sensitivity.analysis.eq(analysis)].sort_values("bin")
        position = (block.bin.to_numpy() - 0.5) / block.bins.to_numpy()
        ax_b.plot(position, block.gain_pct, marker="o", ms=size, color=colour, lw=width, zorder=3)
    ax_b.text(1.04, 20.2, "Weight-lowering\nexperiments", ha="left", va="center", fontsize=TEXT, color=DEVELOPMENTAL_BLUE,
              fontweight="bold", transform=ax_b.get_yaxis_transform(), clip_on=False)
    ax_b.text(1.04, 6.4, "Non-weight-\nlowering", ha="left", va="center", fontsize=TEXT, color=DEVELOPMENTAL_ORANGE, fontweight="bold",
              transform=ax_b.get_yaxis_transform(), clip_on=False)
    ax_b.set_xlim(0, 1)
    ax_b.set_xticks([0.125, 0.875], ["Lightest", "Heaviest"])
    ax_b.set_ylim(-1, 24)
    ax_b.axhline(0, color=INK, lw=0.8)
    ax_b.set_ylabel("Median lifespan increase (%)")
    ax_b.set_xlabel("Pre-treatment weight quartile")
    ax_b.grid(axis="y", color="#ECEAE4", lw=0.7, zorder=0)
    ax_b.text(0.03, 0.97, f"Weight-lowering trend: HR {trend.HR_per_quartile:.2f}\nper heavier quartile, {p_text(trend.p)}",
              transform=ax_b.transAxes, va="top", fontsize=NOTE, color=DEVELOPMENTAL_BLUE)
    panel_title(ax_b, "b", "Gradient mainly in weight-lowering experiments")
    save_figure(fig, output, "Figure_2_matched_quartiles", tight=False)
    table = pd.DataFrame([{"analysis": label, "interaction_HR_per_SD": est.HR, "lower_95": est.lo, "upper_95": est.hi,
                           "p": est.p} for label, est in estimates])
    groups = sensitivity[sensitivity.analysis.isin([line[0] for line in lines])][
        ["analysis", "bins", "bin", "mean_control_weight", "control_median", "treated_median", "gain_pct"]]
    return {"Fig2a": table, "Fig2b": groups,
            "Fig2b_trend": pd.DataFrame({"quantity": ["trend HR per heavier quartile", "lower 95% CI", "upper 95% CI", "P"],
                                         "value": [trend.HR_per_quartile, trend.lo, trend.hi, trend.p]})}


def load_strain_means(path: Path, sex: str) -> dict[str, float]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    rows = raw.get("strainmeans", raw) if isinstance(raw, dict) else raw
    return {row["strain"]: float(row["mean"]) for row in rows
            if str(row.get("sex", "")).lower().startswith(sex) and row.get("mean") is not None}


def cubic(points, n: int = 40) -> np.ndarray:
    """Points along a cubic Bezier curve."""
    t = np.linspace(0, 1, n)[:, None]
    p = np.asarray(points, float)
    return (1 - t) ** 3 * p[0] + 3 * (1 - t) ** 2 * t * p[1] + 3 * (1 - t) * t ** 2 * p[2] + t ** 3 * p[3]


def closed_path(start, segments) -> MPath:
    """Closed path of cubic Bezier segments."""
    verts, codes = [tuple(start)], [MPath.MOVETO]
    for segment in segments:
        verts += [tuple(point) for point in segment]
        codes += [MPath.CURVE4] * 3
    return MPath(verts + [tuple(start)], codes + [MPath.CLOSEPOLY])


def draw_mouse(ax, width_in: float, fat: float = 0.0) -> None:
    """Side view of a laboratory mouse; fat (0 to 1) deepens the trunk at constant length."""
    lw = lambda points: max(points * width_in / 9.6, 0.1)

    def plump(points):
        xy = np.array(points, float).reshape(-1, 2)
        x, y = xy[:, 0], xy[:, 1]
        trunk = np.clip((10.0 - x) / 1.4, 0, 1)
        rise = 1.3 * np.exp(-((x - 6.2) / 2.4) ** 2) * trunk
        sag = 0.7 * np.exp(-((x - 6.3) / 1.9) ** 2)
        y = y + fat * np.where(y > 4.1, rise, -sag)
        x = x - fat * 0.4 * np.clip(4.8 - x, 0, 1)
        return np.c_[x, y]
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
    ax.add_patch(Ellipse(plump([(9.25, 5.98)])[0], 0.56, 0.66, angle=-12, fc="#CF9A90", ec=pink_edge, lw=lw(0.5),
                         zorder=1.5))

    body = [(3.85, 3.15), [
        ((3.55, 4.1), (4.0, 4.95), (4.8, 5.15)), ((5.5, 5.5), (6.4, 5.85), (7.3, 5.8)),
        ((8.2, 5.78), (8.9, 6.0), (9.5, 5.85)), ((10.2, 5.65), (10.85, 4.85), (11.05, 4.47)),
        ((11.1, 4.2), (10.8, 3.98), (10.4, 3.95)), ((9.85, 3.9), (9.25, 3.65), (8.95, 3.3)),
        ((8.72, 3.0), (8.62, 2.72), (8.4, 2.6)), ((7.6, 2.38), (6.6, 2.3), (5.9, 2.35)),
        ((5.0, 2.38), (4.45, 2.5), (4.12, 2.75)), ((3.9, 2.92), (3.8, 3.03), (3.85, 3.15))]]
    body = closed_path(plump(body[0])[0], [plump(segment) for segment in body[1]])
    outline = PathPatch(body, fc="none", ec="none", zorder=2)
    ax.add_patch(outline)
    stops = [(0, "#A39A90"), (0.14, "#877D73"), (0.35, "#6E645B"), (0.7, "#625850"), (1, "#554C45")]
    level = np.linspace(0, 1, 400)
    ramp = np.stack([np.interp(level, [s for s, _ in stops], [to_rgb(c)[k] for _, c in stops]) for k in range(3)], -1)
    image = ax.imshow(np.repeat(ramp[:, None, :], 400, axis=1), extent=(3.5 - 0.4 * fat, 11.2, 2.25 - 0.75 * fat, 6.1 + 1.35 * fat),
                      origin="lower",
                      interpolation="bicubic", aspect="auto", zorder=2)
    image.set_clip_path(outline)
    haunch = dict(xy=(5.25 - 0.2 * fat, 3.75 + 0.2 * fat), width=2.6 * (1 + 0.15 * fat), height=2.45 * (1 + 0.5 * fat),
                  angle=10)
    for shape in (Ellipse(**haunch, fc="#7A7067", ec="none", alpha=0.55, zorder=3),
                  Ellipse(**haunch, fc="none", ec=fur_line, lw=lw(0.6), alpha=0.6, zorder=3),
                  Ellipse((10.0, 4.45), 1.7, 1.0, angle=-20, fc="#80766D", ec="none", alpha=0.4, zorder=3)):
        ax.add_patch(shape)
        shape.set_clip_path(outline)
    ax.plot(*cubic(plump([(8.15, 4.6), (8.0, 4.0), (8.15, 3.35), (8.55, 2.75)])).T, color=fur_line, lw=lw(0.6),
            alpha=0.55, zorder=3, solid_capstyle="round")
    ax.add_patch(PathPatch(body, fc="none", ec=fur_edge, lw=lw(0.6), zorder=4))

    for heel, toe, t in (((5.35, 2.42), (6.95, 2.18), 0.17), ((8.45, 2.75), (9.3, 2.42), 0.12),
                         ((8.9, 2.95), (9.75, 2.62), 0.11)):
        heel, toe = plump([heel, toe])
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

    lift = plump([(8.72, 5.68)])[0, 1] - 5.68
    ax.add_patch(Ellipse((8.72, 5.68 + lift), 0.78, 0.98, angle=-18, fc="#DCA79C", ec=pink_edge, lw=lw(0.6), zorder=6))
    ax.add_patch(Ellipse((8.76, 5.64 + lift), 0.52, 0.7, angle=-18, fc="#C98A80", ec="none", zorder=6))
    ax.add_patch(Ellipse((8.66, 5.7 + lift), 0.32, 0.52, angle=-18, fc="#DCA79C", ec="none", alpha=0.55, zorder=6))
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
    ax.set_ylim(0.8, 7.2)
    ax.set_aspect("equal")
    ax.axis("off")



def figure3(package: Path, output: Path) -> dict:
    """ITP control weight, control-group lifespan and adiposity against reference populations, as labelled strips."""
    inputs, outputs = package / "data" / "inputs", package / "data" / "outputs"
    book = package / "supplementary" / "Supplementary_Data.xlsx"
    weight = load_strain_means(inputs / "mpd_32603_strainmeans.json", "m")
    age = load_strain_means(inputs / "mpd_32680_strainmeans.json", "m")
    matched = {strain: value for strain, value in weight.items()
               if strain in age and age[strain] <= 36.0}
    if len(matched) != 27 or not all(founder in matched for founder in FOUNDERS):
        raise AssertionError(f"Expected 27 age-matched strains including all founders, found {len(matched)}")
    survey_mean = float(np.mean(list(matched.values())))
    quartiles = read_result("matched_quartile_survival").sort_values("quartile")
    lightest, heaviest = quartiles.mean_control_weight.astype(float).iloc[[0, -1]]
    founder_weights = read_result("founder_weight_comparison").set_index("quantity")
    itp_weight = float(founder_weights.loc["itp_male_control_6m_mean_weight_g", "value"])

    body = pd.read_excel(book, sheet_name="SD5_body_composition")
    rapa = body[(body["Intervention"] == "Rapamycin") & (body["strain"] == "UM-HET3") &
                body["dosage"].astype(str).str.contains("42 ppm", regex=False)]
    adiposity = rapa.set_index("sex")["bodyfat_pct_control"].astype(float)
    if not {"Male", "Female"}.issubset(adiposity.index):
        raise AssertionError("Missing UM-HET3 control adiposity values")

    experiments = pd.read_excel(book, sheet_name="SD10_control_lifespan")
    experiments = experiments[experiments.sex.isin(["Female", "Male"])].copy()
    experiments["is_itp"] = experiments.itp_experiment.astype(str).str.upper().eq("ITP")
    groups = experiments.drop_duplicates(["pmid", "sex", "lifespan_control_days"])
    survey_lifespan = pd.read_csv(inputs / "inbred_strain_median_lifespan.csv").set_index("strain")

    DOT, FOUNDER_DOT, ITP_BAR, MUTED = "#2A78D6", "#1BAF7A", "#EB6834", "#52514E"
    rng = np.random.default_rng(3)

    def strip(ax, letter, title, values, founders, bars, x_label, xlim, note):
        """Jittered reference dots, labelled founder strains above and below, and ITP values as bars."""
        ax.scatter(values, rng.uniform(-0.18, 0.18, len(values)), s=72, color=DOT, alpha=0.75, edgecolors="white",
                   linewidths=0.6, zorder=3)
        order = sorted(founders.items(), key=lambda item: item[1])
        gap = 0.17 * (xlim[1] - xlim[0])
        for side, group in ((1, order[0::2]), (-1, order[1::2])):
            xs = [value for _, value in group]
            spread = list(xs)
            for _ in range(2):
                for k in range(1, len(spread)):
                    spread[k] = max(spread[k], spread[k - 1] + gap)
                shift = (sum(spread) - sum(xs)) / len(spread)
                spread = [x - shift for x in spread]
            for k in range(1, len(spread)):
                spread[k] = max(spread[k], spread[k - 1] + gap)
            for (name, value), x in zip(group, spread):
                ax.scatter([value], [side * 0.32], s=88, color=FOUNDER_DOT, edgecolors=INK, linewidths=0.8, zorder=4)
                ax.annotate(name, (value, side * 0.32), xytext=(x, side * 0.95), ha="center",
                            va="bottom" if side > 0 else "top", fontsize=NOTE, color=INK,
                            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=1, shrinkB=5))
        for label, value in bars:
            ax.plot([value, value], [-0.7, 0.7], color=ITP_BAR, lw=3.2, solid_capstyle="round", zorder=5)
            ax.text(value, 1.3, label, ha="center", va="bottom", fontsize=NOTE, fontweight="bold", color=INK)
        ax.set_ylim(-1.8, 1.95)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlim(*xlim)
        ax.set_xlabel(x_label)
        ax.grid(axis="x", color="#E4E3DF", lw=0.7, zorder=0)
        ax.text(0.99, 0.02, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=NOTE, color=MUTED)
        ax.text(-0.01, 1.22, letter, transform=ax.transAxes, ha="right", va="top", fontsize=LETTER, fontweight="bold")
        ax.text(0, 1.22, title, transform=ax.transAxes, ha="left", va="top", fontsize=TITLE)

    def reach_text(count, total, noun, target):
        return f"{count} of {total} {noun} {'reaches' if count == 1 else 'reach'} the {target}"

    base_style()
    # same proportions as the figure this one replaced
    fig = plt.figure(figsize=(7.2, 7.2 * 2121 / 1919))
    key = fig.add_axes([0.0, 0.955, 1.0, 0.04])
    key.axis("off")
    key.set_xlim(0, 1)
    key.set_ylim(-1, 1)
    key.scatter([0.03], [0.5], s=72, color=DOT, alpha=0.75)
    key.text(0.05, 0.5, "Other laboratory strains (a, d, e) or non-ITP control groups (b, c)", va="center", fontsize=TEXT)
    key.scatter([0.03], [-0.6], s=88, color=FOUNDER_DOT, edgecolors=INK, linewidths=0.8)
    key.text(0.05, -0.6, "The four UM-HET3 founder strains", va="center", fontsize=TEXT)
    key.plot([0.6, 0.62], [-0.6, -0.6], color=ITP_BAR, lw=3.2)
    key.text(0.635, -0.6, "ITP controls", va="center", fontsize=TEXT)

    # a: male body weight
    names = list(matched)
    values = np.array([matched[name] for name in names])
    is_founder = np.array([name in FOUNDERS for name in names])
    heavy = int((values >= itp_weight).sum())
    strip(fig.add_axes([0.04, 0.71, 0.94, 0.17]), "a",
          "ITP males were heavier than every founder strain and almost every laboratory strain", values[~is_founder],
          {name: matched[name] for name in FOUNDERS},
          [(f"ITP lightest\nquartile\n{lightest:.1f} g", lightest), (f"ITP mean\n{itp_weight:.1f} g", itp_weight),
           (f"ITP heaviest\nquartile\n{heaviest:.1f} g", heaviest)],
          "Male body weight (g)", (17, 49), reach_text(heavy, len(values), "strains", "ITP mean"))
    summary = [("age_matched_strains", len(matched)), ("survey_mean_weight_g", survey_mean),
               ("male_weight_strains", len(values)),
               ("male_weight_at_or_above_lightest_quarter", int((values >= lightest).sum())),
               ("male_weight_at_or_above_itp", heavy)]
    sheets = {"Fig3a": pd.DataFrame({"strain": names, "body_weight_g": values.round(2),
                                     "age_weeks": [age[n] for n in names], "UM_HET3_founder": is_founder}),
              "Fig3a_reference_lines": pd.DataFrame({
                  "reference": ["ITP lightest quartile", "ITP male controls, six months", "ITP heaviest quartile"],
                  "body_weight_g": np.round([lightest, itp_weight, heaviest], 2)})}

    # b, c: control-group lifespan
    for rect, letter, sex, title in (([0.04, 0.395, 0.44, 0.19], "b", "Male", "Male ITP controls lived as long as\nother control groups"),
                                     ([0.54, 0.395, 0.44, 0.19], "c", "Female", "Female ITP controls outlived every\nnon-ITP control group")):
        block = groups[groups.sex.eq(sex)]
        field = block[~block.is_itp].lifespan_control_days.to_numpy()
        itp_groups = block[block.is_itp].lifespan_control_days.to_numpy()
        itp_median = float(np.median(itp_groups))
        column = f"{sex.lower()}_median_days"
        strip(fig.add_axes(rect), letter, title, field,
              {name: float(survey_lifespan.loc[name, column]) for name in FOUNDERS},
              [(f"ITP median\n{int(itp_median + 0.5)} days", itp_median)], f"{sex} control-group lifespan (days)",
              (380, 1020), reach_text(int((field >= itp_median).sum()), len(field), "non-ITP groups", "ITP median"))
        studies = block[~block.is_itp][["first_author", "year", "pmid", "strain", "lifespan_control_days"]]
        founder_rows = pd.DataFrame({"first_author": "Yuan et al. (strain survey)", "year": 2020, "strain": FOUNDERS,
                                     "lifespan_control_days": [float(survey_lifespan.loc[n, column]) for n in FOUNDERS]})
        sheets[f"Fig3{letter}"] = pd.concat([studies.assign(group="non-ITP control group"),
                                             founder_rows.assign(group="UM-HET3 founder strain"),
                                             pd.DataFrame({"lifespan_control_days": itp_groups, "group": "ITP control group"})],
                                            ignore_index=True)
        sheets[f"Fig3{letter}_ITP_median"] = pd.DataFrame({"quantity": [f"ITP {sex.lower()} control median (days)"],
                                                           "value": [itp_median]})

    # d, e: body fat
    for rect, letter, sex, code in (([0.04, 0.07, 0.44, 0.19], "d", "Male", "m"), ([0.54, 0.07, 0.44, 0.19], "e", "Female", "f")):
        by_strain = load_strain_means(inputs / "mpd_10331_strainmeans.json", code)
        names = list(by_strain)
        values = np.array([by_strain[name] for name in names])
        is_founder = np.array([name in FOUNDERS for name in names])
        itp_value = float(adiposity.loc[sex])
        reaching = int((values >= itp_value).sum())
        strip(fig.add_axes(rect), letter,
              f"ITP {sex.lower()}s carried more fat than\n{len(values) - reaching} of {len(values)} strains",
              values[~is_founder], {name: by_strain[name] for name in FOUNDERS}, [(f"ITP\n{itp_value:.1f}%", itp_value)],
              f"{sex} body fat (%)", (12, 47), reach_text(reaching, len(values), "strains", "ITP value"))
        sheets[f"Fig3{letter}"] = pd.DataFrame({"strain": names, "body_fat_percent": values.round(2),
                                                "UM_HET3_founder": is_founder})
        sheets[f"Fig3{letter}_ITP_value"] = pd.DataFrame({"quantity": [f"ITP {sex.lower()} control body fat (%)"],
                                                          "value": [itp_value]})
        summary += [(f"{sex.lower()}_body_fat_strains", len(values)), (f"{sex.lower()}_body_fat_at_or_above_itp", reaching)]

    pd.DataFrame(summary, columns=["quantity", "value"]).to_csv(outputs / "strain_context_summary.csv", index=False, lineterminator="\n")
    save_figure(fig, output, "Figure_3_strain_context", tight=False)
    return sheets

COMPOUND_NAMES = {"17alpha-estradiol": "17α-estradiol", "2-(2-hydroxyphenyl)benzoxazole": "HBX"}


def signed_percent(value: float) -> str:
    return f"{value:+.1f}%".replace("-", "−")


DRUG_CLASS = {"16-hydroxyestradiol": "estradiol metabolite", "17alpha-estradiol": "non-feminizing estrogen",
              "Acarbose": "α-glucosidase inhibitor", "Canagliflozin": "SGLT2 inhibitor", "Rapamycin": "mTOR inhibitor",
              "MIF098": "MIF antagonist", "Glycine": "amino acid", "Meclizine": "antihistamine, mTORC1 inhibitor",
              "Candesartan": "angiotensin receptor blocker", "NDGA": "antioxidant",
              "Geranylgeranylacetone": "heat-shock protein inducer", "Resveratrol": "plant polyphenol",
              "2-(2-hydroxyphenyl)benzoxazole": "amyloid-binding compound", "Fisetin": "senolytic flavonoid"}


def figure4(package: Path, output: Path) -> dict:
    """Median-lifespan increase in the lightest and heaviest pre-treatment quartiles, one row per compound."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    base_style()
    table = read_result("compound_survival")
    common = read_result("compound_common_interaction").iloc[0]
    is_pooled = table.compound.str.contains("pool", case=False)
    pooled = table[is_pooled].iloc[0]
    compounds = table[~is_pooled].assign(mean=lambda d: (d.lowest_gain_pct + d.highest_gain_pct) / 2)
    kept = compounds.lowest_gain_pct >= 5
    LIGHT, HEAVY, BAND = "#0072B2", "#D55E00", "#F0F0F0"

    rows, y = [], -0.75
    for heading, block in (("Lifespan gain in the lightest mice (≥5%)", compounds[kept]),
                           ("Little or no lifespan gain in the lightest mice", compounds[~kept])):
        y -= 1.0
        rows.append(("heading", heading, y))
        y -= 0.15
        for row in block.sort_values("mean", ascending=False).itertuples(index=False):
            tagged = row.compound in DRUG_CLASS
            y -= 1.0
            rows.append(("compound", row, y))
            y -= 0.25 if tagged else 0
    fig = plt.figure(figsize=(7.2, 6.7))
    left, width = 0.34, 0.635
    ax = fig.add_axes([left, 0.085, width, 0.88])
    edge = (0.012 - left) / width                         # figure left margin in axes units
    ax.set_ylim(y - 0.6, 1.05)
    ax.set_xlim(-14, 36)

    def dumbbell(row, yy, size, labels=NOTE):
        low, high = row.lowest_gain_pct, row.highest_gain_pct
        ax.plot([low, high], [yy, yy], color="#D9D9D9", lw=size / 12, zorder=2,
                solid_capstyle="butt")
        ax.scatter([low], [yy], s=size * 1.3, color=LIGHT, zorder=3)
        ax.scatter([high], [yy], s=size * 1.3, facecolors="white", edgecolors=HEAVY, linewidths=2, zorder=3)
        sign = 1 if high >= low else -1
        for value, colour, side in ((low, LIGHT, -sign), (high, HEAVY, sign)):
            ax.annotate(signed_percent(value), (value, yy), xytext=(side * (6 + size / 14), 0), textcoords="offset points",
                        ha="left" if side > 0 else "right", va="center", fontsize=labels, fontweight="bold", color=colour,
                        bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none") if labels == NOTE else None,
                        zorder=4)

    ax.add_patch(Rectangle((edge, -0.95), 1 - edge + 0.02, 1.7, transform=ax.get_yaxis_transform(), color=BAND,
                           zorder=0, clip_on=False))
    dumbbell(pooled, 0.1, 110, HEAD)
    ax.text(-0.03, 0.1, f"All {len(compounds)} compounds ({int(pooled.experiments)})", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=TEXT, fontweight="bold", color=INK)
    ax.text(12.6, -0.55, f"Shared treatment × weight HR {common.estimate:.2f} ({common.lower_95:.2f}–{common.upper_95:.2f}); "
            f"{p_text(common.p)}", ha="center", va="center", fontsize=NOTE, color=INK)
    for kind, item, yy in rows:
        if kind == "heading":
            ax.text(edge, yy, item, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=TEXT,
                    fontweight="bold", color=INK)
            continue
        ax.text(-0.03, yy, f"{COMPOUND_NAMES.get(item.compound, item.compound)} ({int(item.experiments)})",
                transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=TEXT, color=INK)
        if item.compound in DRUG_CLASS:
            ax.text(-0.03, yy - 0.5, DRUG_CLASS[item.compound], transform=ax.get_yaxis_transform(), ha="right",
                    va="center", fontsize=NOTE, style="italic", color=SLATE)
        dumbbell(item, yy, 62)
    ax.text(edge, 0.85, "Compound (experiments)", transform=ax.get_yaxis_transform(), ha="left", va="bottom",
            fontsize=NOTE, color=SLATE)
    ax.plot([0, 0], [y - 0.6, -1.0], color=INK, lw=0.9, zorder=1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks(range(-10, 36, 10))
    ax.set_xlabel("Increase in median lifespan versus controls (%)")
    ax.grid(axis="x", color="#ECEAE4", lw=0.7, zorder=0)
    key = [Line2D([], [], ls="none", marker="o", ms=8, mfc=LIGHT, mec=LIGHT, label="Lightest quartile"),
           Line2D([], [], ls="none", marker="o", ms=8, mfc="white", mec=HEAVY, mew=2, label="Heaviest quartile")]
    ax.legend(handles=key, loc="lower right", frameon=True, framealpha=1, edgecolor="#D9D6CF", borderpad=0.5)
    save_figure(fig, output, "Figure_4_compound_quartiles", tight=False)
    names = {"lowest_control_median": "lightest_control_median_days", "lowest_treated_median": "lightest_treated_median_days",
             "lowest_gain_pct": "lightest_increase_percent", "highest_control_median": "heaviest_control_median_days",
             "highest_treated_median": "heaviest_treated_median_days", "highest_gain_pct": "heaviest_increase_percent"}
    panel = table[["compound", "experiments", "unique_mice", *names]].rename(columns=names)
    return {"Fig4": panel.round({"lightest_increase_percent": 2, "heaviest_increase_percent": 2}),
            "Fig4_interaction": pd.DataFrame({"quantity": ["common treatment × weight HR", "lower 95% CI", "upper 95% CI", "P"],
                                              "value": [common.estimate, common.lower_95, common.upper_95, common.p]})}


def figure5(package: Path, output: Path) -> dict:
    """Which weight predicts mortality: by age at weighing, and early-life against six-month weight."""
    age = read_result("developmental_estimates")
    control = read_result("control_weight_gradients")
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
    est = decomp.set_index("label")

    base_style()
    # same proportions as the figure this one replaced; both panels share one HR scale
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.09, 7.09 * 1563 / 2815), gridspec_kw={"width_ratios": [1.3, 1]})
    fig.subplots_adjust(left=0.075, right=0.985, top=0.91, bottom=0.15, wspace=0.32)
    limits, ticks = (0.85, 1.7), np.arange(0.9, 1.65, 0.1)

    ages = [1.4, 6, 12]
    for sex, colour, name in (("m", DEVELOPMENTAL_BLUE, "Males"), ("f", DEVELOPMENTAL_ORANGE, "Females")):
        values = np.array([gradients[(sex, value)] for value in ages])
        ax_a.fill_between(ages, values[:, 1], values[:, 2], color=colour, alpha=0.15, lw=0)
        ax_a.plot(ages, values[:, 0], "-o", color=colour, lw=2, ms=6, mec="white", mew=1)
        ax_a.text(12.6, values[-1, 0], name, color=colour, fontsize=TEXT, fontweight="bold", ha="left", va="center")
    ax_a.axhline(1, color=INK, lw=0.9)
    ax_a.text(16.3, 0.99, "no association", ha="right", va="top", fontsize=NOTE, color=INK)
    ax_a.set_xticks(ages, ["1.4", "6", "12"])
    ax_a.set_xlim(0.2, 16.5)
    ax_a.set_ylim(*limits)
    ax_a.set_yticks(ticks, [f"{t:.1f}" for t in ticks])
    ax_a.set_xlabel("Age at weighing (months)")
    ax_a.set_ylabel("Mortality HR per s.d. heavier")
    ax_a.grid(axis="y", color="#ECEAE4", lw=0.7)
    panel_title(ax_a, "a", "Heavier controls die sooner")

    for x, weight, colour in ((0, "1.4-month", DEVELOPMENTAL_BLUE), (1, "6-month", DEVELOPMENTAL_ORANGE)):
        alone, adjusted = est.loc[f"{weight} alone"], est.loc[f"{weight} adjusted"]
        for row, dx in ((alone, -0.15), (adjusted, 0.15)):
            ax_b.plot([x + dx] * 2, [row.lo, row.hi], color=colour, lw=1.1, alpha=0.6)
        ax_b.annotate("", (x + 0.15, adjusted.HR), (x - 0.15, alone.HR),
                      arrowprops=dict(arrowstyle="-|>", color=colour, lw=1.6, shrinkA=6, shrinkB=6))
        ax_b.scatter([x - 0.15], [alone.HR], s=60, facecolors="white", edgecolors=colour, linewidths=1.8, zorder=3)
        ax_b.scatter([x + 0.15], [adjusted.HR], s=60, c=colour, edgecolors="white", zorder=3)
    ax_b.scatter([], [], s=60, facecolors="white", edgecolors=INK, linewidths=1.8, label="Alone")
    ax_b.scatter([], [], s=60, c=INK, edgecolors="white", label="Both in one model")
    ax_b.legend(frameon=False, loc="upper left", handletextpad=0.3, borderaxespad=0.3)
    ax_b.axhline(1, color=INK, lw=0.9)
    ax_b.set_xticks([0, 1], ["1.4-month\nweight", "6-month\nweight"])
    for tick, colour in zip(ax_b.get_xticklabels(), [DEVELOPMENTAL_BLUE, DEVELOPMENTAL_ORANGE]):
        tick.set_color(colour)
        tick.set_fontweight("bold")
    ax_b.set_xlim(-0.6, 1.6)
    ax_b.set_ylim(*limits)
    ax_b.set_yticks(ticks, [f"{t:.1f}" for t in ticks])
    ax_b.set_ylabel("Mortality HR per s.d. heavier (males)")
    ax_b.grid(axis="y", color="#ECEAE4", lw=0.7)
    panel_title(ax_b, "b", "Adult weight carries the risk")
    save_figure(fig, output, "Figure_5_developmental", dpi=600, tight=False)
    columns = ["mortality_HR_per_SD", "lower_95", "upper_95"]
    by_age = pd.DataFrame([{"sex": "male" if sex == "m" else "female", "age_months": months,
                            **dict(zip(columns, gradients[(sex, float(months))]))}
                           for sex in ("m", "f") for months in ages])
    return {"Fig5a": by_age, "Fig5b": decomp[["label", "HR", "lo", "hi"]].set_axis(["model", *columns], axis=1)}


def extended_data_figure1(package: Path, output: Path) -> dict:
    """All-experiment treatment-by-weight estimates."""
    from matplotlib.lines import Line2D

    interactions = pd.read_csv(DERIVED / "experiment_interactions.csv")
    weight = pd.read_csv(DERIVED / "experiment_weight_effects.csv")
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
    interactions[["cohort", "group", "sex", "int_hr", "int_lo", "int_hi", "int_p", "weight_lowering"]].to_csv(
        package / "data" / "outputs" / "extended_data_figure_1_interactions.csv", index=False, lineterminator="\n")
    observed = interactions.groupby("sex").weight_lowering.sum().to_dict()
    if observed != {"m": 42, "f": 40}:
        raise AssertionError(f"Unexpected weight-lowering counts: {observed}")
    pooled = read_result("primary_interaction").set_index("sex")
    experiments = interactions.groupby("sex").size()

    base_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 5.5), sharex=True)
    panels = zip(axes, ["m", "f"], ["a", "b"],
                 [f"Males ({experiments['m']} experiments)", f"Females ({experiments['f']} experiments)"])
    for ax, sex, letter, title in panels:
        frame = interactions[interactions.sex == sex].sort_values("int_hr").reset_index(drop=True)
        y = np.arange(len(frame))
        colours = np.where(frame.weight_lowering, BLUE, AMBER)
        for yi, lo, hi, colour in zip(y, frame.int_lo, frame.int_hi, colours):
            ax.plot([max(lo, 0.35), min(hi, 2.1)], [yi, yi], color=colour, alpha=0.25, lw=0.8)
        ax.scatter(frame.int_hr.clip(0.35, 2.1), y, c=colours, s=8, alpha=0.8, edgecolors="none")
        ax.axvline(1, color=INK, lw=1.0)
        estimate = pooled.loc[sex]
        ax.errorbar(estimate.HR, -7,
                    xerr=[[estimate.HR - estimate.lo], [estimate.hi - estimate.HR]],
                    fmt="D", ms=6, color=INK, capsize=2.5, lw=1.6)
        ax.text(0.37, -11, f"Pooled HR {estimate.HR:.2f} ({estimate.lo:.2f}–{estimate.hi:.2f})",
                fontsize=TEXT, fontweight="bold", zorder=5,
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"))
        ax.set_ylim(len(frame) + 2, -15)
        ax.set_yticks([])
        ax.set_title(letter, loc="left", fontweight="bold", fontsize=LETTER, pad=18)
        ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 5),
                    textcoords="offset points", ha="left", va="bottom",
                    fontsize=TITLE, color="#333333")
        ax.set_xlabel("Treatment × weight interaction HR per +1 s.d.")
        ax.grid(axis="x", color="#E4E1DA", lw=0.7)
    key = [Line2D([], [], ls="none", marker="o", ms=5, mfc=BLUE, mec="none", label="Weight-lowering experiment"),
           Line2D([], [], ls="none", marker="o", ms=5, mfc=AMBER, mec="none", label="Other experiment"),
           Line2D([], [], ls="none", marker="D", ms=6, mfc=INK, mec=INK, label="Pooled estimate (95% CI)")]
    fig.tight_layout(rect=[0, 0.07, 1, 0.98])
    fig.legend(handles=key, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, frameon=False, fontsize=TEXT,
               handletextpad=0.4, columnspacing=2.0)
    save_figure(fig, output, "Extended_Data_Figure_1_all_experiment_interactions")
    table = interactions[["cohort", "group", "sex", "int_hr", "int_lo", "int_hi", "int_p", "weight_lowering"]].set_axis(
        ["cohort", "intervention_group", "sex", "interaction_HR_per_SD", "lower_95", "upper_95", "P", "weight_lowering"], axis=1)
    sheets = {}
    for letter, sex in (("a", "m"), ("b", "f")):
        sheets[f"EDFig1{letter}"] = table[table.sex.eq(sex)].drop(columns="sex").sort_values("interaction_HR_per_SD")
        sheets[f"EDFig1{letter}_pooled"] = pd.DataFrame({"quantity": ["pooled HR per SD", "lower 95% CI", "upper 95% CI"],
                                                         "value": pooled.loc[sex, ["HR", "lo", "hi"]].astype(float).to_numpy()})
    return sheets


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, default=Path(__file__).absolute().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    package = args.package.absolute()
    output = (args.output or package / "figures").absolute()
    for stem, make in (("Fig1", figure1), ("Fig2", figure2), ("Fig3", figure3), ("Fig4", figure4), ("Fig5", figure5),
                       ("Extended_Data_Fig1", extended_data_figure1)):
        save_source_data(output / "source_data" / f"Source_Data_{stem}.xlsx", make(package, output))
    print(f"Regenerated five manuscript figures and one Extended Data figure, with source data, in {output}")


if __name__ == "__main__":
    main()
