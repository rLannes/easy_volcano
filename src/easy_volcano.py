"""Volcano plots for differential-expression results.

Draws volcano plots (-log10 adjusted p-value vs. log2 fold change) from a
pd.DataFrame (DESeq2-style results table). The single public entry point
is `plot_volcano`.

Points are coloured by default according to significance thresholds
(up / down / one criterion met / not significant), but the styling is fully
overridable through a small selector/override system: any subset of points
-- picked by a boolean predicate, a boolean mask, or a list of index labels
-- can be given arbitrary `ax.scatter` keyword arguments (colour, size,
marker, alpha, ...). User-supplied `highlights` are drawn on a top layer so
they are never hidden behind other points.

Two private helpers support this:
    _resolve(df, selector) -- turn a selector into a boolean mask.
    _full(n, val)          -- build a length-n per-point property array,
                              choosing numeric vs. object dtype so colours
                              and markers survive.

matplotlib is configured for editable text in vector output
(pdf.fonttype = 42) so exported PDFs can be edited in Illustrator/Inkscape.
"""

import matplotlib.pyplot as plt
import math
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, MaxNLocator)
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import offset_copy
import numbers
import numpy as np
import pandas as pd
import sys
import adjustText
import functools

import logging
logger = logging.getLogger(__name__)


_MPL_STYLE = { "font.stretch" : "expanded", 'pdf.fonttype':  42 }
# decorator 
def _with_rc_param(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with mpl.rc_context(_MPL_STYLE):
            return func(*args, **kwargs)
    return wrapper

# sentinel 
_UNSET = object()

def _resolve(df, selector, label_col=None):
    """
    private allow to pass multiple type of selection
    """
    if callable(selector):                 
        result = selector(df)
    elif isinstance(selector, pd.Series):    
        result = selector.reindex(df.index, fill_value=False)
    elif isinstance(selector, (list, set, pd.Index, tuple)): 
        key = df.index if label_col is None else df[label_col]
        result =  key.isin(selector)
    else:
        raise TypeError(selector)
    return np.asarray(result, dtype=bool) 

def _full(n, val):
    """Length-n array holding val, with a sensible dtype."""
    if isinstance(val, (int, float)) and not isinstance(val, bool):
        return np.full(n, float(val))
    arr = np.empty(n, dtype=object)   # object dtype holds strings, RGBA tuples, etc.
    arr[:] = [val] * n
    return arr

_ALIASES = {"size": "s", "color": "c", 'a': "alpha"}
def mpl_norm_symbol(override, alias=_ALIASES):
    """ map potential aliases to matplotlib arg"""
    new = {}
    for k, v in override.items():
        key = alias.get(k, k)
        if key in new:
            raise ValueError(f"'{k}' and '{key}' both given : {override}; please only use one") 
        new[key] = v
    return new


def prep_df(df, x, y, pval_floor, clip_pval=None):
        
        n0 = len(df)
        df = df.copy()
        df = df.dropna(subset=[x, y])
        if df.empty:
            raise ValueError("Data Frame is empty after nan removal")
        
        positive = df.loc[df[y] > 0, y]
        if positive.empty:
            raise ValueError(f"all the p-value ({y!r}) are either negative or zeros; there is a likely issue with your data. check your input is your data already -log10 transformed?")

        if len(df) < n0:
            logger.info("Dropped {} NA rows".format(n0 - len(df)))

        if clip_pval:
            if (clip_pval <= 0) or (clip_pval>=1) :
                raise ValueError(f"clip pval must be in range [0-1]. It was {clip_pval}")
            min_ = clip_pval
        elif pval_floor:
            logger.info("using sys.float_info.min as pvlaue floor")
            min_ = sys.float_info.min
        else:
            logger.info("using min(df[pvalue]) as pvalue floor")
            min_ = min(df[df[y] >0][y])

        logger.info("pvalue floor: {}".format(min_))

        df["neglog10"] = -np.log10(df[y].clip(lower=min_))
        df["x"] = df[x]
        df["y"] = df[y]
        return df

def prep_legend(highlights_default, highlights, **scatter_kwargs):
    legend = {}
    #pop the label()
    for (sel, ov)  in highlights_default:
        label = ov.pop("label", None)
        if label:
            legend[label] = {**scatter_kwargs, **ov}

    for (sel, ov)  in highlights:
        label = ov.pop("label", None)
        zorder = ov.get("zorder", 2)
        if label:
            legend[label] = {"zorder": zorder, **scatter_kwargs, **ov}
    return legend


def _extract_label_gene(highlight, seen, res, df, label_column, x, label_format):
    for (sel, ov) in highlight:
        label_f = ov.pop("label_gene", None)

        if not label_f:
            continue

        mask = _resolve(df, sel, label_col=label_column)
        sub = df[mask]

        sub = sub[~sub.index.isin(seen)]
        if sub.empty:
            continue
        
        seen.update(sub.index)

        x_ = sub[x]
        y_ = sub["neglog10"]
        if label_column is not None:
            text = list(sub[label_column])
        else:
            text = list(sub.index)

        label_f = label_f if isinstance(label_f, dict) else {}
        label_f = _parse_label_format({**label_format, **label_f})

        #this_label_format = {**label_format, **label_f}
        res.append((x_, y_, text, label_f))
    return res


def get_index_point_to_labels(df, highlight, label_format, label_genes_format, x, y, x_thr=1,
                            label_genes=None, label_top=None, label_column=None):

    res = []
    seen = set() # priority: hilight > set > top
    # hilight latest will overdraw/erase the previous one if duplicate index
    _extract_label_gene(highlight[::-1], seen, res, df, label_column, x, label_format)

    for (sel, ov) in highlight:
        label_f = ov.pop("label_gene", None)


    index_to_get = []
    label_genes_format = label_genes_format if isinstance(label_genes_format, dict) else {}
    this_label_format = _parse_label_format({**label_format, **label_genes_format})
    #this_label_format = { , **(label_genes_format or {})}
    if label_genes is not None: 
        if isinstance(label_genes, str):
            label_genes = [label_genes]
        label_genes = list(label_genes)

        if label_column is not None:
            missing = set(label_genes) - set(df[label_column])
        else:
            missing = set(label_genes) - set(df.index)

        if missing:
            logger.warning(f"genes not found: {missing}")

        if label_column is not None :
            index_to_get.extend(df[df[label_column].isin(label_genes)].index)
        else:
            index_to_get.extend(df[df.index.isin(label_genes)].index)

        sub = df[df.index.isin(index_to_get)]
        sub = sub[~sub.index.isin(seen)]
        if not sub.empty:
        
            seen.update(sub.index)

            x_ = sub[x]
            y_ = sub["neglog10"]
            if label_column is not None:
                text = list(sub[label_column])
            else:
                text = list(sub.index)

            res.append((x_, y_, text, this_label_format))


    index_to_get = []
    if label_top:
        x_t = x_thr if x_thr is not None else 0
        index_to_get.extend(list(df[(df[x] < -x_t) ].sort_values(y).head(label_top).index))
        index_to_get.extend(list(df[(df[x] > x_t) ].sort_values(y).head(label_top).index))
        sub = df[df.index.isin(index_to_get)]
        
        sub = sub[~sub.index.isin(seen)]
        if not sub.empty:

            seen.update(sub.index)
            x_ = sub[x]
            y_ = sub["neglog10"]
            if label_column is not None:
                text = list(sub[label_column])
            else:
                text = list(sub.index)
            res.append((x_, y_, text, label_format))



    return res



def get_highlight_default(y_thr, x_thr, no_default, x, y, color_sig, color_ns, color_up, color_down, color_up_sig, color_down_sig):
        highlights_default = []
        if not (y_thr is None or x_thr is None or no_default):
            highlights_default = [
                (lambda d: (d[y] <= y_thr) & (abs(d[x]) < x_thr), { "c": color_sig}),
                (lambda d: (d[y] > y_thr) & (abs(d[x]) < x_thr), { "c": color_ns}),

                (lambda d: (d[y] > y_thr) & (d[x] >= x_thr), { "c": color_up}),
                (lambda d: (d[y] > y_thr) & (d[x] <= -x_thr), { "c": color_down}),

                (lambda d: (d[y] <= y_thr) & (d[x] >= x_thr), { "c": color_up_sig, "label": "Overexpressed"}),
                (lambda d: (d[y] <= y_thr) & (d[x] <= -x_thr), { "c": color_down_sig,  "label": "Underexpressed"})
                ]
        return highlights_default


def get_array_style(df, highlights_default, highlights, style, label_column):
        
        n = len(df)
        arrays = {k: _full(len(df), v) for k, v in style.items()}
        for (selector, override) in highlights_default:
            mask = _resolve(df, selector)
            for key, v in override.items():
                if key not in arrays:                 # property not in base → create it
                    arrays[key] = np.full(n, _UNSET, dtype=object) # baseline sentinel       
                a = arrays[key]
                if a.dtype == object:                 # assign as list so tuples don't broadcast
                    #a[mask] = [v] * int(np.count_nonzero(mask))
                    for i in np.flatnonzero(mask):
                        a[i] = v
                else:
                    a[mask] = v
                    
        #arrays = {k: _full(len(df), v) for k, v in style.items()}
        # apply custom style if any
        drawn_last = np.zeros(len(df), dtype=bool)
        for (selector, override) in highlights:
            mask = _resolve(df, selector, label_col=label_column)
            drawn_last |= mask
            for key, v in override.items():
                if key not in arrays:                 # property not in base → create it
                    arrays[key] = np.full(n, _UNSET, dtype=object) # baseline sentinel
                    #arrays[key] = _full(n, v)         # baseline = its own value; see note
                a = arrays[key]
                if a.dtype == object:                 # assign as list so tuples don't broadcast
                    #a[mask] = [v] * int(np.count_nonzero(mask))
                    for i in np.flatnonzero(mask):
                        a[i] = v
                else:
                    a[mask] = v
        return arrays, drawn_last

CALL_LEVEL = {"marker", "hatch", "zorder", "rasterized"}
# more exist but I don't think those heve their place : cmap, norm, vmin, vmax, 
# may be in future: rasterized, hatch no sure hatch can be applied to scatter plot mauy end up using collection patch?

_LABEL_FORMAT_DEFAULT = {"fontsize": 10, "fontstyle": "normal", "ha": "center", "va": "center"}
_BLACK_LIST_LABEL_FORMAT = ["x", "y", "s", "position", "text"]
_LINES_FORMAT_DEFAULT = {"color": "black", "linestyle": "--"}

def _check_lines_param(lines):
    if not (isinstance(lines, list) or isinstance(lines, dict)):
        raise ValueError("h_lines parameter must be a list of dict or a dict")
    
    lines = lines if isinstance(lines, list) else [lines]
    return lines 

def _parse_label_format(label_format, _BLACK_LIST=_BLACK_LIST_LABEL_FORMAT):

    label_format = {**_LABEL_FORMAT_DEFAULT, **(label_format or {})}
    for e in _BLACK_LIST:

        if e in  label_format:
            raise ValueError(f"label_format cannot contain {e!r} as key")

    return label_format

_LEGEND_PRESETS = {
    "right":  {"bbox_to_anchor": (1.0, 0.5), "fancybox": False, "loc": "center left",  "ncols": 1, "borderaxespad": 0},
    "top":    {"bbox_to_anchor": (0.5, 1.0), "fancybox": False, "loc": "lower center", "borderaxespad": 0},
    "bottom": {"bbox_to_anchor": (0.5, 0.0), "fancybox": False, "loc": "upper center", "borderaxespad": 0},
}

# 72 points to inches
def _px_to_pt(fig, px):
    return px * 72 / fig.dpi

# used for legend placement such a pain
def _gap_above(ax, renderer):
    """Points between the axes top and the top of everything drawn so far."""
    bb = ax.get_tightbbox(renderer)
    return _px_to_pt(ax.figure, bb.y1 - ax.bbox.y1)

def _gap_below(ax, renderer):
    bb = ax.get_tightbbox(renderer)
    return _px_to_pt(ax.figure, ax.bbox.y0 - bb.y0)

@_with_rc_param
def plot_volcano(df,
                 *,
# turn pval_floor to false and clip column y to this value. (for plotting)
                 color_bg = "white",
                 x="log2FoldChange", y = "padj",
                 x_thr=1, y_thr=0.01, 
                 pval_floor=False,
                 clip_pval = None,
                 highlights=None, 

                 label_genes = None,
                 label_genes_format = None, # if set apply to label_genes
                 label_genes_line = False,
                 label_min_arrow_len = 100,
                 label_top = 5, # take int label the top X gene over and under expressed (so *2 labels)
                 label_column = None, # instead of using the index look into this column for label name # label_genes must be from this columns
                 label_format = _LABEL_FORMAT_DEFAULT,
                 label_avoid_self = False,
                
                
                color_sig = "darkgrey",
                color_ns = "lightgrey",
                color_up = "#FDAE61",
                color_up_sig = "#D7191C",
                color_down =  "#ABD9E9",
                color_down_sig = "#2C7BB6",
                marker = ".",
                alpha=0.55,

                show_thresholds_line = True,
                line_format_default = _LINES_FORMAT_DEFAULT,
                h_lines = None, # take dict / list of dict: dict {"y": position, **line_kargs} # (y -> transform to neglog or directly neglog)
                v_lines = None, # take  dict / list of dict:  {"x": position, **line_kargs}

                display_grid=True, 
                no_default=False,
                xlabel='$Log_{2}$ fold change', ylabel='$-Log_{10}$ adj.P',
                size=100, color_base="grey",
                width=6, height=6,
                ax = None,
                rasterize = False, # if set override any other arg and rasterize all
                savefig=None,

                title="",
                title_pad=15,
                title_fontsize=28,
                subtitle="",
                subtitle_pad=15,
                subtitle_fontsize=18,


                show_legend=True,
                legend_loc="bottom",
                legend_kwargs = None,
                legend_nrows = None, # default
                padding_legend=15,

                show_text_summary = True,
                font_size_text_summary=16,
                padding_text_summary = 15,
                **scatter_kwargs):


    """Draw a volcano plot from a differential-expression results table.

    Each gene is a point: log2 fold change on the x-axis, -log10 of the adjusted
    p-value on the y-axis. Points are first coloured by significance category
    (see "Default colours"), then any `highlights` are drawn on top of them.
    Rows with a missing value in the `x` or `y` column are dropped.

    Example
    -------
    >>> fig, ax = plot_volcano(df, y_thr=1e-6, label_top=5,
    ...     highlights=[(["FKBP5", "PER1"], {"color": "purple", "size": 200, "label": "my genes"})])

    Parameters
    ----------
    df : pandas.DataFrame
        Results table with one row per gene (e.g. DESeq2, edgeR or limma output).
        Gene names are read from the index, or from `label_column`.

    Data and thresholds
    ~~~~~~~~~~~~~~~~~~~
    x : str, default "log2FoldChange"
        Column with the log2 fold change (x-axis). edgeR and limma: "logFC".
    y : str, default "padj"
        Column with the adjusted p-value, plotted as -log10 (y-axis).
        edgeR: "FDR"; limma: "adj.P.Val".
    x_thr : float or None, default 1
        Absolute log2 fold-change threshold. Used for the default colours, the
        vertical threshold lines, the summary and `label_top`.
    y_thr : float or None, default 0.01
        Adjusted p-value threshold (a p-value, not its -log10). Used for the
        default colours, the horizontal threshold line and the summary.
        If `x_thr` or `y_thr` is None, the default colouring is turned off.
    clip_pval : float or None, default None
        If set (strictly between 0 and 1), every p-value below it is raised to it,
        so a few extreme genes do not squash the rest of the plot.
        Takes priority over `pval_floor`.
    pval_floor : bool, default False
        Where p-values of exactly 0 are drawn (log10(0) is undefined).
        False: at the smallest non-zero p-value of the table.
        True: at the smallest positive float (about 1e-308), i.e. far above the rest.

    Default colours
    ~~~~~~~~~~~~~~~
    Each gene falls in one category, from its p-value (p, column `y`) and its
    log2 fold change (lfc, column `x`):

    color_up_sig : str, default "#D7191C"
        p <= y_thr and lfc >= x_thr. Shown in the legend as "Overexpressed".
    color_down_sig : str, default "#2C7BB6"
        p <= y_thr and lfc <= -x_thr. Shown in the legend as "Underexpressed".
    color_sig : str, default "darkgrey"
        p <= y_thr but |lfc| < x_thr.
    color_up : str, default "#FDAE61"
        p > y_thr but lfc >= x_thr.
    color_down : str, default "#ABD9E9"
        p > y_thr but lfc <= -x_thr.
    color_ns : str, default "lightgrey"
        p > y_thr and |lfc| < x_thr.
    no_default : bool, default False
        Turn the default colouring off: every point gets `color_base`.
        Combine with `highlights` for a "spotlight" plot.

    Points
    ~~~~~~
    color_base : str, default "grey"
        Colour of points that no category or highlight colours.
    size : float, default 100
        Point size (matplotlib's `s`, an area in points squared).
    marker : str, default "."
        Point shape: any matplotlib marker ("o", "D", "^", ...).
    alpha : float, default 0.55
        Point transparency, from 0 (invisible) to 1 (opaque).
    rasterize : bool, default False
        When saving to PDF or SVG, draw the points as an image while text, axes
        and lines stay vector: much smaller files that open quickly in
        Illustrator. Overrides any "rasterized" set in `highlights`.
    **scatter_kwargs
        Any other keyword is passed to matplotlib's `ax.scatter` for every point,
        e.g. edgecolors="none" or linewidths=0.

    Highlights
    ~~~~~~~~~~
    highlights : list of (selector, style) pairs, default None
        Restyle chosen genes. Highlighted genes are drawn on top of all other
        points. Pairs are applied in order: when a gene matches several pairs,
        the later pair wins for the options it sets.

        selector (which genes), one of:
          - a list, tuple or set of names, matched against the index, or against
            `label_column` if set. A single gene is still a list: ["FKBP5"].
          - a function of the table, e.g. lambda d: d["baseMean"] < 100.
            `d` is `df` without the NA rows, so any column can be used.
          - a boolean pandas Series indexed like `df`.

        style (how to draw them): a dict of `ax.scatter` options ("color", "size",
        "alpha", "marker", "edgecolors", "linewidths", "zorder", ...). "color",
        "size" and "a" are short for matplotlib's "c", "s" and "alpha".
        Two extra keys:
          - "label": str, adds the group to the legend.
          - "label_gene": True to write the gene names on the plot, or a dict of
            text options (like `label_format`) to also style those names.

        Example: [(lambda d: d.index.str.startswith("COL"),
                   {"color": "green", "size": 150, "label": "collagens"})]

    Gene labels
    ~~~~~~~~~~~
    Labels are positioned with adjustText to limit overlaps. Each gene is labelled
    once, with priority: highlight "label_gene" > `label_genes` > `label_top`.

    label_top : int or None, default 5
        Label the N genes with the smallest p-value on each side (lfc > x_thr and
        lfc < -x_thr), so up to 2 * N labels. None: no automatic labels.
    label_genes : str or list of str, default None
        Genes to label, matched against the index or `label_column`.
        Names that are not found are reported in a warning.
    label_column : str or None, default None
        Column holding the gene names to use, in labels and in name selectors,
        instead of the index (e.g. "gene_name" when the index holds Ensembl IDs).
    label_format : dict, default {"fontsize": 10, "fontstyle": "normal", "ha": "center", "va": "center"}
        Text options for all labels (matplotlib's `ax.text`), merged on the default.
        The keys "x", "y", "s", "position" and "text" are reserved.
    label_genes_format : dict or None, default None
        Text options applied on top of `label_format` to the `label_genes` only,
        e.g. {"fontweight": "bold", "color": "purple"}.
    label_genes_line : bool, default False
        Draw a line between each label and its point when the label was moved
        further than `label_min_arrow_len` (applies to all labels).
    label_min_arrow_len : float, default 100
        Shortest line drawn by `label_genes_line`, in pixels. With the default,
        most lines are too short to be drawn: use about 5 to see them.
    label_avoid_self : bool, default False
        Also push each label away from its own point.

    Lines and grid
    ~~~~~~~~~~~~~~
    show_thresholds_line : bool, default True
        Draw the threshold lines at `y_thr` and at -`x_thr` and +`x_thr`.
    line_format_default : dict, default {"color": "black", "linestyle": "--"}
        Style of the threshold lines (matplotlib's `ax.axhline` / `ax.axvline`
        options), merged on the default.
    h_lines : dict or list of dict, default None
        Extra horizontal lines. Each dict needs "y", a p-value (drawn at
        -log10(y)), plus any line option, e.g. {"y": 1e-20, "color": "red"}.
        Black dashed by default.
    v_lines : dict or list of dict, default None
        Extra vertical lines. Each dict needs "x", a log2 fold change, plus any
        line option, e.g. {"x": 2, "linestyle": ":"}. Black dashed by default.
    display_grid : bool, default True
        Draw major and minor grid lines.

    Title, legend and summary
    ~~~~~~~~~~~~~~~~~~~~~~~~~
    The elements around the plot are stacked so they never overlap. Above the
    plot, from the axes up: legend (if legend_loc="top"), subtitle, title. Below
    the plot, from the x-axis label down: legend (if legend_loc="bottom"),
    summary. Paddings are gaps in points, so they do not depend on figure size.

    title : str, default ""
        Plot title. Not drawn if empty.
    title_fontsize : float, default 28
    title_pad : float, default 15
        Gap (points) between the title and the element below it.
    subtitle : str, default ""
        Text between the title and the plot. Not drawn if empty.
    subtitle_fontsize : float, default 18
    subtitle_pad : float, default 15
        Gap (points) between the subtitle and the element below it.
    xlabel : str, default "$Log_{2}$ fold change"
    ylabel : str, default "$-Log_{10}$ adj.P"
        Axis labels. Text between $ signs is rendered as math.
    show_legend : bool, default True
        Draw the legend. False is the same as legend_loc=None.
    legend_loc : {"bottom", "top", "right"} or None, default "bottom"
        Legend position. None: no legend. Entries are "Overexpressed",
        "Underexpressed" and every highlight with a "label"; with no entries,
        no legend is drawn.
    legend_nrows : int or None, default None
        Number of legend rows. None: one row for "top" and "bottom", one column
        for "right". An explicit "ncols" in `legend_kwargs` wins.
    legend_kwargs : dict or None, default None
        Options passed to matplotlib's `ax.legend`, applied on top of the preset,
        e.g. {"fontsize": 12, "frameon": False, "title": "groups"}.
    padding_legend : float, default 15
        Gap (points) between the legend and the plot ("top") or the x-axis
        label ("bottom"). Not used for "right".
    show_text_summary : bool, default True
        Write "up N; down N; ns N" under the plot: genes with p <= y_thr and
        lfc >= x_thr, genes with p <= y_thr and lfc <= -x_thr, and all others.
        Needs both `x_thr` and `y_thr`.
    font_size_text_summary : float, default 16
    padding_text_summary : float, default 15
        Gap (points) between the summary and the element above it.

    Figure and output
    ~~~~~~~~~~~~~~~~~
    ax : matplotlib Axes or None, default None
        Draw into this existing Axes, e.g. one panel of plt.subplots(); `width`
        and `height` are then ignored. None: create a new figure.
    width : float, default 6
        Figure width in inches.
    height : float, default 6
        Figure height in inches.
        Tick and axis-label font sizes follow the size of the plot.
    color_bg : str, default "white"
        Background colour of the plot area.
    savefig : str or path, default None
        Save the figure to this file (format from the extension: .png, .pdf,
        .svg, ...), cropped to its content. Text stays editable in PDF files.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
        Keep customising them with matplotlib (annotations, shading, limits, ...).

    Raises
    ------
    ValueError
        If no row is left after dropping missing values, if no p-value is above 0
        (is `y` already -log10 transformed?), if `clip_pval` is not between 0 and 1,
        or if `legend_loc` is not one of the presets.
    """

    df = prep_df(df, x, y, pval_floor, clip_pval)


    label_format = _parse_label_format(label_format=label_format)
    if not highlights:
        highlights = []

    highlights = [(sel, mpl_norm_symbol(over)) for (sel, over) in highlights]
    scatter_kwargs = mpl_norm_symbol(scatter_kwargs)

    if label_genes is not None:
        if isinstance(label_genes, str):
            label_genes = [label_genes]
        label_genes = list(label_genes)


    
    highlights_default = get_highlight_default(y_thr, x_thr, no_default, x, y, color_sig, color_ns, color_up, color_down, color_up_sig, color_down_sig)
    legend = prep_legend(highlights_default, highlights, marker=marker, **scatter_kwargs)
    # prep label (now wotk with  to [(set(label), format_dict)])
    labels_hold = []
    

    style = {"c": color_base, "zorder": 1, "s": size,
              "alpha" :alpha, "marker": marker, **scatter_kwargs}


    # may tweak this part
    labels_container = get_index_point_to_labels(df, highlight=highlights,
                                                    label_format=label_format,
                                                    label_genes_format=label_genes_format,
                                                    x_thr=x_thr,
                                                    x=x, y=y,
                                                    label_genes=label_genes,
                                                    label_top=label_top, 
                                                    label_column=label_column)
    
    (arrays, drawn_last) =  get_array_style(df, highlights_default, highlights, style, label_column)
    # mpl is annoying those cn only be passed per call not per array
    # use this so modified point are plotted on top.
    
    unset = {k : np.array([e is _UNSET for e in val])
             for k, val in arrays.items() if val.dtype == object}
    unset = {k: v for k, v in unset.items() if v.any()}

    sub = dict((k, v)  for (k, v) in arrays.items() if k in CALL_LEVEL)
    sub.update({f"_unset_{k}": v for k, v in unset.items()})
    sub["_top_layer"] = drawn_last
    group_df = pd.DataFrame(sub)
    indices_same_call_level = group_df.groupby(list(sub), sort=False, dropna=False).indices

    # start plotting
    if ax is None:
        fig , ax = plt.subplots(figsize=(width, height))        
    else:
        fig = ax.figure

    ax.set_facecolor(color_bg)

    for index in sorted(indices_same_call_level.values(), key=lambda p: drawn_last[p[0]]):
        first = index[0]
        keep = lambda k: arrays[k][first] is not _UNSET
        per_call = {k: arrays[k][index[0]] for k in CALL_LEVEL if k in arrays and keep(k)}
        per_array = {k: v[index] for (k,v) in arrays.items() if k not in CALL_LEVEL and keep(k)}
        if rasterize:
            per_call["rasterized"] = True
        ax.scatter(df[x].values[index], df["neglog10"].values[index],
                   **per_array, **per_call)


    texts_ = []
    label_box_container = []
    for (X, Y, text, fmt) in labels_container:
        for xi, yi, s in zip(X, Y , text):
            box = fmt.pop("boxes_label", None)
            t = ax.text(x=xi, y=yi, s=s,  **fmt)
            #if box is not None:
            #    label_box_container.append()
            texts_.append(t)


    # add boxed_label option
    if texts_:
        if label_genes_line:
            adjustText.adjust_text(texts_, expand=(1.1, 1.2), ax=ax,
                                   min_arrow_len=label_min_arrow_len,
                                   avoid_self=label_avoid_self,
                                   arrowprops=dict(arrowstyle='-', color='black'))
        else:
            adjustText.adjust_text(texts_, expand=(1.1, 1.2), ax=ax,
                                   min_arrow_len=label_min_arrow_len,
                                   avoid_self=label_avoid_self)

            #boxes_label
    
    if display_grid:
        
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 2, 5, 10]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))   # 2 minors per major
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.grid(True, which="major", linewidth=0.6, alpha=0.5)

        ax.grid(True, which="minor", linewidth=0.3, alpha=0.3)

    # may refactor line stuff in three small f()
    line_format_default = {**_LINES_FORMAT_DEFAULT, **line_format_default}
    if show_thresholds_line:
        line_format_default = line_format_default if isinstance(line_format_default, dict) else {}
        
        
        if y_thr is not None:
            h_line = -math.log(y_thr, 10)
            ax.axhline(h_line, **line_format_default)

        if x_thr is not None:
            v_line = x_thr
            ax.axvline(v_line, **line_format_default)
            ax.axvline(-v_line, **line_format_default)

    if h_lines:

        h_lines = _check_lines_param(h_lines)
        for h_line in h_lines:

            if not isinstance(h_line, dict):
                raise ValueError("h_lines parameter must be a list of dict or a dict")

            copy_fmt = dict(h_line)
            y_p = copy_fmt.pop("y", None)
            # neglog =  copy_line.pop("neglog", False) not sure about that option; to niche to be usefull
            if (y_p is None) or not (isinstance(y_p, numbers.Real) or (isinstance(y_p, bool))):
                raise ValueError("h_lines dict(s) must contain a 'y' key with a numerical value (pvalue, plot position)")
            if y_p <= 0:
                raise ValueError("'y' cannot be equal to zeros or negative, undefined log10")
            
            copy_fmt = {**_LINES_FORMAT_DEFAULT, **copy_fmt}
            ax.axhline(-math.log10(y_p), **copy_fmt)

    if v_lines:

        v_lines = _check_lines_param(v_lines)
        for v_line in v_lines:

            if not isinstance(v_line, dict):
                raise ValueError("v_lines parameter must be a list of dict or a dict")
            
            copy_fmt = dict(v_line)
            x_p = copy_fmt.pop("x", None)

            if (x_p is None) or not (isinstance(x_p, numbers.Real) or (isinstance(x_p, bool))):
                raise ValueError("v_lines dict(s) must contain a 'x' key with a numerical value (pvalue, plot position)")
            
            copy_fmt = {**_LINES_FORMAT_DEFAULT, **copy_fmt}
            ax.axvline(x_p, **copy_fmt)

    #xlim = ax.get_xlim()
    #ylim = ax.get_ylim()

    w, h = ax.bbox.width / fig.dpi, ax.bbox.height / fig.dpi    # axes size in inches
    scale = min(w / 6.5, h / 3.5)    
    ax.tick_params(axis='y', labelsize=20 * scale)
    ax.tick_params(axis='x', labelsize=20 * scale)
    ax.margins(x=0.1, tight=True)
    ax.margins(y=0.1, tight=True)  
    y_lim = ax.get_ylim()
    ax.set_ylim((-2, y_lim[1]))

    ax.spines[['right', 'top']].set_visible(False)


    ax.set_ylabel(ylabel, fontsize=24 * scale)
    ax.set_xlabel(xlabel , fontsize=24 * scale)


    # I use a tight box approach to place the titel and extar text + legend
    # small performance hit but should be the cleaner
    legend_elements = None
    legend_args = None
    # make legend handle:
    if show_legend and legend_loc is not None:
        if legend_loc not in _LEGEND_PRESETS:
            raise ValueError(f"legend must be one of {list(_LEGEND_PRESETS)} or None, got {legend_loc!r}")
        legend_args = {**_LEGEND_PRESETS[legend_loc], **(legend_kwargs or {})}

        legend_elements = []
        for k, v in legend.items():
            c = v.get("c")
            if not c: 
                c = v.get("color")
            m = v.get("marker", ".")
            m = "o" if m == "." else m 
            a = v.get("alpha", alpha)
            legend_elements.append(Line2D([0], [0],  marker=m, #alpha = a,
                                        markerfacecolor=c, color='w',
                                          label=k, markersize=10))

        n = len(legend_elements)
        if legend_nrows is not None:
            if not isinstance(legend_nrows, int) or legend_nrows < 1:
                raise ValueError(f"legend_nrows must be a positive int or None, got {legend_nrows!r}")
            ncols = math.ceil(n / legend_nrows)
        elif legend_loc in ("top", "bottom"):
            ncols = n            # one row
        else:
            ncols = 1            # one column (right)
        legend_args.setdefault("ncols", ncols)


    renderer = fig.canvas.get_renderer()

    if legend_loc == "bottom" and show_legend:
        
        bottom = offset_copy(ax.transAxes,
                              fig=fig,
                               y=-(_gap_below(ax, renderer) + padding_legend),
                              units="points")
        if legend_elements:
            ax.legend(handles=legend_elements, bbox_transform=bottom, **legend_args)


    if show_text_summary:
        bottom = offset_copy(ax.transAxes,
                              fig=fig,
                               y=-(_gap_below(ax, renderer) + padding_text_summary),
                              units="points")

    # plot a quick summary of the data
        txt = ''
        if x_thr is not None and y_thr is not None:
            up = df[(df[y] <= y_thr) & (df[x] >= x_thr)].shape[0]
            down =  df[(df[y] <= y_thr) & (df[x] <= -x_thr)].shape[0]
            ns = df.shape[0] - (up + down)
            txt = "up {}; down {}; ns {}".format(up, down, ns )
        if txt:
            ax.figure.text(0.5, 0, txt, wrap=True,
                            va="top",
                            horizontalalignment='center',
                            fontsize=font_size_text_summary,
                            transform=bottom)

            

    #borderaxespad=0

    # add legend, subtitle, title
    if legend_loc == "top" and show_legend:
        up = offset_copy(ax.transAxes, fig=fig,
                        y=_gap_above(ax, renderer) + padding_legend,
                            units="points")
        if legend_elements:
            ax.legend(handles=legend_elements,
                   bbox_transform=up,
                   #borderaxespad=0,
                   **legend_args)


    if subtitle:
        up = offset_copy(ax.transAxes, fig=fig,
                                y=_gap_above(ax, renderer) + subtitle_pad,
                                    units="points")
        ax.text(0.5, 1.0, subtitle, transform=up, ha="center",
                 va="bottom", fontsize=subtitle_fontsize)


    if title:
        ax.set_title(title, fontsize=title_fontsize,
                      pad=_gap_above(ax, renderer) + title_pad)

    if legend_loc == "right" and show_legend:
        if legend_elements:
            ax.legend(handles=legend_elements,
                           **legend_args)
        
    plt.tight_layout()
    if savefig:
        fig.savefig(savefig, bbox_inches="tight")

    return (fig, ax)







## add boxed name 
"""
t1 = axs.text(0.4,0.6, 'Hello world line 1', ha='center', color='red', weight='bold', transform=axs.transAxes)
t2 = axs.text(0.5,0.5, 'Hello world line 2', ha='center', color='green', weight='bold', transform=axs.transAxes)
t3 = axs.text(0.6,0.4, 'Hello world line 3', ha='center', color='blue', weight='bold', transform=axs.transAxes)

fig.canvas.draw()

textobjs = [t1,t2,t3]

xmin = min([t.get_window_extent().xmin for t in textobjs])
xmax = max([t.get_window_extent().xmax for t in textobjs])
ymin = min([t.get_window_extent().ymin for t in textobjs])
ymax = max([t.get_window_extent().ymax for t in textobjs])

xmin, ymin = fig.transFigure.inverted().transform((xmin, ymin))
xmax, ymax = fig.transFigure.inverted().transform((xmax, ymax))

rect = patches.Rectangle((xmin,ymin),xmax-xmin,ymax-ymin, facecolor='grey', alpha=0.2, transform=fig.transFigure)

axs.add_patch(rect)

"""