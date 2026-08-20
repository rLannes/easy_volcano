"""Volcano plots for differential-expression results.

Draws volcano plots (-log10 adjusted p-value vs. log2 fold change) from a
DESeq2-style results table stored as a TSV. The single public entry point
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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, MaxNLocator)
import matplotlib as mpl
import numpy as np
import pandas as pd
import sys
import adjustText

mpl.rcParams["font.stretch"] = "expanded"
plt.rcParams['pdf.fonttype'] = 42

def _resolve(df, selector):
    if callable(selector):                 # lambda d: d["value"] > 3
        result = selector(df)
    elif isinstance(selector, pd.Series):    # a precomputed boolean mask
        result = selector.reindex(df.index, fill_value=False)
    elif isinstance(selector, (list, set, pd.Index, tuple)):  # list of index labels
        result =  df.index.isin(selector)
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
    
def plot_volcano(df, alpha_value=1, width=12, height=7.5,
                 p_float_min=False, color_bg = "white",
                 display_gene_name = set(), display_gene_name_line=False,
                 x="log2FoldChange", y = "padj",
                 x_thr=1, y_thr=0.01, 
                 display_grid=True, 
                 highlights={}, 
                 color_ns="grey",
                 color_up="#F4AB5C",
                 color_down="#ACC8BE",
                 color_one_crit="lightgrey",
                 no_default=False,
                 size=5, color="grey",
                 output=None, title="Title", return_fig=False, **defaults):

    """Draw a volcano plot from a differential-expression results table.

    Reads a tab-separated results file, computes -log10 of the adjusted
    p-value, and scatters it against the log2 fold change. Points are first
    coloured by the default significance categories, then any user-supplied
    `highlights` are overlaid on top. Rows containing any NA are dropped
    before plotting.

    
        plot_volcano(df,
             y="padj", x="log2FoldChange",
             highlights=[(["kl-3", "kl-5", "kl-2"], {"c": "red", "s":30}),],
             display_grid=True, display_gene_name=["kl-3", "kl-5", "kl-2"])#

    Parameters
    ----------
    df : pd.DataFrame
        expected to have no NA
    alpha_value : float, default 1
        Point transparency passed to scatter.
    width, height : float, default 12, 7.5
        Figure size in inches.
    p_float_min : bool, default False
        Floor used to avoid log10(0). If True, clip p-values at the smallest
        positive float (sys.float_info.min); if False, clip at the smallest
        positive value observed in column `y`.
    color_bg : str, default "white"
        Figure background colour.
    display_gene_name : set of str, default empty
        Index labels to annotate. Labels are repositioned with adjustText to
        reduce overlap.
    display_gene_name_line : bool
        add a line from gene label to gene
    x : str, default "log2FoldChange"
        Column holding the log2 fold change (x-axis).
    y : str, default "padj"
        Column holding the adjusted p-value (converted to -log10 for the
        y-axis).
    x_thr : float, default 1
        Absolute log2 fold-change threshold used by the default colouring.
    y_thr : float, default 0.01
        Adjusted p-value threshold used by the default colouring (compared on
        the raw p-value, not on -log10).
    display_grid : bool, default True
        Whether to draw major and minor grid lines.
    highlights : list of (selector, override) tuples, default empty
        Custom styling layer, drawn on top of everything else. `selector` is
        anything accepted by _resolve (a predicate `lambda d: ...`, a boolean
        Series, or a list of index labels); `override` is a dict of scatter
        keyword arguments, e.g. {"c": "red", "s": 20}, applied to the
    color_ns : str, default "grey"
        Colour for points meeting neither criterion (not significant and
        within the fold-change threshold).
    color_up : str, default "#F4AB5C"
        Colour for up-regulated points (y < y_thr and x > x_thr).
    color_down : str, default "#ACC8BE"
        Colour for down-regulated points (y < y_thr and x < -x_thr).
    color_one_crit : str, default "lightgrey"
        Colour for points meeting exactly one of the two criteria.
    no_default : bool, default False
        If True (or if `x_thr` or `y_thr` is None), skip the default
        significance colouring; all points keep the base `color`.
    size : float, default 5
        Base point size.
    color : str, default "grey"
        Base point colour, used before any override.
    output : str or path, optional
        If given, save the figure to this path with bbox_inches="tight".
    title : str, default "Title"
        Plot title.
    return_fig : bool, default False
        return (fig, ax)
    **defaults
        Extra keyword arguments forwarded to scatter for every point
        (e.g. marker, linewidths).

    Returns
    -------
    None
        The figure is drawn and optionally saved and or returned.
    """

    df = df.dropna()
    df = df.copy()

    highlights_default = {}
    if not (y_thr is None or x_thr is None or no_default):
        highlights_default= [(lambda d: (d[y] >=y_thr) & (abs(d[x])<=x_thr), { "c": color_ns}),
             (lambda d: (d[y] >=y_thr) ^ (abs(d[x])<=x_thr), { "c": color_one_crit}),
             (lambda d: (d[y] < y_thr) & (d[x] > x_thr), { "c": color_up}),
             (lambda d: (d[y] < y_thr) & (d[x] < -x_thr), { "c": color_down})]



    
    fig , ax = plt.subplots(figsize=(width, height))
    fig.patch.set_facecolor(color_bg)
    
    
    if p_float_min:
        min_ = sys.float_info.min
    else:
        min_ = min(df[df[y] >0][y])
        
    df["neglog10"] = -np.log10(df[y].clip(lower=min_))


    style = {"c": color, "s": size, "alpha" :alpha_value, **defaults}

    n = len(df)
    arrays = {k: _full(len(df), v) for k, v in style.items()}
    for (selector, override) in highlights_default:
        mask = _resolve(df, selector)
        for key, v in override.items():
            if key not in arrays:                 # property not in base → create it
                arrays[key] = _full(n, v)         # baseline = its own value; see note
            a = arrays[key]
            if a.dtype == object:                 # assign as list so tuples don't broadcast
                a[mask] = [v] * int(np.count_nonzero(mask))
            else:
                a[mask] = v
                
    #arrays = {k: _full(len(df), v) for k, v in style.items()}
    drawn_last = np.zeros(len(df), dtype=bool)
    for (selector, override) in highlights:
        mask = _resolve(df, selector)
        drawn_last |= mask
        for key, v in override.items():
            if key not in arrays:                 # property not in base → create it
                arrays[key] = _full(n, v)         # baseline = its own value; see note
            a = arrays[key]
            if a.dtype == object:                 # assign as list so tuples don't broadcast
                a[mask] = [v] * int(np.count_nonzero(mask))
            else:
                a[mask] = v

    base = ~drawn_last
    
    ax.scatter(x = df[x].values[base],
               y = df["neglog10"].values[base],
               zorder=1,
               **{k: v[base] for k, v in arrays.items()})

    if drawn_last.any():
        ax.scatter(x = df[x].values[drawn_last],
               y = df["neglog10"].values[drawn_last],
               zorder=3,
               **{k: v[drawn_last] for k, v in arrays.items()})

    
    texts_ = []
    if display_gene_name: 
        sub = df[df.index.isin(display_gene_name)]
        text = list(sub.index)
        x_axis= x 
        X = sub.apply(lambda x: x[x_axis], axis=1)
        Y = sub.apply(lambda x: x["neglog10"], axis=1)
        obj = []
        for i, v in enumerate(text):
            texts_.append(ax.text(x=X.iloc[i], y=Y.iloc[i], s=v, zorder=15, fontsize=10))
        if display_gene_name_line:
            adjustText.adjust_text(texts_, expand=(1, 1),  arrowprops=dict(arrowstyle='-', color='black'))
        else:
            adjustText.adjust_text(texts_, expand=(1, 1))
    
    if display_grid:
        
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 2, 5, 10]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))   # 2 minors per major
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.grid(True, which="major", linewidth=0.6, alpha=0.5)

        ax.grid(True, which="minor", linewidth=0.3, alpha=0.3)

    if y_thr is None:
        h_line = -math.log(y_thr, 10)
        ax.axhline(h_line, color="black", linestyle="--")

    if x_thr is None:
        v_line = x_thr
        ax.axvline(v_line, color="black", linestyle="--")
        ax.axvline(-v_line, color="black", linestyle="--")
        
    ax.tick_params(axis='y', labelsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.margins(x=0.02, tight=True)
    ax.margins(y=0.02, tight=True)  
    
    ax.set_title(title, fontsize=28)
    ax.spines[['right', 'top']].set_visible(False)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()

    ax.set_ylabel('$-Log_{10}$ adj.P', fontsize=24)
    ax.set_xlabel('$Log_{2}$ fold change', fontsize=24)
    
    if output:
        fig.savefig(output, bbox_inches="tight")
    if return_fig:
        return fig, ax


