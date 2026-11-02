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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, MaxNLocator)
from matplotlib.lines import Line2D
import matplotlib as mpl
import numpy as np
import pandas as pd
import sys
import adjustText

mpl.rcParams["font.stretch"] = "expanded"
plt.rcParams['pdf.fonttype'] = 42

def _resolve(df, selector):
    """
    private allow to pass multiple type of selection
    """
    if callable(selector):                 
        result = selector(df)
    elif isinstance(selector, pd.Series):    
        result = selector.reindex(df.index, fill_value=False)
    elif isinstance(selector, (list, set, pd.Index, tuple)): 
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

CALL_LEVEL = {"marker", "hatch", "zorder"}
# more exist but I don't think those heve thei place : cmap, norm, vmin, vmax, 
# may be in future: rasterized, hatch no sure hatch can be applied to scatter plot mauy end up using collection patch?


def plot_volcano(df,
                 pval_floor=False,
                 color_bg = "white",

                 label_genes = None,
                 label_genes_line = False,
                 label_top = None, # take int label the top X gene over and under expressed (so *2 labels)
                 x="log2FoldChange", y = "padj",
                 x_thr=1, y_thr=0.01, 
                 display_grid=True, 
                 highlights=None, 
                
                color_sig = "#F7DDAA",
                color_ns = "lightgrey",
                color_up = "#FDAE61",
                color_up_sig = "#D7191C",
                color_down =  "#ABD9E9",
                color_down_sig = "#2C7BB6",
                alpha=1,

                 show_thresholds_line=True,
                 no_default=False,
                 xlabel='$Log_{2}$ fold change', ylabel='$-Log_{10}$ adj.P',
                 size=5, color_base="grey",
                 width=12, height=7.5,
                 custom_ax = None,
                 savefig=None, title="Title",
                return_fig=False, legend_bbox=(1.01, 0.5), **scatter_kwargs):


    """Draw a volcano plot from a differential-expression results table.

    Reads a tab-separated results file, computes -log10 of the adjusted
    p-value, and scatters it against the log2 fold change. Points are first
    coloured by the default significance categories, then any user-supplied
    `highlights` are overlaid on top. 
    Rows containing any NA in the required column (lfc and padj) are dropped
    before plotting.

    
        plot_volcano(df,
             y="padj", x="log2FoldChange",
             highlights=[(["kl-3", "kl-5", "kl-2"], {"color": "red", "size":30}),],
             display_grid=True, display_gene_name=["kl-3", "kl-5", "kl-2"])#

    Parameters
    ----------
    df : pd.DataFrame
        Rows containing any NA in the required column (lfc and padj) are dropped
            before plotting.
    alpha : float, default 1
        Point transparency passed to scatter.
    width, height : float, default 12, 7.5
        Figure size in inches.
    pval_floor : bool, default False
        Floor used to avoid log10(0). If True, clip p-values at the smallest
        positive float (sys.float_info.min); if False, clip at the smallest
        positive value observed in column `y`.
    color_bg : str, default "white"
        Figure background colour.
    label_genes : set of str, default empty
        Index labels to annotate. Labels are repositioned with adjustText to
        reduce overlap.
    label_genes_line : bool
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
    highlights : list of (selector, override) tuples, default None
        Custom styling layer, drawn on top of everything else. `selector` is
        anything accepted by _resolve (a predicate `lambda d: ...`, a boolean
        Series, or a list of index labels); `override` is a dict of scatter
        keyword arguments, e.g. {"color": "red", "size": 20} 

        d[x] and d[y] refer to the "x" argument and the "y" argument (usualy lfc and pvalue).
        you can also use any coluns present in the dataframe df.
        [(lambda d: (d[y] >= 2) ^ (abs(d[x]) <= 1), { "color": "blue"})]

        USing a label key work in the override, will add it to the legends

    color_sig : str, default "#F7DDAA"
        gene that pass the  p-value threshold but not the log folchange threshold.
        (y > y_thr AND abs(x) > x_thr)
    
    color_ns : str, default "lightgrey"
        gene that do not pass any threshold.
        (y < y_thr AND abs(x) < x_thr)

    color_up : str, default "#FDAE61"
        Colour for up-regulated points (y < y_thr and abs(x) > x_thr).

    color_up_sig : str, default "#D77191C"
        Colour for significantly up-regulated points (y > y_thr and abs(x) > x_thr).

    color_down : str, default "#ABD9E9"
        Colour for down-regulated points (y < y_thr and abs(x) < -x_thr).
    
    color_down_sig : str, default "#2C7BB6"
        Colour for significantly down-regulated points (y > y_thr and abs(x) < -x_thr).


    no_default : bool, default False
        If True (or if `x_thr` or `y_thr` is None), skip the default
        significance colouring; all points keep the base `color`.
    size : float, default 5
        Base point size.
    marker: string '.'
     default marker 
    color_base : str, default "grey"
        Base point colour, used before any override. if no default and no highlight all the point will be this color.
    savefig : str or path, optional
        If given, save the figure to this path with bbox_inches="tight".
    title : str, default "Title"
        Plot title.
    return_fig : bool, default False
        return (fig, ax)
    xlabel: str default '$Log_{2}$ fold change',
    ylabel: str '$-Log_{10}$ adj.P'
    custom_ax: if set function will not create a figure but plot in this axes directly not compatible with:
      return fig and savefig. 
    **scatter_kwargs
        Extra keyword arguments forwarded to scatter for every point
        (e.g. marker, linewidths).

    Returns
    -------
    None
        The figure is drawn
        optional:
            if return_fig == True  return the (fig, ax)
        if savefig is not None save the figure
    """


    if custom_ax is not None and (return_fig or savefig is not None):
        raise ValueError("custom_ax argument is not compatible with return_fig or savefig argument")
    
    df = df.dropna(subset=[x, y])
    df = df.copy()

    if not highlights:
        highlights = []

    highlights = [(sel, mpl_norm_symbol(over)) for (sel, over) in highlights]
    scatter_kwargs = mpl_norm_symbol(scatter_kwargs)

    if not label_genes:
        label_genes = None

    highlights_default = []
    if not (y_thr is None or x_thr is None or no_default):
        highlights_default = [
             (lambda d: (d[y] <= y_thr) & (abs(d[x]) <= x_thr), { "c": color_sig, "s": size, "alpha": alpha, "marker": ".", "zorder": 1}),
             (lambda d: (d[y] >= y_thr) ^ (abs(d[x]) <= x_thr), { "c": color_ns, "s": size, "alpha": alpha, "marker": ".", "zorder": 1}),

             (lambda d: (d[y] > y_thr) & (d[x] > x_thr), { "c": color_up, "s": size, "alpha": alpha, "marker": ".", "zorder": 1}),
             (lambda d: (d[y] > y_thr) & (d[x] < -x_thr), { "c": color_down, "s": size, "alpha": alpha, "marker": ".", "zorder": 1}),

             (lambda d: (d[y] <= y_thr) & (d[x] > x_thr), { "c": color_up_sig, "s": size, "alpha": alpha, "marker": ".", "label": "Overexpressed", "zorder": 1}),
             (lambda d: (d[y] <= y_thr) & (d[x] < -x_thr), { "c": color_down_sig, "s": size, "alpha": alpha, "marker": ".", "label": "Underexpressed", "zorder": 1})
             ]


    if custom_ax:
        ax = custom_ax
    else:
        fig , ax = plt.subplots(figsize=(width, height))
        fig.patch.set_facecolor(color_bg)
    
    if pval_floor:
        min_ = sys.float_info.min
    else:
        min_ = min(df[df[y] >0][y])
        
    df["neglog10"] = -np.log10(df[y].clip(lower=min_))
    df["x"] = df[x]
    df["y"] = df[y]

    legend = {}
    #pop the label()
    for (sel, ov)  in highlights_default:
        label = ov.pop("label", None)
        if label:
            legend[label] = {**scatter_kwargs, **ov}

    for (sel, ov)  in highlights:
        label = ov.pop("label", None)
        zorder = ov.pop("zorder", 2)
        if label:
            legend[label] = {"zorder": zorder, **scatter_kwargs, **ov}
            

    style = {"c": color_base, "zorder": 1, "s": size, "alpha" :alpha, "marker": '.', **scatter_kwargs}
    # apply default style
    n = len(df)
    arrays = {k: _full(len(df), v) for k, v in style.items()}
    for (selector, override) in highlights_default:
        mask = _resolve(df, selector)
        for key, v in override.items():
            if key not in arrays:                 # property not in base → create it
                arrays[key] = _full(n, v)         
            a = arrays[key]
            if a.dtype == object:                 # assign as list so tuples don't broadcast
                a[mask] = [v] * int(np.count_nonzero(mask))
            else:
                a[mask] = v
                
    #arrays = {k: _full(len(df), v) for k, v in style.items()}
    # apply custom style if any
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

    # mpl is annoying those cn only be passed per call not per array
    # use this so modified point are plotted on top.
    
    sub = dict((k, v)  for (k, v) in arrays.items() if k in CALL_LEVEL)
    sub["_top_layer"] = drawn_last
    group_df = pd.DataFrame(sub)
    indices_same_call_level = group_df.groupby(list(sub), sort=False, dropna=False).indices

    
    for index in sorted(indices_same_call_level.values(), key=lambda p: drawn_last[p[0]]):
        per_call = {k: arrays[k][index[0]] for k in CALL_LEVEL if k in arrays }
        per_array = {k: v[index] for (k,v) in arrays.items() if k not in CALL_LEVEL }
        ax.scatter(df[x].values[index], df["neglog10"].values[index],
                   **per_array, **per_call)

    
    texts_ = []
    genes = []

    if label_top:
        genes.extend(list(df[(df["x"] < -x_thr) ].sort_values("y").head(label_top).index))
        genes.extend(list(df[(df["x"] > x_thr) ].sort_values("y").head(label_top).index))

    if label_genes: 
        if isinstance(label_genes, str):
            label_genes = [label_genes]

        missing = set(label_genes) - set(df.index)

        if missing:
            # TODO LOG
            print(f"WARNING: genes not found: {missing}", file=sys.stderr)

        genes.extend(label_genes)

    if genes:

        sub = df[df.index.isin(genes)]
        text = list(sub.index)
        x_axis = x 
        X = sub.apply(lambda x: x[x_axis], axis=1)
        Y = sub.apply(lambda x: x["neglog10"], axis=1)
        for i, v in enumerate(text):
            texts_.append(ax.text(x=X.iloc[i], y=Y.iloc[i], s=v, fontsize=10))


    if texts_:
        if label_genes_line:
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

    if show_thresholds_line:
        if y_thr is not None:
            h_line = -math.log(y_thr, 10)
            ax.axhline(h_line, color="black", linestyle="--")

        if x_thr is not None:
            v_line = x_thr
            ax.axvline(v_line, color="black", linestyle="--")
            ax.axvline(-v_line, color="black", linestyle="--")
            
    ax.tick_params(axis='y', labelsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.margins(x=0.02, tight=True)
    ax.margins(y=0.02, tight=True)  
    
    ax.set_title(title, fontsize=28)
    ax.spines[['right', 'top']].set_visible(False)

    # plot a quick summary of the data
    txt = ''
    if x_thr and y_thr:
        up = df[(df["neglog10"]< -math.log(y_thr, 10)) & (df[x] > x_thr)].shape[0]
        down =  df[(df["neglog10"]< -math.log(y_thr, 10)) & (df[x] < -x_thr)].shape[0]
        ns = df.shape[0] - (up + down)
        txt = "\n\nup {}; down {}; ns {}".format(up, down, ns )

    ax.set_ylabel(ylabel, fontsize=24)
    ax.set_xlabel(xlabel , fontsize=24)

    if txt:
        plt.figtext(0.5, -0.05, txt, wrap=True, horizontalalignment='center', fontsize=16)

    ##
    legend_elements = []
    for k, v in legend.items():
        c = v.get("c")
        if not c: 
            c = v.get("color")
        m = v.get("marker", ".")
        legend_elements.append(Line2D([0], [0],  marker=m,
                                       markerfacecolor=c, color='w', label=k, markersize=25))
    if legend_elements:
        ax.legend(handles=legend_elements, bbox_to_anchor=legend_bbox)

    # can't be used if custom_ax is set
    if savefig:
        fig.savefig(savefig, bbox_inches="tight")
    if return_fig:
        return fig, ax
    if custom_ax: 
        return ax


