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
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib as mpl
import numbers
import numpy as np
import pandas as pd
import sys
import adjustText
import functools

import logging
logger = logging.getLogger(__name__)

mpl.rcParams["font.stretch"] = "expanded"
plt.rcParams['pdf.fonttype'] = 42

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
    
    lines = lines if isinstance(lines, dict) else [lines]
    return lines 
def _parse_label_format(label_format, _BLACK_LIST=_BLACK_LIST_LABEL_FORMAT):

    label_format = {**_LABEL_FORMAT_DEFAULT, **(label_format or {})}
    for e in _BLACK_LIST:

        if e in  label_format:
            raise ValueError(f"label_format cannot contain {e!r} as key")

    return label_format

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
                 label_top = None, # take int label the top X gene over and under expressed (so *2 labels)
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
                size=50, color_base="grey",
                width=8, height=4.8,
                ax = None,
                rasterize = False, # if set override any other arg and rasterize all
                savefig=None,

                title="",
                title_pad=30,
                title_fontsize=28,

                show_legend= True,

                legend_bbox=(1.0, 0.50),
                legend_loc = "right center",
                legend_ncols = 2,
                legend_fancybox=True,

                show_text_summary = True,
                **scatter_kwargs):


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

             highlights special keyword: 
                - boxes_label : boolean (in future may support dict)
                - label : str -> add str to legend
                - label_gene : None or dict indicate that those genes must have their label set, the dict can take any text formatting option from mpl. 

    Parameters
    ----------
    df : pd.DataFrame
        Rows containing any NA are dropped before plotting.
    color_bg : str, default "white",
    x : str , default "log2FoldChange"
        column used for the X-axis
    y : str , default "padj",
        column used for the Y-axis
    x_thr : float, default 1,
        X-axis threshold used for coloring, summarizing count and plotting lines
    y_thr : float, default 0.01, 
        Y-axis threshold used for coloring, summarizing count and plotting lines
    
    ----------
    ## p-val cap
        By default easyvolcano shrink the Y axis to the smallest value present in df.
        p-value of zeros are transformed to this values (needed for the log transform). 
    clip_pval : float: default None
        if set superseed pval_floor, clip all pvalue smaller than 'clip_pval' to 'clip_pval'
    pval_floor : boolean, default False,
        if set  p-value of zeros are transformed to sys.float_info.min. 
    ----------

    ----------
    ## label
    label_genes = None,
    label_genes_line = False,
    label_top = None, # take int label the top X gene over and under expressed (so *2 labels)
    label_column = None, # inseatd of using the index look into this column for label name # label_genes must be from this columns
    label_format = _LABEL_FORMAT_DEFAULT,
    ----------

    alpha : float, default 1
        Point transparency passed to scatter.
    width, float, default 8
        Figure width in inches.
    height : float, default, 4.5
        Figure height in inches.

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
        [(lambda d: (d[y] >= 2) & (abs(d[x]) <= 1), { "color": "blue"})]

        USing a label key work in the override, will add it to the legends
        list selectors match label_column if set, otherwise the index.
s
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
    (fig ax)
        The figure is drawn
        optional:
            if savefig is not None and is a vlaid file save the figure
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
    legend = prep_legend(highlights_default, highlights, **scatter_kwargs)
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
            adjustText.adjust_text(texts_, expand=(1.1, 1.2),
                                   min_arrow_len=label_min_arrow_len,
                                   avoid_self=label_avoid_self,
                                   arrowprops=dict(arrowstyle='-', color='black'))
        else:
            adjustText.adjust_text(texts_, expand=(1.1, 1.2),
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
            y = copy_fmt.pop("y", None)
            # neglog =  copy_line.pop("neglog", False) not sure about that option; to niche to be usefull
            if (y is None) or not (isinstance(y, numbers.Real) or (isinstance(y, bool))):
                raise ValueError("h_lines dict(s) must contain a 'y' key with a numerical value (pvalue, plot position)")
            if y <= 0:
                raise ValueError("'y' cannot be equal to zeros or negative, undefined log10")
            
            copy_fmt = {**_LINES_FORMAT_DEFAULT, **copy_fmt}
            ax.axhline(-math.log10(y), **copy_fmt)

    if v_lines:

        v_lines = _check_lines_param(v_lines)
        for v_line in v_lines:

            if not isinstance(v_line, dict):
                raise ValueError("v_lines parameter must be a list of dict or a dict")
            
            copy_fmt = dict(v_line)
            x = copy_fmt.pop("x", None)

            if (x is None) or not (isinstance(x, numbers.Real) or (isinstance(x, bool))):
                raise ValueError("v_lines dict(s) must contain a 'x' key with a numerical value (pvalue, plot position)")
            
            copy_fmt = {**_LINES_FORMAT_DEFAULT, **copy_fmt}
            ax.axvline(x, **copy_fmt)

            
    ax.tick_params(axis='y', labelsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.margins(x=0.02, tight=True)
    ax.margins(y=0.02, tight=True)  
    if title:
        ax.set_title(title, fontsize=title_fontsize, pad=title_pad)
    ax.spines[['right', 'top']].set_visible(False)

    if show_text_summary:
    # plot a quick summary of the data
        txt = ''
        if x_thr is not None and y_thr is not None:
            up = df[(df[y] <= y_thr) & (df[x] >= x_thr)].shape[0]
            down =  df[(df[y] <= y_thr) & (df[x] <= -x_thr)].shape[0]
            ns = df.shape[0] - (up + down)
            txt = "\n\nup {}; down {}; ns {}".format(up, down, ns )
        if txt:
            ax.figure.text(0.5, -0.05, txt, wrap=True, horizontalalignment='center', fontsize=16)
    
    ax.set_ylabel(ylabel, fontsize=24)
    ax.set_xlabel(xlabel , fontsize=24)

        ##
    bb = ax.get_tightbbox(fig.canvas.get_renderer())    # axes + ticks + xlabel
    y_bottom = ax.transAxes.inverted().transform((0, bb.y0))[1]
    if show_legend:
        legend_elements = []
        for k, v in legend.items():
            c = v.get("c")
            if not c: 
                c = v.get("color")
            m = v.get("marker", ".")
            legend_elements.append(Line2D([0], [0],  marker=m,
                                        markerfacecolor=c, color='w', label=k, markersize=18))
        if legend_elements:
            ax.legend(handles=legend_elements, bbox_to_anchor=legend_bbox, fancybox=legend_fancybox,
                       loc=legend_loc, ncols=legend_ncols)

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