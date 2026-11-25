"""Tests for easy_volcano.

Tests describe the *intended* behaviour (docstring + common sense for users),
so a failing test means either a bug or a design decision to revisit.
"""
import logging
import math

import matplotlib
matplotlib.use("Agg")                     # no display needed, works in CI
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import easy_volcano as ev

X, Y = "log2FoldChange", "padj"

COLORS = dict(color_sig="sig", color_ns="ns", color_up="up", color_down="down",
              color_up_sig="up_sig", color_down_sig="down_sig")


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #

@pytest.fixture
def df():
    """One gene per category, plus edge cases (padj == 0, NA)."""
    return pd.DataFrame(
        {X:           [3,    -3,   0.1, 2,   np.nan, 0.2,  4,   -4],
         Y:           [1e-5, 1e-4, 0.5, 0.0, 0.01,   1e-3, 0.2, 0.3],
         "gene_name": ["A",  "B",  "C", "D", "E",    "F",  "G", "H"]},
        index=[f"ENS{i}" for i in range(1, 9)])
    # ENS1 up_sig | ENS2 down_sig | ENS3 ns | ENS4 up_sig (p=0) | ENS5 NA (dropped)
    # ENS6 sig (p ok, small lfc) | ENS7 up (big lfc, p not ok) | ENS8 down


@pytest.fixture
def prepped(df):
    return ev.prep_df(df, X, Y, pval_floor=False)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def spy_ax(monkeypatch):
    """An Axes whose scatter() calls are recorded (kwargs only)."""
    fig, ax = plt.subplots()
    calls = []
    original = ax.scatter

    def scatter(*args, **kwargs):
        calls.append({"x": np.asarray(args[0]), **kwargs})
        return original(*args, **kwargs)

    monkeypatch.setattr(ax, "scatter", scatter)
    return ax, calls


def _default_colors(prepped):
    hd = ev.get_highlight_default(0.01, 1, False, X, Y, **COLORS)
    style = {"c": "base", "zorder": 1, "s": 5, "alpha": 1, "marker": "."}
    arrays, _ = ev.get_array_style(prepped, hd, [], style, None)
    return dict(zip(prepped.index, arrays["c"]))


# --------------------------------------------------------------------------- #
# prep_df
# --------------------------------------------------------------------------- #

def test_prep_df_drops_na_and_does_not_mutate(df):
    before = df.copy()
    out = ev.prep_df(df, X, Y, pval_floor=False)
    assert "ENS5" not in out.index
    assert len(out) == 7
    pd.testing.assert_frame_equal(df, before)


def test_prep_df_zero_pvalue_floored_to_min_observed(prepped):
    # smallest positive padj is 1e-5 -> padj=0 becomes 5
    assert prepped.loc["ENS4", "neglog10"] == pytest.approx(5)
    assert np.isfinite(prepped["neglog10"]).all()


def test_prep_df_float_floor(df):
    out = ev.prep_df(df, X, Y, pval_floor=True)
    assert out.loc["ENS4", "neglog10"] == pytest.approx(-math.log10(np.finfo(float).tiny))


@pytest.mark.parametrize("pval_floor", [False, True])
def test_prep_df_no_positive_pvalue_raises(df, pval_floor):
    # no positive p-value is never meaningful biology: refuse to plot
    df[Y] = 0.0
    with pytest.raises(ValueError, match="padj"):
        ev.prep_df(df, X, Y, pval_floor=pval_floor)


def test_prep_df_logs_dropped_rows(df, caplog):
    with caplog.at_level(logging.INFO, logger="easy_volcano"):
        ev.prep_df(df, X, Y, pval_floor=False)
    assert "Dropped 1" in caplog.text


def test_plot_all_na_raises_clear_error(df):
    df[X] = np.nan
    with pytest.raises(ValueError, match="empty after"):
        ev.plot_volcano(df)


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #

def test_full_numeric_is_float():
    a = ev._full(3, 5)
    assert a.dtype == float and a.tolist() == [5.0] * 3


def test_full_keeps_tuples_whole():
    a = ev._full(3, (1, 0, 0))
    assert a.dtype == object and all(v == (1, 0, 0) for v in a)


def test_full_bool_is_not_numeric():
    assert ev._full(2, True).dtype == object


def test_norm_symbol_aliases():
    assert ev.mpl_norm_symbol({"size": 3, "color": "r", "a": 0.5, "marker": "x"}) == \
        {"s": 3, "c": "r", "alpha": 0.5, "marker": "x"}


def test_norm_symbol_conflict_raises():
    with pytest.raises(ValueError, match="only use one"):
        ev.mpl_norm_symbol({"size": 3, "s": 4})


def test_norm_symbol_does_not_mutate():
    d = {"size": 3}
    ev.mpl_norm_symbol(d)
    assert d == {"size": 3}


# --------------------------------------------------------------------------- #
# _resolve
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("sel, label_col", [
    (lambda d: d["gene_name"] == "A", None),
    (["ENS1"], None),
    (("ENS1",), None),
    ({"ENS1"}, None),
    (pd.Index(["ENS1"]), None),
    (["A"], "gene_name"),
    (pd.Series([True], index=["ENS1"]), None),            # partial index -> rest False
])
def test_resolve_selector_types(prepped, sel, label_col):
    mask = ev._resolve(prepped, sel, label_col=label_col)
    assert mask.dtype == bool
    assert prepped.index[mask].tolist() == ["ENS1"]


def test_resolve_series_built_on_raw_df_aligns_after_dropna(df, prepped):
    # users build masks on their own df, which still has the NA row
    mask = ev._resolve(prepped, df[Y] < 1e-3)
    assert prepped.index[mask].tolist() == ["ENS1", "ENS2", "ENS4"]


def test_resolve_bare_string_is_rejected(prepped):
    # a str is iterable; make sure "ENS1" isn't silently treated as chars
    with pytest.raises(TypeError):
        ev._resolve(prepped, "ENS1")


def test_resolve_unknown_names_match_nothing(prepped):
    assert not ev._resolve(prepped, ["nope"]).any()


# --------------------------------------------------------------------------- #
# default colouring
# --------------------------------------------------------------------------- #

def test_default_categories_follow_docstring(prepped):
    assert _default_colors(prepped) == {
        "ENS1": "up_sig", "ENS2": "down_sig", "ENS4": "up_sig",
        "ENS3": "ns",     # passes nothing
        "ENS6": "sig",    # passes p-value only
        "ENS7": "up",     # passes lfc only
        "ENS8": "down",
    }


def test_threshold_boundaries_are_consistent():
    # exactly on a threshold counts as passing it (p <= y_thr, |lfc| >= x_thr)
    d = pd.DataFrame({X: [1, -1, 0.5, 1, 0.5],
                      Y: [0.5, 0.5, 0.01, 0.01, 0.5]},
                     index=["up", "down", "sig", "up_sig", "ns"])
    colors = _default_colors(ev.prep_df(d, X, Y, pval_floor=False))
    assert colors == {k: k for k in d.index}


def test_every_point_gets_a_category_colour(prepped):
    assert "base" not in _default_colors(prepped).values()


@pytest.mark.parametrize("kw", [dict(no_default=True), dict(x_thr=None), dict(y_thr=None)])
def test_no_default_colouring(kw):
    args = dict(y_thr=0.01, x_thr=1, no_default=False) | kw
    assert ev.get_highlight_default(x=X, y=Y, **args, **COLORS) == []


# --------------------------------------------------------------------------- #
# highlights / style arrays
# --------------------------------------------------------------------------- #

def _style(prepped, highlights, label_column=None):
    style = {"c": "base", "zorder": 1, "s": 5, "alpha": 1, "marker": "."}
    highlights = [(s, ev.mpl_norm_symbol(o)) for s, o in highlights]
    return ev.get_array_style(prepped, [], highlights, style, label_column)


def test_highlight_overrides_only_selected(prepped):
    arrays, drawn_last = _style(prepped, [(["ENS1"], {"color": "red", "size": 30})])
    assert arrays["c"].tolist() == ["red"] + ["base"] * 6
    assert arrays["s"].tolist() == [30.0] + [5.0] * 6
    assert drawn_last.tolist() == [True] + [False] * 6


def test_highlight_by_label_column(prepped):
    arrays, _ = _style(prepped, [(["B"], {"color": "red"})], label_column="gene_name")
    assert arrays["c"][prepped.index.get_loc("ENS2")] == "red"


def test_later_highlight_wins(prepped):
    arrays, _ = _style(prepped, [(["ENS1", "ENS2"], {"color": "red"}),
                                 (["ENS2"], {"color": "blue"})])
    assert arrays["c"][:2].tolist() == ["red", "blue"]


def test_new_property_only_applies_to_selected_points(prepped):
    # 'edgecolors' isn't in the base style: unselected points must not get it
    arrays, _ = _style(prepped, [(["ENS1"], {"edgecolors": "black"})])
    assert arrays["edgecolors"][0] == "black"
    assert all(v != "black" for v in arrays["edgecolors"][1:])


def test_rgba_tuple_colour(prepped):
    arrays, _ = _style(prepped, [(["ENS1"], {"color": (1, 0, 0, 1)})])
    assert arrays["c"][0] == (1, 0, 0, 1)


# --------------------------------------------------------------------------- #
# which points get a text label
# --------------------------------------------------------------------------- #

def _labels(prepped, highlight=(), label_format=None, label_genes_format=None, **kw):
    """Run get_index_point_to_labels and flatten to {row index: (text, format)}."""
    groups = ev.get_index_point_to_labels(
        prepped,
        highlight=[(sel, dict(ov)) for sel, ov in highlight],   # function pops keys
        label_format=ev._parse_label_format(label_format),
        label_genes_format=label_genes_format,
        x=X, y=Y, **kw)
    return {i: (t, fmt)
            for (xs, _, text, fmt) in groups
            for i, t in zip(xs.index, text)}


def test_label_top_with_label_column(prepped):                  # regression
    rows = _labels(prepped, x_thr=1, label_top=1, label_column="gene_name")
    assert set(rows) == {"ENS2", "ENS4"}


def test_label_genes_by_index(prepped):
    rows = _labels(prepped, label_genes=["ENS1"])
    assert list(rows) == ["ENS1"]


def test_label_genes_single_string(prepped):
    rows = _labels(prepped, label_genes="A", label_column="gene_name")
    assert list(rows) == ["ENS1"]


@pytest.mark.parametrize("container", [list, set, tuple, pd.Index, pd.Series, np.array])
def test_label_genes_accepts_common_containers(prepped, container):
    rows = _labels(prepped, label_genes=container(["A", "B"]), label_column="gene_name")
    assert sorted(rows) == ["ENS1", "ENS2"]


def test_label_genes_missing_warns(prepped, caplog):
    with caplog.at_level(logging.WARNING, logger="easy_volcano"):
        rows = _labels(prepped, label_genes=["A", "ZZZ"], label_column="gene_name")
    assert list(rows) == ["ENS1"]
    assert "ZZZ" in caplog.text


def test_label_top_without_x_thr_does_not_crash(prepped):
    rows = _labels(prepped, x_thr=None, label_top=1)
    assert len(rows) == 2


def test_label_text_uses_label_column(prepped):
    rows = _labels(prepped, label_genes=["A"], label_column="gene_name")
    assert rows["ENS1"][0] == "A"


# --------------------------------------------------------------------------- #
# per-highlight labels (label_gene) and label_genes_format
# --------------------------------------------------------------------------- #

def test_highlight_label_gene_true_labels_selection(prepped):
    rows = _labels(prepped, highlight=[(["ENS1", "ENS2"], {"color": "red", "label_gene": True})])
    assert sorted(rows) == ["ENS1", "ENS2"]


def test_highlight_without_label_gene_adds_no_label(prepped):
    rows = _labels(prepped, highlight=[(["ENS1"], {"color": "red"})])
    assert rows == {}


@pytest.mark.parametrize("falsy", [False, None, {}])
def test_highlight_label_gene_falsy_adds_no_label(prepped, falsy):
    rows = _labels(prepped, highlight=[(["ENS1"], {"label_gene": falsy})])
    assert rows == {}


def test_highlight_label_gene_callable_selector(prepped):
    rows = _labels(prepped, highlight=[(lambda d: d[X] < -2, {"label_gene": True})])
    assert sorted(rows) == ["ENS2", "ENS8"]


def test_highlight_label_gene_true_uses_global_format(prepped):
    rows = _labels(prepped, label_format={"fontsize": 14},
                   highlight=[(["ENS1"], {"label_gene": True})])
    assert rows["ENS1"][1]["fontsize"] == 14


def test_highlight_label_gene_dict_merges_on_global_format(prepped):
    rows = _labels(prepped, label_format={"fontsize": 14},
                   highlight=[(["ENS1"], {"label_gene": {"fontweight": "bold"}})])
    fmt = rows["ENS1"][1]
    assert (fmt["fontsize"], fmt["fontweight"]) == (14, "bold")


def test_highlight_label_gene_reserved_key(prepped):
    with pytest.raises(ValueError, match="'x'"):
        _labels(prepped, highlight=[(["ENS1"], {"label_gene": {"x": 1}})])


def test_label_genes_format_applies_to_label_genes_only(prepped):
    rows = _labels(prepped, label_format={"fontsize": 14},
                   label_genes=["ENS1"], label_genes_format={"color": "red"},
                   label_top=1)
    assert rows["ENS1"][1]["color"] == "red"
    assert rows["ENS1"][1]["fontsize"] == 14          # still built on the global format
    assert "color" not in rows["ENS2"][1]             # ENS2 comes from label_top


def test_label_genes_format_reserved_key(prepped):
    with pytest.raises(ValueError, match="'text'"):
        _labels(prepped, label_genes=["ENS1"], label_genes_format={"text": "z"})


def test_label_priority_highlight_over_label_genes_over_top(prepped):
    # ENS1 is picked by all three sources, ENS4 by label_genes and label_top
    rows = _labels(prepped, label_top=1,
                   label_genes=["ENS1", "ENS4"], label_genes_format={"color": "green"},
                   highlight=[(["ENS1"], {"label_gene": {"color": "red"}})])
    assert rows["ENS1"][1]["color"] == "red"
    assert rows["ENS4"][1]["color"] == "green"


def test_label_each_gene_once_across_sources(prepped):
    groups = ev.get_index_point_to_labels(
        prepped, highlight=[(["ENS1"], {"label_gene": True})],
        label_format=ev._parse_label_format(None), label_genes_format=None,
        x=X, y=Y, label_genes=["ENS1", "ENS2"], label_top=1)
    labelled = [i for (xs, *_ ) in groups for i in xs.index]
    assert sorted(labelled) == sorted(set(labelled))


def test_later_highlight_label_format_wins(prepped):
    rows = _labels(prepped, highlight=[(["ENS1"], {"label_gene": {"color": "red"}}),
                                       (["ENS1"], {"label_gene": {"color": "blue"}})])
    assert rows["ENS1"][1]["color"] == "blue"


def test_label_position_matches_point(prepped):
    groups = ev.get_index_point_to_labels(
        prepped, highlight=[], label_format=ev._parse_label_format(None),
        label_genes_format=None, x=X, y=Y, label_genes=["ENS1"])
    xs, ys, _, _ = groups[0]
    assert xs.iloc[0] == prepped.loc["ENS1", X]
    assert ys.iloc[0] == prepped.loc["ENS1", "neglog10"]


# --------------------------------------------------------------------------- #
# plot_volcano: API contract
# --------------------------------------------------------------------------- #

def test_plot_returns_fig_and_ax(df):
    fig, ax = ev.plot_volcano(df)
    assert ax.figure is fig


def test_plot_given_ax_is_used_and_returned(df):
    user_fig, user_ax = plt.subplots()
    fig, ax = ev.plot_volcano(df, ax=user_ax)
    assert ax is user_ax
    assert fig is user_fig
    assert len(user_ax.collections) > 0


def test_options_are_keyword_only(df):
    with pytest.raises(TypeError):
        ev.plot_volcano(df, False)


@pytest.mark.parametrize("ext", ["png", "pdf", "svg"])
def test_savefig(df, tmp_path, ext):
    path = tmp_path / f"volcano.{ext}"
    ev.plot_volcano(df, savefig=path)
    assert path.stat().st_size > 0


def test_savefig_with_given_ax(df, tmp_path):
    path = tmp_path / "volcano.png"
    _, ax = plt.subplots()
    ev.plot_volcano(df, ax=ax, savefig=path)
    assert path.stat().st_size > 0


def test_plot_does_not_mutate_input(df):
    before = df.copy()
    ev.plot_volcano(df, label_top=1, highlights=[(["ENS1"], {"color": "red", "label": "x"})])
    pd.testing.assert_frame_equal(df, before)


def test_plot_does_not_mutate_highlights(df):
    hl = [(["ENS1"], {"color": "red", "label": "mine", "zorder": 3})]
    ev.plot_volcano(df, highlights=hl)
    assert hl == [(["ENS1"], {"color": "red", "label": "mine", "zorder": 3})]


def test_custom_column_names(df):
    df = df.rename(columns={X: "logFC", Y: "adj.P.Val"})
    fig, ax = ev.plot_volcano(df, x="logFC", y="adj.P.Val")
    assert len(ax.collections) > 0


def test_duplicate_index(df):
    df.index = ["g1", "g1", "g2", "g3", "g4", "g5", "g6", "g6"]
    _, ax = ev.plot_volcano(df, label_genes=["g1"], label_top=None)
    assert [t.get_text() for t in ax.texts] == ["g1", "g1"]


# --------------------------------------------------------------------------- #
# plot_volcano: what ends up on the axes
# --------------------------------------------------------------------------- #

def test_every_point_plotted_once(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, ax=ax, highlights=[(["ENS1"], {"color": "red"})])
    xs = np.concatenate([c["x"] for c in calls])
    assert sorted(xs) == sorted(df[X].dropna())


def test_highlights_drawn_last(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, ax=ax, highlights=[(["ENS1"], {"color": "red"})])
    assert calls[-1]["x"].tolist() == [3]


def test_highlight_zorder_is_applied(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, ax=ax, highlights=[(["ENS1"], {"color": "red", "zorder": 5})])
    assert any(c.get("zorder") == 5 for c in calls)


def test_new_property_not_passed_to_other_points(df, spy_ax):
    ax, calls = spy_ax
    # ENS2 is also on the top layer: only the "unset" grouping keeps it edge-free
    ev.plot_volcano(df, ax=ax, highlights=[(["ENS1"], {"edgecolors": "black"}),
                                                  (["ENS2"], {"color": "red"})])
    with_edges = [c for c in calls if "edgecolors" in c]
    assert [c["x"].tolist() for c in with_edges] == [[3]]
    assert list(with_edges[0]["edgecolors"]) == ["black"]
    # every other point is still drawn, just without the argument
    assert sum(len(c["x"]) for c in calls) == 7


def test_marker_kwarg_applies_to_all_points(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, ax=ax, marker="^")
    assert {c["marker"] for c in calls} == {"^"}


def test_labels_use_label_column(df):                           # regression
    _, ax = ev.plot_volcano(df, label_genes=["A", "B"], label_column="gene_name",
                            label_top=None)
    assert sorted(t.get_text() for t in ax.texts) == ["A", "B"]


@pytest.mark.parametrize("container", [pd.Index, pd.Series, np.array])
def test_plot_label_genes_containers(df, container):
    _, ax = ev.plot_volcano(df, label_genes=container(["A", "B"]),
                            label_column="gene_name", label_top=None)
    assert sorted(t.get_text() for t in ax.texts) == ["A", "B"]


def test_labels_default_to_index(df):
    _, ax = ev.plot_volcano(df, label_genes=["ENS1"], label_top=None)
    assert [t.get_text() for t in ax.texts] == ["ENS1"]


def test_label_top_and_label_genes_no_duplicates(df):
    _, ax = ev.plot_volcano(df, label_top=1, label_genes=["B", "F"],
                            label_column="gene_name")
    assert sorted(t.get_text() for t in ax.texts) == ["B", "D", "F"]


def test_no_labels_requested(df):
    _, ax = ev.plot_volcano(df, label_top=None)
    assert len(ax.texts) == 0


def test_default_legend(df):
    _, ax = ev.plot_volcano(df)
    assert [t.get_text() for t in ax.get_legend().get_texts()] == \
        ["Overexpressed", "Underexpressed"]


def test_highlight_label_added_to_legend(df):
    _, ax = ev.plot_volcano(df,
                            highlights=[(["ENS1"], {"color": "red", "label": "mine"})])
    assert "mine" in [t.get_text() for t in ax.get_legend().get_texts()]


def test_no_default_has_no_legend(df):
    _, ax = ev.plot_volcano(df, no_default=True)
    assert ax.get_legend() is None


def test_show_legend_off(df):
    _, ax = ev.plot_volcano(df, show_legend=False,
                            highlights=[(["ENS1"], {"color": "red", "label": "mine"})])
    assert ax.get_legend() is None


@pytest.mark.parametrize("rasterize", [False, True])
def test_rasterize(df, rasterize):
    _, ax = ev.plot_volcano(df, rasterize=rasterize,
                            highlights=[(["ENS1"], {"color": "red"})])
    assert [c.get_rasterized() for c in ax.collections] == \
        [rasterize] * len(ax.collections)


def test_rasterize_overrides_highlight(df):
    _, ax = ev.plot_volcano(df, rasterize=True,
                            highlights=[(["ENS1"], {"color": "red", "rasterized": False})])
    assert all(c.get_rasterized() for c in ax.collections)


def test_rasterized_per_highlight(df):
    _, ax = ev.plot_volcano(df, highlights=[(["ENS1"], {"color": "red", "rasterized": True})])
    flags = [c.get_rasterized() for c in ax.collections]
    assert True in flags and False in flags


def test_threshold_lines(df):
    _, ax = ev.plot_volcano(df, x_thr=1, y_thr=0.01)
    h = [l.get_ydata()[0] for l in ax.lines if l.get_ydata()[0] == l.get_ydata()[1]]
    v = [l.get_xdata()[0] for l in ax.lines if l.get_xdata()[0] == l.get_xdata()[1]]
    assert h == [pytest.approx(2)]
    assert sorted(v) == [-1, 1]


def test_threshold_lines_off(df):
    _, ax = ev.plot_volcano(df, show_thresholds_line=False)
    assert len(ax.lines) == 0


def test_title_and_axis_labels(df):
    _, ax = ev.plot_volcano(df, title="T", xlabel="LFC", ylabel="P")
    assert (ax.get_title(), ax.get_xlabel(), ax.get_ylabel()) == ("T", "LFC", "P")


# --------------------------------------------------------------------------- #
# summary text under the plot
# --------------------------------------------------------------------------- #

def _summary(fig):
    return " ".join(t.get_text() for t in fig.texts)


def test_summary_counts(df):
    # after dropna: up_sig = ENS1, ENS4; down_sig = ENS2; 4 others
    fig, _ = ev.plot_volcano(df)
    assert "up 2; down 1; ns 4" in _summary(fig)


def test_summary_with_zero_lfc_threshold(df):
    fig, _ = ev.plot_volcano(df, x_thr=0)
    assert "up" in _summary(fig)


def test_summary_goes_on_given_ax_figure(df):
    fig1, ax1 = plt.subplots()
    fig2 = plt.figure()                        # now the "current" figure
    ev.plot_volcano(df, ax=ax1)
    assert "up" in _summary(fig1)
    assert _summary(fig2) == ""


def test_no_summary_without_thresholds(df):
    fig, _ = ev.plot_volcano(df, x_thr=None, y_thr=None)
    assert _summary(fig) == ""


def test_show_text_summary_off(df):
    fig, _ = ev.plot_volcano(df, show_text_summary=False)
    assert _summary(fig) == ""


def test_axis_labels_without_summary(df):
    # turning off the summary should not affect the axis labels
    _, ax = ev.plot_volcano(df, show_text_summary=False, xlabel="LFC", ylabel="P")
    assert (ax.get_xlabel(), ax.get_ylabel()) == ("LFC", "P")


# --------------------------------------------------------------------------- #
# label_format
# --------------------------------------------------------------------------- #

def test_parse_label_format_defaults():
    assert ev._parse_label_format(None) == ev._LABEL_FORMAT_DEFAULT


def test_parse_label_format_partial_keeps_defaults():
    lf = ev._parse_label_format({"color": "red"})
    assert lf == {**ev._LABEL_FORMAT_DEFAULT, "color": "red"}


def test_parse_label_format_does_not_mutate_default():
    before = dict(ev._LABEL_FORMAT_DEFAULT)
    ev._parse_label_format({"fontsize": 20})
    assert ev._LABEL_FORMAT_DEFAULT == before


@pytest.mark.parametrize("key", ["x", "y", "s", "position", "text"])
def test_parse_label_format_reserved_keys(key):
    with pytest.raises(ValueError, match=key):
        ev._parse_label_format({key: None})


def test_label_format_applied_to_labels(df):
    _, ax = ev.plot_volcano(df, label_genes=["ENS1"],
                            label_format={"fontsize": 14, "color": "red"})
    t = ax.texts[0]
    assert (t.get_fontsize(), t.get_color()) == (14, "red")


def test_label_format_validated_without_labels(df):
    with pytest.raises(ValueError):
        ev.plot_volcano(df, label_format={"x": 1})


def test_plot_highlight_label_gene_draws_formatted_labels(df):
    _, ax = ev.plot_volcano(df, label_format={"fontsize": 14}, label_top=None,
                            highlights=[(["ENS1"], {"color": "red",
                                                    "label_gene": {"fontweight": "bold"}})])
    [t] = ax.texts
    assert (t.get_text(), t.get_fontsize(), t.get_fontweight()) == ("ENS1", 14, "bold")


def test_plot_label_gene_not_passed_to_scatter(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, ax=ax, highlights=[(["ENS1"], {"color": "red", "label_gene": True})])
    assert all("label_gene" not in c for c in calls)


def test_plot_label_gene_does_not_mutate_highlights(df):
    highlights = [(["ENS1"], {"color": "red", "label_gene": {"color": "blue"}})]
    ev.plot_volcano(df, highlights=highlights)
    assert highlights == [(["ENS1"], {"color": "red", "label_gene": {"color": "blue"}})]


def test_plot_label_genes_format(df):
    _, ax = ev.plot_volcano(df, label_genes=["ENS1"], label_genes_format={"color": "red"},
                            label_top=1)
    colors = {t.get_text(): t.get_color() for t in ax.texts}
    assert colors["ENS1"] == "red"
    assert colors["ENS2"] != "red"                  # from label_top, global format


def test_plot_label_gene_with_label_column(df):
    _, ax = ev.plot_volcano(df, label_column="gene_name", label_top=None,
                            highlights=[(["A", "B"], {"label_gene": True})])
    assert sorted(t.get_text() for t in ax.texts) == ["A", "B"]
