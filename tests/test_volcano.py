"""Tests for easyvolcano.

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

import easyvolcano as ev

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
    with caplog.at_level(logging.INFO, logger="easyvolcano"):
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

def test_label_top_with_label_column(prepped):                  # regression
    rows = ev.get_index_point_to_labels(prepped, x_thr=1, label_top=1,
                                        label_genes=None, label_column="gene_name")
    assert set(rows) == {"ENS2", "ENS4"}


def test_label_genes_by_index(prepped):
    rows = ev.get_index_point_to_labels(prepped, label_genes=["ENS1"])
    assert list(rows) == ["ENS1"]


def test_label_genes_single_string(prepped):
    rows = ev.get_index_point_to_labels(prepped, label_genes="A", label_column="gene_name")
    assert list(rows) == ["ENS1"]


@pytest.mark.parametrize("container", [list, set, tuple, pd.Index, pd.Series, np.array])
def test_label_genes_accepts_common_containers(prepped, container):
    rows = ev.get_index_point_to_labels(prepped, label_genes=container(["A", "B"]),
                                        label_column="gene_name")
    assert sorted(rows) == ["ENS1", "ENS2"]


def test_label_genes_missing_warns(prepped, caplog):
    with caplog.at_level(logging.WARNING, logger="easyvolcano"):
        rows = ev.get_index_point_to_labels(prepped, label_genes=["A", "ZZZ"],
                                            label_column="gene_name")
    assert list(rows) == ["ENS1"]
    assert "ZZZ" in caplog.text


def test_label_top_without_x_thr_does_not_crash(prepped):
    rows = ev.get_index_point_to_labels(prepped, x_thr=None, label_top=1)
    assert len(rows) == 2


# --------------------------------------------------------------------------- #
# plot_volcano: API contract
# --------------------------------------------------------------------------- #

def test_plot_returns_none_by_default(df):
    assert ev.plot_volcano(df) is None


def test_plot_return_fig(df):
    fig, ax = ev.plot_volcano(df, return_fig=True)
    assert ax.figure is fig


def test_plot_custom_ax_returns_ax(df):
    _, ax = plt.subplots()
    assert ev.plot_volcano(df, custom_ax=ax) is ax


@pytest.mark.parametrize("kw", [dict(return_fig=True), dict(savefig="out.png")])
def test_custom_ax_incompatible_args(df, kw):
    _, ax = plt.subplots()
    with pytest.raises(ValueError):
        ev.plot_volcano(df, custom_ax=ax, **kw)


@pytest.mark.parametrize("ext", ["png", "pdf", "svg"])
def test_savefig(df, tmp_path, ext):
    path = tmp_path / f"volcano.{ext}"
    ev.plot_volcano(df, savefig=path)
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
    fig, ax = ev.plot_volcano(df, x="logFC", y="adj.P.Val", return_fig=True)
    assert len(ax.collections) > 0


def test_duplicate_index(df):
    df.index = ["g1", "g1", "g2", "g3", "g4", "g5", "g6", "g6"]
    _, ax = ev.plot_volcano(df, label_genes=["g1"], return_fig=True)
    assert [t.get_text() for t in ax.texts] == ["g1", "g1"]


# --------------------------------------------------------------------------- #
# plot_volcano: what ends up on the axes
# --------------------------------------------------------------------------- #

def test_every_point_plotted_once(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, custom_ax=ax, highlights=[(["ENS1"], {"color": "red"})])
    xs = np.concatenate([c["x"] for c in calls])
    assert sorted(xs) == sorted(df[X].dropna())


def test_highlights_drawn_last(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, custom_ax=ax, highlights=[(["ENS1"], {"color": "red"})])
    assert calls[-1]["x"].tolist() == [3]


def test_highlight_zorder_is_applied(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, custom_ax=ax, highlights=[(["ENS1"], {"color": "red", "zorder": 5})])
    assert any(c.get("zorder") == 5 for c in calls)


def test_new_property_not_passed_to_other_points(df, spy_ax):
    ax, calls = spy_ax
    # ENS2 is also on the top layer: only the "unset" grouping keeps it edge-free
    ev.plot_volcano(df, custom_ax=ax, highlights=[(["ENS1"], {"edgecolors": "black"}),
                                                  (["ENS2"], {"color": "red"})])
    with_edges = [c for c in calls if "edgecolors" in c]
    assert [c["x"].tolist() for c in with_edges] == [[3]]
    assert list(with_edges[0]["edgecolors"]) == ["black"]
    # every other point is still drawn, just without the argument
    assert sum(len(c["x"]) for c in calls) == 7


def test_marker_kwarg_applies_to_all_points(df, spy_ax):
    ax, calls = spy_ax
    ev.plot_volcano(df, custom_ax=ax, marker="^")
    assert {c["marker"] for c in calls} == {"^"}


def test_labels_use_label_column(df):                           # regression
    _, ax = ev.plot_volcano(df, label_genes=["A", "B"], label_column="gene_name",
                            return_fig=True)
    assert sorted(t.get_text() for t in ax.texts) == ["A", "B"]


@pytest.mark.parametrize("container", [pd.Index, pd.Series, np.array])
def test_plot_label_genes_containers(df, container):
    _, ax = ev.plot_volcano(df, label_genes=container(["A", "B"]),
                            label_column="gene_name", return_fig=True)
    assert sorted(t.get_text() for t in ax.texts) == ["A", "B"]


def test_labels_default_to_index(df):
    _, ax = ev.plot_volcano(df, label_genes=["ENS1"], return_fig=True)
    assert [t.get_text() for t in ax.texts] == ["ENS1"]


def test_label_top_and_label_genes_no_duplicates(df):
    _, ax = ev.plot_volcano(df, label_top=1, label_genes=["B", "F"],
                            label_column="gene_name", return_fig=True)
    assert sorted(t.get_text() for t in ax.texts) == ["B", "D", "F"]


def test_no_labels_requested(df):
    _, ax = ev.plot_volcano(df, return_fig=True)
    assert len(ax.texts) == 0


def test_default_legend(df):
    _, ax = ev.plot_volcano(df, return_fig=True)
    assert [t.get_text() for t in ax.get_legend().get_texts()] == \
        ["Overexpressed", "Underexpressed"]


def test_highlight_label_added_to_legend(df):
    _, ax = ev.plot_volcano(df, return_fig=True,
                            highlights=[(["ENS1"], {"color": "red", "label": "mine"})])
    assert "mine" in [t.get_text() for t in ax.get_legend().get_texts()]


def test_no_default_has_no_legend(df):
    _, ax = ev.plot_volcano(df, no_default=True, return_fig=True)
    assert ax.get_legend() is None


def test_threshold_lines(df):
    _, ax = ev.plot_volcano(df, x_thr=1, y_thr=0.01, return_fig=True)
    h = [l.get_ydata()[0] for l in ax.lines if l.get_ydata()[0] == l.get_ydata()[1]]
    v = [l.get_xdata()[0] for l in ax.lines if l.get_xdata()[0] == l.get_xdata()[1]]
    assert h == [pytest.approx(2)]
    assert sorted(v) == [-1, 1]


def test_threshold_lines_off(df):
    _, ax = ev.plot_volcano(df, show_thresholds_line=False, return_fig=True)
    assert len(ax.lines) == 0


def test_title_and_axis_labels(df):
    _, ax = ev.plot_volcano(df, title="T", xlabel="LFC", ylabel="P", return_fig=True)
    assert (ax.get_title(), ax.get_xlabel(), ax.get_ylabel()) == ("T", "LFC", "P")


# --------------------------------------------------------------------------- #
# summary text under the plot
# --------------------------------------------------------------------------- #

def _summary(fig):
    return " ".join(t.get_text() for t in fig.texts)


def test_summary_counts(df):
    # after dropna: up_sig = ENS1, ENS4; down_sig = ENS2; 4 others
    fig, _ = ev.plot_volcano(df, return_fig=True)
    assert "up 2; down 1; ns 4" in _summary(fig)


def test_summary_with_zero_lfc_threshold(df):
    fig, _ = ev.plot_volcano(df, x_thr=0, return_fig=True)
    assert "up" in _summary(fig)


def test_summary_goes_on_custom_ax_figure(df):
    fig1, ax1 = plt.subplots()
    fig2 = plt.figure()                        # now the "current" figure
    ev.plot_volcano(df, custom_ax=ax1)
    assert "up" in _summary(fig1)
    assert _summary(fig2) == ""


def test_no_summary_without_thresholds(df):
    fig, _ = ev.plot_volcano(df, x_thr=None, y_thr=None, return_fig=True)
    assert _summary(fig) == ""
