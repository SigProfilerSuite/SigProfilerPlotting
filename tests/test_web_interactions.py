"""End-to-end tests of the local dashboard in a real browser.

Each test starts the real server on a free port, opens the page in headless
Chromium through Playwright, uploads a synthetic matrix through the file
picker and drives the page as a user would.

Needs the test-only dependency Playwright and a browser for it:

    pip install playwright
    playwright install chromium

If Playwright's own Chromium is missing, an installed Google Chrome is used
instead; with neither, these tests are skipped.
"""

import random
import threading
from contextlib import contextmanager

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from sigProfilerPlotting.web.matrix import REFERENCE_CONTEXTS

PAGE = 15          # plots drawn before the reader asks for more
PAGE_4608 = 4

# SBS4608 panel checkboxes, and the layout axes each panel owns.
PANELS = {
    "transcription_bias": "#opt-bias",
    "total_heatmap": "#opt-total",
    "detailed_heatmaps": "#opt-detailed",
}
PANEL_AXES = {
    "transcription_bias": ["yaxis5", "yaxis6", "xaxis2", "xaxis3"],
    "total_heatmap": ["yaxis4"],
    "detailed_heatmaps": ["yaxis2", "yaxis3"],
}
VARIANTS = ("default", "strand_bias", "genic_bias")

XSS_NAMES = [
    "<img src=x onerror=window.__pwned=1>",
    "<svg onload=window.__pwned=2>",
    "</title><script>window.__pwned=3</script>",
]


def _labels(context):
    return list(REFERENCE_CONTEXTS[context])


def _write_matrix(path, context, names, low=5, high=60, cell=None):
    rng = random.Random(11)
    rows = ["MutationType\t" + "\t".join(names)]
    for index, label in enumerate(_labels(context)):
        values = [str(rng.randint(low, high)) for _ in names]
        if cell is not None and index == 0:
            values[0] = cell
        rows.append(label + "\t" + "\t".join(values))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _launch(playwright):
    try:
        return playwright.chromium.launch()
    except Exception as first:
        try:
            return playwright.chromium.launch(channel="chrome")
        except Exception:
            pytest.skip(f"no browser for Playwright; run `playwright install chromium` ({first})")


@contextmanager
def _dashboard():
    """The real server on a free port, and a page open on it."""
    from sigProfilerPlotting.web.server import make_server

    server = make_server("127.0.0.1", 0)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        with sync_api.sync_playwright() as playwright:
            browser = _launch(playwright)
            try:
                page = browser.new_page(viewport={"width": 1400, "height": 900},
                                        accept_downloads=True)
                errors = []
                page.on("pageerror", lambda exc: errors.append(str(exc)))
                page.goto(f"http://127.0.0.1:{server.server_address[1]}/")
                yield page
                assert not errors, errors
            finally:
                browser.close()
    finally:
        server.shutdown()
        server.server_close()


def _upload(page, path, context, ready):
    page.set_input_files("#matrixFile", str(path))
    page.wait_for_function("c => document.getElementById('matrixSummary').textContent.startsWith(c)",
                           arg=context)
    _wait_ready(page, ready)


def _wait_ready(page, n):
    page.wait_for_function(
        "n => document.querySelectorAll('.card[data-state=ready]').length >= n", arg=n)


def _plot(page, index=0):
    return page.locator(".card .js-plotly-plot").nth(index)


def _x_range(plot):
    return plot.evaluate("gd => gd._fullLayout.xaxis.range.slice()")


def _drag(page, plot, start, end, steps=12):
    box = plot.locator(".nsewdrag").first.bounding_box()
    page.mouse.move(box["x"] + box["width"] * start[0], box["y"] + box["height"] * start[1])
    page.mouse.down()
    page.mouse.move(box["x"] + box["width"] * end[0], box["y"] + box["height"] * end[1],
                    steps=steps)
    page.mouse.up()


def _sbs4608_page(page, tmp_path):
    """One SBS4608 sample drawn, and the figures the server sent for it."""
    path = _write_matrix(tmp_path / "sbs4608.txt", "SBS4608", ["S0"], low=0, high=20)
    _upload(page, path, "SBS4608", 1)
    return page.evaluate(
        "() => fetch('/api/render?sample=S0&unit=counts').then(r => r.json())")["variants"]


def _height(plot):
    return plot.evaluate("gd => gd.layout.height")


def _wait_height(page, test, full):
    page.wait_for_function(
        f"full => document.querySelector('.card .js-plotly-plot').layout.height {test} full",
        arg=full)


def _reference_compact_plot():
    """plotInteractive's own compactPlot, from the page it exports."""
    from sigProfilerPlotting import sbs4608_interactive

    page = sbs4608_interactive._SBS4608_PAGE_TEMPLATE.replace("{{", "{").replace("}}", "}")
    start = page.index("const SECTION_AXES")
    end = page.index("function applyPlotOptions")
    return f"(() => {{ {page[start:end]}\n return compactPlot; }})()"


# What a reader can see change when panels are hidden: which traces show,
# where every panel sits, how tall the figure is, and the legend.
_PROJECT = """(traces, layout) => ({
  visible: traces.map(t => t.visible === undefined ? null : t.visible),
  domains: Object.keys(layout).filter(k => /^[xy]axis\\d*$/.test(k)).sort()
    .map(k => [k, layout[k].domain || null, layout[k].visible === undefined ? null : layout[k].visible]),
  height: layout.height,
  showlegend: layout.showlegend === undefined ? null : layout.showlegend,
  legend: layout.legend ? [layout.legend.orientation || null, layout.legend.x ?? null,
    layout.legend.y ?? null, layout.legend.xanchor || null, layout.legend.yanchor || null] : null,
  annotations: (layout.annotations || []).map(a => [a.y, a.visible === undefined ? null : a.visible]),
})"""


def _render_requests(page):
    seen = []
    page.on("request", lambda r: "/api/render" in r.url and seen.append(r.url))
    return seen


# ------------------------------------------------------------------ loading

def test_only_the_first_page_of_plots_is_fetched(tmp_path):
    names = [f"PD{i:04d}a" for i in range(40)]
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", names)
    with _dashboard() as page:
        requests = _render_requests(page)
        _upload(page, path, "SBS96", PAGE)
        page.wait_for_timeout(500)
        assert page.locator(".card").count() == PAGE
        assert len(requests) == PAGE
        assert page.text_content("#moreText") == f"Showing {PAGE} of 40 selected."
        assert page.locator("#sampleList li").count() == 40

        page.click("#loadMore")
        _wait_ready(page, 2 * PAGE)
        assert page.locator(".card").count() == 2 * PAGE
        assert len(requests) == 2 * PAGE
        titles = page.locator(".card-title").all_text_contents()
        assert titles == names[:2 * PAGE]


def test_two_near_simultaneous_load_more_triggers_only_advance_one_page(tmp_path):
    # #sentinel sits immediately after the "Load more" button in the DOM, so
    # scrolling the button into view -- a real click, or Playwright's own
    # actionability check -- can itself bring the sentinel into the
    # IntersectionObserver's rootMargin and fire loadMore() a second time via
    # the auto-scroll path, on top of the button's own click handler. Call
    # loadMore() twice back-to-back (simulating that race directly, rather
    # than depending on exact viewport/font metrics to reproduce it) and
    # confirm only one page's worth of new cards actually appears.
    names = [f"PD{i:04d}a" for i in range(40)]
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", names)
    with _dashboard() as page:
        requests = _render_requests(page)
        _upload(page, path, "SBS96", PAGE)
        page.evaluate("loadMore(); loadMore();")
        _wait_ready(page, 2 * PAGE)
        page.wait_for_timeout(200)
        assert page.locator(".card").count() == 2 * PAGE
        assert len(requests) == 2 * PAGE


def test_scrolling_to_the_end_fetches_the_next_page(tmp_path):
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", [f"S{i}" for i in range(20)])
    with _dashboard() as page:
        requests = _render_requests(page)
        _upload(page, path, "SBS96", PAGE)
        page.locator("#sentinel").scroll_into_view_if_needed()
        _wait_ready(page, 20)
        assert len(requests) == 20
        assert page.is_hidden("#more")


@pytest.mark.parametrize("context", ["ID28", "DBS1248", "SBS384"])
def test_contexts_without_a_dedicated_renderer_still_draw(tmp_path, context):
    path = _write_matrix(tmp_path / f"{context}.txt", context, ["A", "B"])
    with _dashboard() as page:
        _upload(page, path, context, 2)
        assert page.is_hidden("#variantPanel")
        assert page.is_visible("#unitPanel")
        assert _plot(page).evaluate("gd => gd.data.length") > 0


def test_an_invalid_matrix_shows_the_servers_message(tmp_path):
    path = tmp_path / "bad.txt"
    path.write_text("MutationType\tA\tA\n" + "".join(f"{label}\t1\t2\n" for label in _labels("SBS96")))
    with _dashboard() as page:
        page.set_input_files("#matrixFile", str(path))
        page.wait_for_selector("#matrixError:not(.hidden)")
        assert "Duplicate sample/signature column 'A'" in page.text_content("#matrixError")
        assert page.locator(".card").count() == 0


# -------------------------------------------------------------- the list

def test_select_all_clear_all_and_search(tmp_path):
    names = [f"PD{i:04d}a" for i in range(40)]
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", names)
    with _dashboard() as page:
        _upload(page, path, "SBS96", PAGE)
        assert page.text_content("#selectionCount") == "40 of 40 selected"

        page.click("#clearAll")
        assert page.text_content("#selectionCount") == "0 of 40 selected"
        assert page.locator(".card").count() == 0
        assert page.is_visible("#emptyState")

        page.fill("#sampleSearch", "pd001")
        visible = page.locator("#sampleList li:not([hidden])")
        assert visible.count() == 10
        assert page.text_content("#selectAll") == "Select 10 shown"
        page.click("#selectAll")
        assert page.text_content("#selectionCount") == "10 of 40 selected"
        _wait_ready(page, 10)
        assert page.locator(".card-title").all_text_contents() == names[10:20]

        page.fill("#sampleSearch", "no such sample")
        assert page.is_visible("#noMatch")

        page.fill("#sampleSearch", "")
        page.click("#sampleList li:nth-child(1) input")
        assert page.text_content("#selectionCount") == "11 of 40 selected"
        _wait_ready(page, 11)
        assert page.locator(".card-title").first.text_content() == names[0]

        page.click("#selectAll")
        assert page.text_content("#selectionCount") == "40 of 40 selected"


def test_more_than_100_samples_defaults_to_just_the_first(tmp_path):
    # A "select everything" default makes both the sidebar and the plot
    # area unusable before the reader has looked at anything on a matrix
    # with hundreds of columns -- past the threshold, start with only the
    # first sample selected; the reader can still select/search for more.
    names = [f"PD{i:04d}a" for i in range(101)]
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", names)
    with _dashboard() as page:
        _upload(page, path, "SBS96", 1)
        assert page.text_content("#selectionCount") == "1 of 101 selected"
        assert page.locator(".card-title").all_text_contents() == [names[0]]
        assert page.is_checked("#sampleList li:nth-child(1) input")
        assert not page.is_checked("#sampleList li:nth-child(2) input")

        # Still fully selectable/searchable from there.
        page.click("#selectAll")
        assert page.text_content("#selectionCount") == "101 of 101 selected"


def test_exactly_100_samples_still_selects_everything(tmp_path):
    names = [f"PD{i:04d}a" for i in range(100)]
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", names)
    with _dashboard() as page:
        _upload(page, path, "SBS96", PAGE)
        assert page.text_content("#selectionCount") == "100 of 100 selected"


# ---------------------------------------------------- plot interactions

def test_hovering_a_bar_shows_its_value(tmp_path):
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", ["Hover_me"], low=40, high=60)
    with _dashboard() as page:
        _upload(page, path, "SBS96", 1)
        bar = _plot(page).locator("g.point path").nth(10)
        box = bar.bounding_box()
        page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
        tip = page.locator(".hoverlayer .hovertext").first
        tip.wait_for(state="visible")
        assert "Hover_me" in tip.text_content()


def test_zoom_pan_and_reset(tmp_path):
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", ["S1"])
    with _dashboard() as page:
        _upload(page, path, "SBS96", 1)
        plot = _plot(page)
        full = _x_range(plot)

        _drag(page, plot, (0.2, 0.5), (0.5, 0.6))           # box zoom, the default
        zoomed = _x_range(plot)
        assert zoomed[1] - zoomed[0] < 0.5 * (full[1] - full[0])

        plot.hover()
        plot.locator('.modebar-btn[data-attr="dragmode"][data-val="pan"]').click()
        _drag(page, plot, (0.6, 0.5), (0.3, 0.5))
        panned = _x_range(plot)
        assert panned[0] > zoomed[0]
        assert panned[1] - panned[0] == pytest.approx(zoomed[1] - zoomed[0], rel=1e-6)

        plot.hover()
        plot.locator('.modebar-btn[data-attr="zoom"][data-val="reset"]').click()
        assert _x_range(plot) == pytest.approx(full)


def test_the_export_button_downloads_an_image(tmp_path):
    path = _write_matrix(tmp_path / "sbs96.txt", "SBS96", ["S1"])
    with _dashboard() as page:
        _upload(page, path, "SBS96", 1)
        plot = _plot(page)
        plot.hover()
        # The server's config is used as sent: no Plotly logo, SVG export.
        assert plot.locator(".modebar-btn--logo").count() == 0
        with page.expect_download() as info:
            plot.locator('.modebar-btn[data-title^="Download plot"]').click()
        download = info.value
        assert download.suggested_filename.endswith(".svg")
        saved = tmp_path / download.suggested_filename
        download.save_as(saved)
        assert saved.read_text(encoding="utf-8").lstrip().startswith("<svg")


@pytest.mark.parametrize("context", ["SBS96", "DBS78"])
@pytest.mark.parametrize("unit, whole", [("percent", 100.0), ("proportion", 1.0)])
def test_a_unit_redraws_each_sample_as_shares_of_its_total(tmp_path, context, unit, whole):
    # Waits on the values themselves, not the axis title: each renderer is
    # free to word its own title.
    path = _write_matrix(tmp_path / f"{context}.txt", context, ["S1", "S2"])
    with _dashboard() as page:
        _upload(page, path, context, 2)
        total = "gd => gd.data.reduce((s, t) => s + t.y.reduce((a, b) => a + b, 0), 0)"
        assert _plot(page).evaluate(total) > 100
        page.check(f'input[name="unit"][value="{unit}"]', force=True)
        page.wait_for_function(
            "whole => [...document.querySelectorAll('.card .js-plotly-plot')]"
            f".every(gd => Math.abs(({total})(gd) - whole) < 1e-6 * whole)",
            arg=whole)
        for index in range(2):
            assert _plot(page, index).evaluate(total) == pytest.approx(whole)


def test_sbs4608_unit_toggle_changes_the_plot(tmp_path):
    path = _write_matrix(tmp_path / "sbs4608.txt", "SBS4608", ["S0"], low=1, high=20)
    titles = {
        "counts": "Number of Single Base Substitutions",
        "percent": "Percentage of Single Base Substitutions",
        "proportion": "Proportion of Single Base Substitutions",
    }
    with _dashboard() as page:
        _upload(page, path, "SBS4608", 1)
        title = "gd => gd.layout.yaxis.title.text"
        assert _plot(page).evaluate(title) == titles["counts"]  # the default unit

        for unit in ("percent", "proportion", "counts"):
            page.check(f'input[name="unit"][value="{unit}"]', force=True)
            page.wait_for_function(
                "want => { const gd = document.querySelectorAll('.card .js-plotly-plot')[0];"
                " return gd.layout.yaxis.title.text === want; }",
                arg=titles[unit])


def test_sbs1536_has_no_variant_or_transcription_bias_controls(tmp_path):
    sbs1536_path = _write_matrix(tmp_path / "sbs1536.txt", "SBS1536", ["S0"], low=1, high=20)
    sbs4608_path = _write_matrix(tmp_path / "sbs4608.txt", "SBS4608", ["S0"], low=1, high=20)
    with _dashboard() as page:
        _upload(page, sbs1536_path, "SBS1536", 1)
        assert page.is_visible("#variantPanel")
        assert page.is_visible("#unitPanel")
        assert page.is_hidden("#topPanelGroup")
        assert page.is_hidden("#biasCheckboxRow")
        assert page.is_visible("#opt-total")
        assert page.is_visible("#opt-detailed")

        plot = _plot(page)
        shape = "gd => gd.data.map(t => [t.type, t.yaxis || 'y'])"
        before = plot.evaluate(shape)
        assert ["bar", "y"] in before
        assert ["heatmap", "y4"] in before

        page.uncheck("#opt-detailed", force=True)
        page.wait_for_function(
            "() => { const gd = document.querySelectorAll('.card .js-plotly-plot')[0];"
            " return !gd.data.some(t => (t.yaxis === 'y2' || t.yaxis === 'y3') && t.visible !== false); }"
        )

        # Switching to an SBS4608 matrix afterward must restore the controls
        # this context hides -- nothing should carry over as stale state.
        _upload(page, sbs4608_path, "SBS4608", 1)
        assert page.is_visible("#topPanelGroup")
        assert page.is_visible("#biasCheckboxRow")
        assert page.is_checked("#opt-bias")
        assert page.is_checked('input[name="variant"][value="default"]')


def test_sbs4608_switches_the_whole_figure(tmp_path):
    path = _write_matrix(tmp_path / "sbs4608.txt", "SBS4608", [f"S{i}" for i in range(6)],
                         low=0, high=20)
    with _dashboard() as page:
        requests = _render_requests(page)
        _upload(page, path, "SBS4608", PAGE_4608)
        assert len(requests) == PAGE_4608
        assert page.is_visible("#variantPanel")
        assert page.is_visible("#unitPanel")

        served = page.evaluate(
            "() => fetch('/api/render?sample=S0&unit=counts').then(r => r.json())")
        plot = _plot(page)
        shape = ("gd => ({traces: gd.data.map(t => [t.type, t.name || null]),"
                 " height: gd.layout.height})")

        def expected(variant):
            v = served["variants"][variant]
            return {"traces": [[t["type"], t.get("name") or None] for t in v["traces"]],
                    "height": v["layout"]["height"]}

        assert plot.evaluate(shape) == expected("default")
        for variant in ("strand_bias", "genic_bias", "default"):
            page.check(f'input[name="variant"][value="{variant}"]', force=True)
            page.wait_for_function(
                "want => { const gd = document.querySelectorAll('.card .js-plotly-plot')[0];"
                " return JSON.stringify(gd.data.map(t => t.name || null))"
                " === JSON.stringify(want); }",
                arg=[name for _, name in expected(variant)["traces"]])
            assert plot.evaluate(shape) == expected(variant)

        other = _write_matrix(tmp_path / "sbs96.txt", "SBS96", ["A"])
        _upload(page, other, "SBS96", 1)
        assert page.is_visible("#unitPanel")
        assert page.is_hidden("#variantPanel")


def test_each_sbs4608_panel_can_be_hidden_and_restored(tmp_path):
    with _dashboard() as page:
        variants = _sbs4608_page(page, tmp_path)
        served = variants["default"]
        plot = _plot(page)
        full = _height(plot)
        assert full == served["layout"]["height"]
        for panel, box in PANELS.items():
            page.uncheck(box)
            _wait_height(page, "<", full)
            start, end = served["sections"][panel]["trace_range"]
            assert end > start
            visible = plot.evaluate("gd => gd.data.map(t => t.visible)")
            assert all(v is False for v in visible[start:end]), panel
            assert all(v is not False for i, v in enumerate(visible)
                       if not start <= i < end), panel
            for axis in PANEL_AXES[panel]:
                if axis in served["layout"]:
                    assert plot.evaluate(f"gd => gd.layout.{axis}.visible") is False, axis

            page.check(box)
            _wait_height(page, "===", full)
            assert all(v is not False
                       for v in plot.evaluate("gd => gd.data.map(t => t.visible)")), panel


def test_sbs4608_panels_are_packed_exactly_as_plot_interactive_packs_them(tmp_path):
    """The dashboard's copy of compactPlot against the export's own, for every
    combination of panels in every top-panel view."""
    with _dashboard() as page:
        variants = _sbs4608_page(page, tmp_path)
        page.evaluate(f"v => {{ window.__served = v; window.__reference = {_reference_compact_plot()};"
                      f" window.__project = {_PROJECT}; }}", variants)
        names = list(PANELS)
        for variant in VARIANTS:
            page.check(f'input[name="variant"][value="{variant}"]', force=True)
            for step in range(8):
                gray = step ^ (step >> 1)              # one box changes per step
                toggles = {name: not gray & (1 << i) for i, name in enumerate(names)}
                for name, on in toggles.items():
                    page.set_checked(PANELS[name], on)
                compare = ("([variant, toggles]) => {"
                           " const gd = document.querySelector('.card .js-plotly-plot');"
                           " const ref = window.__reference(window.__served[variant], toggles);"
                           " return [JSON.stringify(window.__project(gd.data, gd.layout)),"
                           "         JSON.stringify(window.__project(ref.traces, ref.layout))]; }")
                try:
                    page.wait_for_function(f"a => {{ const [got, want] = ({compare})(a);"
                                           " return got === want; }", arg=[variant, toggles])
                except sync_api.TimeoutError:
                    got, want = page.evaluate(compare, [variant, toggles])
                    assert got == want, (variant, toggles)


def test_hiding_every_panel_keeps_the_strand_legend_clear_of_the_plot(tmp_path):
    with _dashboard() as page:
        variants = _sbs4608_page(page, tmp_path)
        plot = _plot(page)
        page.check('input[name="variant"][value="strand_bias"]', force=True)
        full = variants["strand_bias"]["layout"]["height"]
        _wait_height(page, "===", full)
        for box in PANELS.values():
            page.uncheck(box)
        page.wait_for_function(
            "() => document.querySelectorAll('.card .js-plotly-plot')[0]"
            ".data.every(t => t.visible !== false || t.xaxis !== 'x' || t.yaxis !== 'y')"
            " && document.querySelectorAll('#opt-bias:checked, #opt-total:checked,"
            " #opt-detailed:checked').length === 0")
        _wait_height(page, "<", full)
        page.wait_for_timeout(300)

        legend = plot.evaluate("gd => gd.layout.legend")
        assert plot.evaluate("gd => gd.layout.showlegend") is not False
        assert (legend["orientation"], legend["xanchor"], legend["yanchor"]) == \
            ("h", "center", "bottom")
        # In the drawn figure: the legend has entries, sits inside the figure,
        # and overlaps neither the bars nor their x labels.
        boxes = plot.evaluate("""gd => {
          const r = el => { const b = el.getBoundingClientRect();
                            return {top: b.top, bottom: b.bottom, left: b.left, right: b.right}; };
          return {
            legend: r(gd.querySelector('g.legend')),
            entries: gd.querySelectorAll('g.legend g.traces').length,
            bars: r(gd.querySelector('.nsewdrag[data-subplot="xy"]')),
            labels: r(gd.querySelector('g.xaxislayer-above')),
            figure: r(gd.querySelector('svg.main-svg')),
          }; }""")

        def apart(a, b):
            return (a["bottom"] <= b["top"] or b["bottom"] <= a["top"]
                    or a["right"] <= b["left"] or b["right"] <= a["left"])

        assert boxes["entries"] > 0
        assert apart(boxes["legend"], boxes["bars"])
        assert apart(boxes["legend"], boxes["labels"])
        assert boxes["legend"]["bottom"] <= boxes["figure"]["bottom"]

        # Another view keeps the panels hidden; the default view has no
        # strand legend to keep, so it shows none.
        page.check('input[name="variant"][value="default"]', force=True)
        page.wait_for_function(
            "() => document.querySelector('.card .js-plotly-plot').layout.showlegend === false")
        start, end = variants["default"]["sections"]["total_heatmap"]["trace_range"]
        visible = plot.evaluate("gd => gd.data.map(t => t.visible)")
        assert all(v is False for v in visible[start:end])


# --------------------------------------------------------------- safety

def test_markup_in_sample_names_is_shown_as_text(tmp_path):
    path = _write_matrix(tmp_path / "xss.txt", "SBS96", XSS_NAMES)
    with _dashboard() as page:
        dialogs = []
        page.on("dialog", lambda d: (dialogs.append(d.message), d.dismiss()))
        _upload(page, path, "SBS96", len(XSS_NAMES))
        page.wait_for_timeout(300)
        assert page.locator("#sampleList span").all_text_contents() == XSS_NAMES
        assert page.locator(".card-title").all_text_contents() == XSS_NAMES
        assert page.evaluate("window.__pwned") is None
        assert page.evaluate("document.querySelectorAll('img, [onerror], [onload]').length") == 0
        assert page.locator("script").count() == 2          # plotly.js and app.js
        assert dialogs == []

        page.fill("#sampleSearch", "<svg")
        assert page.locator("#sampleList li:not([hidden]) span").all_text_contents() == [XSS_NAMES[1]]


def test_markup_in_a_server_error_is_shown_as_text(tmp_path):
    bad = "<img src=x onerror=window.__pwned=4>"
    path = _write_matrix(tmp_path / "bad.txt", "SBS96", ["A"], cell=bad)
    with _dashboard() as page:
        page.set_input_files("#matrixFile", str(path))
        page.wait_for_selector("#matrixError:not(.hidden)")
        assert bad in page.text_content("#matrixError")
        assert page.evaluate("window.__pwned") is None
        assert page.evaluate("document.querySelectorAll('img').length") == 0
