'use strict';
// SigProfilerPlotting dashboard: upload a matrix, pick samples, draw each one.
//
// Every string that comes from the user's file or from the server (sample
// names, error messages) reaches the page through textContent or value, never
// innerHTML, so a sample named like markup is shown as the text it is.
//
// Plots are drawn lazily: only the first page of selected samples is fetched,
// and more are fetched when the reader scrolls to the end or asks for them.

const PAGE = 15;
const PAGE_4608 = 4;      // each SBS4608 sample is three full figures, ~0.7 MB
const MAX_PARALLEL = 3;   // render requests in flight at once
const SELECT_ALL_THRESHOLD = 100; // above this many samples/signatures, default to just the first
// Contexts cosmic.collapse_to_sbs96() can reduce to the 96 trinucleotide
// contexts COSMIC's own catalogue uses -- matches sigProfilerPlotting/web/
// renderers.py's COSMIC_COMPATIBLE_CONTEXTS.
const COSMIC_COMPATIBLE_CONTEXTS = new Set(['SBS96', 'SBS288', 'SBS1536', 'SBS4608']);

const $ = (id) => document.getElementById(id);

const state = {
  context: null,
  samples: [],            // column order, as in the matrix
  items: [],              // {name, lower, li, input}
  selected: new Set(),
  unit: 'counts',
  variant: 'default',
  // SBS4608 panels shown, as in plotInteractive: all on to begin with.
  toggles: { transcription_bias: true, total_heatmap: true, detailed_heatmaps: true },
  cosmic: { builds: [], catalogues: [], build: null, enabled: false },
  limit: 0,               // how many selected samples may be drawn
  cards: new Map(),       // sample name -> card
  generation: 0,          // bumped per matrix, so late answers are dropped
  sentinelVisible: false,
};

const pageSize = () => (state.context === 'SBS4608' ? PAGE_4608 : PAGE);

/* ------------------------------------------------------------ requests */

let active = 0;
const queue = [];

function enqueue(task) {
  return new Promise((resolve, reject) => {
    queue.push({ task, resolve, reject });
    pump();
  });
}

function pump() {
  while (active < MAX_PARALLEL && queue.length) {
    const { task, resolve, reject } = queue.shift();
    active++;
    task().then(resolve, reject).finally(() => {
      active--;
      pump();
      maybeLoadMore();
    });
  }
}

const idle = () => active === 0 && queue.length === 0;

async function readJson(response) {
  try {
    return await response.json();
  } catch (e) {
    return { error: `The server answered ${response.status} ${response.statusText} with no details.` };
  }
}

/* -------------------------------------------------------------- upload */

function showError(message) {
  const box = $('matrixError');
  box.textContent = message;
  box.classList.remove('hidden');
}

function clearError() {
  $('matrixError').textContent = '';
  $('matrixError').classList.add('hidden');
}

async function upload(file) {
  clearError();
  $('fileName').textContent = file.name;
  const before = $('matrixSummary').textContent;
  $('matrixSummary').textContent = 'Reading the matrix';
  const form = new FormData();
  form.append('matrix', file, file.name);
  let response;
  let body;
  try {
    response = await fetch('/api/matrix', { method: 'POST', body: form });
    body = await readJson(response);
  } catch (e) {
    $('matrixSummary').textContent = before;
    showError(`The dashboard server could not be reached (${e.message}). Check that it is still running.`);
    return;
  }
  if (!response.ok || !Array.isArray(body.samples)) {
    $('matrixSummary').textContent = before;
    showError(body.error || `The upload failed: the server answered ${response.status}.`);
    return;
  }
  loadMatrix(body);
}

function loadMatrix(info) {
  state.generation++;
  for (const name of [...state.cards.keys()]) dropCard(name);
  state.context = String(info.context);
  state.samples = info.samples.map(String);
  // Past a threshold a "select everything" default makes the sidebar and
  // the plot area equally unusable before the reader has looked at
  // anything -- start with just the first one, selectable up from there.
  state.selected = state.samples.length > SELECT_ALL_THRESHOLD
    ? new Set(state.samples.slice(0, 1))
    : new Set(state.samples);
  state.limit = pageSize();

  const n = state.samples.length;
  $('matrixSummary').textContent =
    `${state.context} · ${info.labels_count} mutation contexts · ` +
    `${n} ${n === 1 ? 'sample' : 'samples'}`;
  // SBS1536 reuses the SBS4608 panel renderer (same "sbs4608" wire kind --
  // a 1536-context matrix is exactly what SBS4608Data collapses to once
  // strand is summed away) but has no strand/transcription data at all, so
  // it only ever has a "default" view and no transcription-bias panel.
  const hasPanels = state.context === 'SBS4608' || state.context === 'SBS1536';
  const hasVariants = state.context === 'SBS4608';
  state.variant = 'default';
  // SBS1536 has no transcription-bias panel to show at all (no checkbox is
  // rendered for it below) -- force that section off so compactPlot()
  // reclaims the vertical space it would otherwise leave blank, and keep
  // state.toggles/the checkboxes' checked state in sync either way so a
  // later matrix in the same session starts from a consistent default
  // rather than whatever was left over from the one before it.
  state.toggles = { transcription_bias: hasVariants, total_heatmap: true, detailed_heatmaps: true };
  $('opt-bias').checked = state.toggles.transcription_bias;
  $('opt-total').checked = state.toggles.total_heatmap;
  $('opt-detailed').checked = state.toggles.detailed_heatmaps;
  $('variantPanel').classList.toggle('hidden', !hasPanels);
  $('variantPanelTitle').textContent = `${state.context} view`;
  $('topPanelGroup').classList.toggle('hidden', !hasVariants);
  $('biasCheckboxRow').classList.toggle('hidden', !hasVariants);
  if (!hasVariants) {
    const defaultRadio = document.querySelector('input[name="variant"][value="default"]');
    if (defaultRadio) defaultRadio.checked = true;
  }
  // COSMIC comparison only makes sense for the contexts collapse_to_sbs96()
  // can reduce to the 96 COSMIC trinucleotide contexts; reset per matrix so
  // a comparison left on for one file doesn't silently carry into the next.
  state.cosmic.enabled = false;
  $('cosmicEnabled').checked = false;
  $('cosmicPanel').classList.toggle('hidden', !COSMIC_COMPATIBLE_CONTEXTS.has(state.context));
  $('sampleSearch').value = '';
  for (const id of ['sampleSearch', 'selectAll', 'clearAll']) $(id).disabled = false;
  renderList();
  syncCosmicControls();
  sync();
}

/* ------------------------------------------------------------- COSMIC */

function cosmicEligible() {
  return (
    state.cosmic.builds.length > 0 &&
    COSMIC_COMPATIBLE_CONTEXTS.has(state.context) &&
    state.selected.size === 1
  );
}

// Keeps the Compare-with-COSMIC controls consistent with the current
// selection -- called after every selection change (via updateCount()) and
// on each matrix load, so the checkbox can never stay checked once it is
// no longer meaningful (more than one, or zero, samples selected).
function syncCosmicControls() {
  const eligible = cosmicEligible();
  if (!eligible && state.cosmic.enabled) {
    state.cosmic.enabled = false;
    $('cosmicEnabled').checked = false;
    // A card already drawn in compare mode (e.g. the reader just selected a
    // second sample, making the comparison no longer meaningful) does not
    // otherwise get re-fetched just because the selection count changed --
    // without this it would keep showing the now-stale comparison.
    for (const card of state.cards.values()) loadCard(card);
  }
  $('cosmicEnabled').disabled = !eligible;
  $('cosmicBuild').disabled = !eligible || !state.cosmic.enabled;
  $('cosmicHint').textContent = state.cosmic.builds.length === 0
    ? 'No COSMIC catalogue is configured. SigProfilerPlotting does not bundle this ' +
      'data -- launch the dashboard with --cosmic-dir pointing at your own ' +
      'downloaded COSMIC SBS reference files (see the docs) to enable this panel.'
    : 'Select exactly one sample or signature to compare it.';
  $('cosmicHint').classList.toggle('hidden', eligible);
}

async function loadCosmicBuilds() {
  const select = $('cosmicBuild');
  try {
    const response = await fetch('/api/cosmic/builds');
    const body = await readJson(response);
    state.cosmic.builds = (response.ok && Array.isArray(body.builds)) ? body.builds : [];
    state.cosmic.catalogues = (response.ok && Array.isArray(body.catalogues)) ? body.catalogues : [];
  } catch (e) {
    state.cosmic.builds = [];
    state.cosmic.catalogues = [];
  }
  select.replaceChildren();
  for (const build of state.cosmic.builds) {
    const option = document.createElement('option');
    option.value = build;
    option.textContent = build;
    const info = state.cosmic.catalogues.find((c) => c.build === build);
    if (info) {
      option.title = `COSMIC v${info.version} -- ${info.signature_count} signatures\n` +
        `${info.source_path}\nsha256:${info.checksum}`;
    }
    select.append(option);
  }
  state.cosmic.build = state.cosmic.builds[0] || null;
}

function onCosmicChange() {
  state.cosmic.enabled = $('cosmicEnabled').checked;
  state.cosmic.build = $('cosmicBuild').value || state.cosmic.build;
  syncCosmicControls();
  for (const card of state.cards.values()) loadCard(card);
}

/* --------------------------------------------------------- sample list */

function renderList() {
  const list = $('sampleList');
  list.replaceChildren();
  const frag = document.createDocumentFragment();
  state.items = state.samples.map((name) => {
    const li = document.createElement('li');
    const label = document.createElement('label');
    const input = document.createElement('input');
    input.type = 'checkbox';
    input.value = name;
    input.checked = state.selected.has(name);
    const text = document.createElement('span');
    text.textContent = name;
    label.append(input, text);
    li.append(label);
    frag.append(li);
    return { name, lower: name.toLowerCase(), li, input };
  });
  list.append(frag);
  applySearch();
}

function shownItems() {
  return state.items.filter((item) => !item.li.hidden);
}

function applySearch() {
  const q = $('sampleSearch').value.trim().toLowerCase();
  let shown = 0;
  for (const item of state.items) {
    item.li.hidden = q !== '' && !item.lower.includes(q);
    if (!item.li.hidden) shown++;
  }
  $('noMatch').classList.toggle('hidden', shown > 0 || state.items.length === 0);
  $('selectAll').textContent = q ? `Select ${shown} shown` : 'Select all';
  $('clearAll').textContent = q ? 'Clear shown' : 'Clear all';
  updateCount();
}

function updateCount() {
  const total = state.samples.length;
  $('selectionCount').textContent = total
    ? `${state.selected.size} of ${total} selected` : '';
  syncCosmicControls();
}

function setSelected(items, on) {
  for (const item of items) {
    item.input.checked = on;
    if (on) state.selected.add(item.name);
    else state.selected.delete(item.name);
  }
  if (!on && state.selected.size === 0) state.limit = pageSize();
  updateCount();
  sync();
}

/* --------------------------------------------------------------- cards */

function makeCard(name) {
  const el = $('cardTemplate').content.firstElementChild.cloneNode(true);
  el.querySelector('.card-title').textContent = name;
  const plot = el.querySelector('.plot');
  plot.setAttribute('aria-label', `Plot of ${name}`);
  el.dataset.state = 'loading';
  return {
    name, el, plot, status: el.querySelector('.card-status'),
    data: null, token: 0, alive: true, drawn: false, pending: Promise.resolve(),
  };
}

function dropCard(name) {
  const card = state.cards.get(name);
  if (!card) return;
  card.alive = false;
  if (card.drawn && window.Plotly) Plotly.purge(card.plot);
  if (window.Plotly) {
    for (const el of card.el.querySelectorAll('.compare-scroll .plot')) Plotly.purge(el);
  }
  card.el.remove();
  state.cards.delete(name);
}

function setStatus(card, text, kind) {
  card.status.textContent = text;
  card.status.classList.toggle('hidden', !text);
  card.status.classList.toggle('error', kind === 'error');
  card.el.dataset.state = kind || (text ? 'loading' : 'ready');
}

function loadCard(card) {
  const token = ++card.token;
  const unit = state.unit;
  const generation = state.generation;
  const compareCosmic = (state.cosmic.enabled && cosmicEligible() && state.cosmic.build) || null;
  const stale = () => !card.alive || card.token !== token || generation !== state.generation;
  setStatus(card, card.drawn ? 'Updating' : 'Loading', 'loading');
  enqueue(async () => {
    if (stale()) return null;
    const params = { sample: card.name, unit };
    if (compareCosmic) params.compare_cosmic = state.cosmic.build;
    const query = new URLSearchParams(params);
    const response = await fetch(`/api/render?${query}`);
    const body = await readJson(response);
    if (!response.ok) throw new Error(body.error || `The server answered ${response.status}.`);
    return body;
  }).then((body) => {
    if (body === null || stale()) return;
    card.data = body;
    draw(card);
  }).catch((e) => {
    if (stale()) return;
    setStatus(card, `This plot could not be drawn: ${e.message}`, 'error');
  });
}

function figureOf(card) {
  const d = card.data;
  if (d.kind === 'standard') return d.figure;
  if (d.kind === 'sbs4608') return d.variants && d.variants[state.variant];
  return null;
}

function draw(card) {
  const figure = figureOf(card);
  if (!figure) {
    setStatus(card, `This page cannot draw a plot of kind “${card.data.kind}”.`, 'error');
    return;
  }
  // A copy each time: Plotly writes into the objects it is given, and a
  // variant drawn again later must start from what the server sent.
  // compactPlot copies the traces and layout itself.
  const { traces, layout, config } = card.data.kind === 'sbs4608'
    ? { ...compactPlot(figure, state.toggles), config: structuredClone(figure.config) }
    : structuredClone(figure);
  const plot = card.drawn ? Plotly.react : Plotly.newPlot;
  card.pending = card.pending
    .then(() => (card.alive ? plot(card.plot, traces, layout, config) : null))
    .then(() => {
      if (!card.alive) return;
      card.drawn = true;
      setStatus(card, '', null);
      drawCompareMatches(card);
    })
    .catch((e) => setStatus(card, `This plot could not be drawn: ${e.message}`, 'error'));
}

function drawCompareMatches(card) {
  const compareEl = card.el.querySelector('.compare');
  const scrollEl = compareEl.querySelector('.compare-scroll');
  // Purge before rebuilding either way: cheap (at most 5 full-size plots,
  // and only ever reachable with exactly one sample selected), and avoids
  // tracking per-match Plotly state across redraws (unit change, a
  // different COSMIC build, or leaving compare mode entirely).
  for (const el of scrollEl.querySelectorAll('.plot')) {
    if (window.Plotly) Plotly.purge(el);
  }
  const compare = card.data.cosmic_compare;
  if (!compare) {
    compareEl.classList.add('hidden');
    scrollEl.replaceChildren();
    return;
  }
  compareEl.classList.remove('hidden');
  const frag = document.createDocumentFragment();
  const toRender = [];
  for (const match of compare.matches) {
    const el = $('compareMatchTemplate').content.firstElementChild.cloneNode(true);
    const plotEl = el.querySelector('.plot');
    plotEl.setAttribute('aria-label', `${match.name}, cosine similarity ${match.score.toFixed(3)}`);
    frag.append(el);
    toRender.push({ plotEl, figure: match.figure });
  }
  scrollEl.replaceChildren(frag);
  for (const { plotEl, figure } of toRender) {
    const { traces, layout, config } = structuredClone(figure);
    Plotly.newPlot(plotEl, traces, layout, config);
  }
  setupCompareNav(compareEl, scrollEl, compare.matches.length);
}

// One match panel is visible at a time (full width, scroll-snapped); the
// prev/next buttons and the "N of 5" readout stay in sync with whichever
// one a scroll (trackpad, wheel, or the buttons themselves) brings to rest.
function setupCompareNav(compareEl, scrollEl, count) {
  const positionEl = compareEl.querySelector('.compare-position');
  const updatePosition = () => {
    const width = Math.max(scrollEl.clientWidth, 1);
    const index = Math.min(Math.round(scrollEl.scrollLeft / width), count - 1);
    positionEl.textContent = count ? `${index + 1} of ${count}` : '';
  };
  scrollEl.onscroll = updatePosition;
  compareEl.querySelector('.compare-prev').onclick = () =>
    scrollEl.scrollBy({ left: -scrollEl.clientWidth, behavior: 'smooth' });
  compareEl.querySelector('.compare-next').onclick = () =>
    scrollEl.scrollBy({ left: scrollEl.clientWidth, behavior: 'smooth' });
  updatePosition();
}

/* ------------------------------------------------- SBS4608 panels */

// Ported unchanged from plotInteractive's page (sbs4608_interactive.py), so
// the dashboard hides and repacks panels exactly as that export does.
const SECTION_AXES = {
  total_heatmap: ['yaxis4'],
  detailed_heatmaps: ['yaxis2', 'yaxis3'],
  transcription_bias: ['yaxis5', 'yaxis6', 'xaxis2', 'xaxis3'],
};

function compactPlot(variant, toggles) {
  // Start from immutable full-layout data so repeated toggles cannot drift.
  const traces = JSON.parse(JSON.stringify(variant.traces));
  const layout = JSON.parse(JSON.stringify(variant.layout));
  const removed = [];
  for (const [key, sec] of Object.entries(variant.sections)) {
    const on = toggles[key];
    const [start, end] = sec.trace_range;
    if (!on && end > start) removed.push([...sec.paper_band]);
    for (let i = start; i < end; i++) traces[i].visible = on;
    SECTION_AXES[key].forEach(axis => {
      if (layout[axis]) layout[axis].visible = on && end > start;
    });
    (sec.annotation_indices || []).forEach(index => {
      layout.annotations[index].visible = on;
    });
  }
  const margin = layout.margin || {};
  const verticalMargins = (margin.t === undefined ? 100 : margin.t) +
    (margin.b === undefined ? 80 : margin.b);
  const originalPlotPixels = layout.height - verticalMargins;
  const detailedRange = variant.sections.detailed_heatmaps.trace_range;
  const anyHeatmap = toggles.total_heatmap ||
    (toggles.detailed_heatmaps && detailedRange[1] > detailedRange[0]);
  const needsMainLegend = !toggles.transcription_bias &&
    variant.main_plot_type !== 'default';
  if (needsMainLegend) {
    // Reserve a footer sized in PIXELS (legend height + the x-axis tick-label
    // gutter below the bars), not a fixed paper fraction -- a fixed fraction
    // (e.g. 0.055) is only big enough at the full figure height; once enough
    // other panels are hidden and the figure shrinks, that same fraction is
    // too few pixels and the legend rides up over the tick labels.
    const legendPixels = 50;
    const labelGutterPixels = 40;
    const biasBand = removed.find(band => band[0] === 0);
    if (biasBand) {
      const footerFraction = (legendPixels + labelGutterPixels) / originalPlotPixels;
      biasBand[0] = Math.min(biasBand[1], footerFraction);
    }
    layout.legend = {
      ...(layout.legend || {}),
      orientation: 'h',
      x: 0.5,
      xanchor: 'center',
      y: 0.012,
      yanchor: 'bottom',
    };
  } else if (!toggles.transcription_bias) {
    layout.showlegend = false;
  }
  if (toggles.transcription_bias && !anyHeatmap && removed.length) {
    // Keep a small label gutter, not a missing panel, below the main x labels.
    const topBand = removed.reduce((a, b) => a[1] > b[1] ? a : b);
    topBand[1] -= 0.04;
  }
  const remaining = 1 - removed.reduce((sum, band) => sum + band[1] - band[0], 0);
  const mapY = y => (y - removed.reduce((sum, [lo, hi]) =>
    sum + Math.max(0, Math.min(y, hi) - lo), 0)) / remaining;
  for (const [key, axis] of Object.entries(layout)) {
    if (!/^yaxis\d*$/.test(key) || !axis.domain) continue;
    if (axis.visible === false) axis.domain = [0, 0.001];
    else axis.domain = axis.domain.map(mapY);
  }
  (layout.annotations || []).forEach(annotation => {
    if (annotation.yref === 'paper') annotation.y = mapY(annotation.y);
  });
  (layout.shapes || []).forEach(shape => {
    if (shape.yref === 'paper') {
      shape.y0 = mapY(shape.y0);
      shape.y1 = mapY(shape.y1);
    }
  });
  traces.forEach(trace => {
    if (trace.visible === false || !trace.colorbar) return;
    if (trace.colorbar.y !== undefined) trace.colorbar.y = mapY(trace.colorbar.y);
    if (trace.colorbar.len !== undefined && trace.colorbar.lenmode !== 'pixels')
      trace.colorbar.len /= remaining;
  });
  if (layout.legend && layout.legend.y !== undefined)
    layout.legend.y = mapY(layout.legend.y);
  layout.height = Math.round(verticalMargins + (layout.height - verticalMargins) * remaining);
  return {traces, layout};
}

/* ------------------------------------------------------ what is drawn */

function sync() {
  const order = state.samples.filter((name) => state.selected.has(name));
  const target = order.slice(0, state.limit);
  const keep = new Set(target);
  for (const name of [...state.cards.keys()]) if (!keep.has(name)) dropCard(name);

  const plots = $('plots');
  let prev = null;
  for (const name of target) {
    let card = state.cards.get(name);
    if (!card) {
      card = makeCard(name);
      state.cards.set(name, card);
      loadCard(card);
    }
    const spot = prev ? prev.el.nextSibling : plots.firstChild;
    if (card.el !== spot) plots.insertBefore(card.el, spot);
    prev = card;
  }

  const rest = order.length - target.length;
  $('more').classList.toggle('hidden', rest <= 0);
  if (rest > 0) {
    $('moreText').textContent = `Showing ${target.length} of ${order.length} selected.`;
    $('loadMore').textContent = `Load ${Math.min(rest, pageSize())} more`;
  }
  const empty = $('emptyState');
  empty.classList.toggle('hidden', target.length > 0);
  if (state.samples.length && !target.length) {
    empty.querySelector('h2').textContent = 'Nothing selected';
    empty.querySelector('p').textContent = 'Select samples or signatures in the list to draw them.';
  }
}

let loadMoreCooldownUntil = 0;

function loadMore() {
  // #sentinel sits immediately after the "Load more" button in the DOM, so
  // scrolling the button into view (as a real click, or Playwright's own
  // actionability check, does) can itself bring the sentinel into the
  // IntersectionObserver's rootMargin and fire the auto-load path too --
  // one logical "load more" action must not advance state.limit twice just
  // because both triggers landed within the same moment.
  const now = Date.now();
  if (now < loadMoreCooldownUntil) return;
  loadMoreCooldownUntil = now + 1000;
  state.limit += pageSize();
  sync();
}

function maybeLoadMore() {
  // Only once the current page has arrived, so a visible end of the page does
  // not fetch every selected sample in one go.
  if (!state.sentinelVisible || !idle()) return;
  if ($('more').classList.contains('hidden')) return;
  loadMore();
}

/* -------------------------------------------------------------- wiring */

function onUnit(unit) {
  state.unit = unit;
  for (const card of state.cards.values()) loadCard(card);
}

function redrawSbs4608() {
  for (const card of state.cards.values()) {
    if (card.data && card.data.kind === 'sbs4608') draw(card);
  }
}

function onVariant(variant) {
  state.variant = variant;
  redrawSbs4608();
}

function onToggle(section, on) {
  state.toggles[section] = on;
  redrawSbs4608();
}

async function quit() {
  if (!window.confirm('Stop the SigProfilerPlotting dashboard server? This closes the session for every open tab.')) {
    return;
  }
  const button = $('quitButton');
  button.disabled = true;
  button.textContent = 'Stopping...';
  try {
    let token = '';
    try {
      const identity = await (await fetch('/api/identity')).json();
      token = identity.token || '';
    } catch (e) {
      // If /api/identity itself is unreachable the shutdown call below will
      // fail the same way -- fall through with no token rather than block
      // the user from even trying to stop the server.
    }
    await fetch('/api/shutdown', { method: 'POST', headers: { 'X-Dashboard-Token': token } });
  } catch (e) {
    // The server closes its socket as part of shutting down, so the fetch
    // itself failing (rather than resolving) is the expected success case.
  }
  document.body.innerHTML = '';
  document.title = 'SigProfilerPlotting dashboard (stopped)';
  const message = document.createElement('p');
  message.className = 'stopped-message';
  message.textContent = 'The dashboard server has stopped. You can close this tab.';
  document.body.appendChild(message);
}

function init() {
  if (!window.Plotly) {
    showError('The plotting library did not load. Reload the page; if this persists, restart the dashboard.');
  }
  $('matrixFile').addEventListener('change', (ev) => {
    const file = ev.target.files && ev.target.files[0];
    ev.target.value = '';           // the same file can be chosen again after a fix
    if (file) upload(file);
  });
  $('sampleList').addEventListener('change', (ev) => {
    const input = ev.target;
    if (input.type !== 'checkbox') return;
    if (input.checked) state.selected.add(input.value);
    else state.selected.delete(input.value);
    updateCount();
    sync();
  });
  $('sampleSearch').addEventListener('input', applySearch);
  $('selectAll').addEventListener('click', () => setSelected(shownItems(), true));
  $('clearAll').addEventListener('click', () => setSelected(shownItems(), false));
  $('loadMore').addEventListener('click', loadMore);
  document.querySelectorAll('input[name="unit"]').forEach((r) =>
    r.addEventListener('change', () => r.checked && onUnit(r.value)));
  document.querySelectorAll('input[name="variant"]').forEach((r) =>
    r.addEventListener('change', () => r.checked && onVariant(r.value)));
  document.querySelectorAll('input[name="section"]').forEach((c) =>
    c.addEventListener('change', () => onToggle(c.value, c.checked)));
  $('cosmicEnabled').addEventListener('change', onCosmicChange);
  $('cosmicBuild').addEventListener('change', onCosmicChange);
  $('quitButton').addEventListener('click', quit);

  loadCosmicBuilds();

  new IntersectionObserver((entries) => {
    state.sentinelVisible = entries.some((e) => e.isIntersecting);
    maybeLoadMore();
  }, { rootMargin: '300px' }).observe($('sentinel'));
}

init();
