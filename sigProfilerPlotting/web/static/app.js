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

const $ = (id) => document.getElementById(id);

const state = {
  context: null,
  samples: [],            // column order, as in the matrix
  items: [],              // {name, lower, li, input}
  selected: new Set(),
  unit: 'counts',
  variant: 'default',
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
  state.selected = new Set(state.samples);
  state.limit = pageSize();

  const n = state.samples.length;
  $('matrixSummary').textContent =
    `${state.context} · ${info.labels_count} mutation contexts · ` +
    `${n} ${n === 1 ? 'sample' : 'samples'}`;
  $('variantPanel').classList.toggle('hidden', state.context !== 'SBS4608');
  // SBS4608 is always drawn as a percentage, so a unit choice would do nothing.
  $('unitPanel').classList.toggle('hidden', state.context === 'SBS4608');
  $('sampleSearch').value = '';
  for (const id of ['sampleSearch', 'selectAll', 'clearAll']) $(id).disabled = false;
  renderList();
  sync();
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
  const stale = () => !card.alive || card.token !== token || generation !== state.generation;
  setStatus(card, card.drawn ? 'Updating' : 'Loading', 'loading');
  enqueue(async () => {
    if (stale()) return null;
    const query = new URLSearchParams({ sample: card.name, unit });
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
  const { traces, layout, config } = structuredClone(figure);
  const plot = card.drawn ? Plotly.react : Plotly.newPlot;
  card.pending = card.pending
    .then(() => (card.alive ? plot(card.plot, traces, layout, config) : null))
    .then(() => {
      if (!card.alive) return;
      card.drawn = true;
      setStatus(card, '', null);
    })
    .catch((e) => setStatus(card, `This plot could not be drawn: ${e.message}`, 'error'));
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

function loadMore() {
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

function onVariant(variant) {
  state.variant = variant;
  for (const card of state.cards.values()) {
    if (card.data && card.data.kind === 'sbs4608') draw(card);
  }
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

  new IntersectionObserver((entries) => {
    state.sentinelVisible = entries.some((e) => e.isIntersecting);
    maybeLoadMore();
  }, { rootMargin: '300px' }).observe($('sentinel'));
}

init();
