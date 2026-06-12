// File: js/plot.js
import { state, subscribe } from './state.js';
import { parseRational } from './rational.js';

const X_SAMPLES_MIN = 1200;
const SPAN_MIN = 5;
const SPAN_MUL = 1.5;
const Y_MUL = 1.3;
const Y_FALLBACK_MUL = 0.5;
const ZOOM_EXTENT = [1e-10, 1e12];

const COLORS = {
  free_generator:    '#d62728',
  combination:       '#ff9e3d',
  torsion_generator: '#1f77b4',
  torsion_other:     '#88c0d0',
  mixed:             '#9467bd',
};

const SVG_SYMBOL = {
  free_generator:    { type: d3.symbolCircle,   size: 140 },
  combination:       { type: d3.symbolCross,    size:  60 },
  torsion_generator: { type: d3.symbolDiamond,  size: 110 },
  torsion_other:     { type: d3.symbolSquare,   size:  90 },
  mixed:             { type: d3.symbolTriangle, size: 100 },
};

// Cardano 3 分岐に対応する色 (参考コードの tab:blue/orange/green)
const BAND_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c'];
const BAND_OPACITY = 0.15;

let svg, gAxes, gBands, gCurve, gPoints;
let xScale = null, yScale = null;
let baseDomain = { x: [-1, 1], y: [-1, 1] };
let symlogConst = { x: 1, y: 1 };
let zoomBehavior;

export function initPlot() {
  const el = document.getElementById('plot');
  svg = d3.select(el);
  // 描画順: 軸 → 帯 → 曲線 → 点
  gAxes   = svg.append('g').attr('class', 'ax');
  gBands  = svg.append('g').attr('class', 'bd');
  gCurve  = svg.append('g').attr('class', 'cv');
  gPoints = svg.append('g').attr('class', 'pt');

  state.set({ zoom: d3.zoomIdentity });

  zoomBehavior = d3.zoom()
    .scaleExtent(ZOOM_EXTENT)
    .on('zoom', (e) => state.set({ zoom: e.transform }));

  svg.call(zoomBehavior)
     .on('dblclick.zoom', null);

  window.addEventListener('resize', renderPlot);
}

export function resetView() {
  if (!state.get('currentCurve')) return;
  svg.call(zoomBehavior.transform, d3.zoomIdentity);
  state.set({ zoom: d3.zoomIdentity });
  computeBaseDomain(state.get('currentCurve'));
  rebuildScales();
  renderPlot();
}

function computeBaseDomain(curve) {
  if (!curve || !curve.a_invariants) return;
  const a = curve.a_invariants.map(Number);
  const [a1, a2, a3, a4, a6] = a;
  const pts = curve.points || [];
  const xs = pts.map(p => parseRational(p.xy[0])).filter(Number.isFinite);
  const ys = pts.map(p => parseRational(p.xy[1])).filter(Number.isFinite);

  let xmin, xmax;
  if (xs.length > 0) {
    const finite = xs.filter(v => Math.abs(v) < 1e15);
    if (finite.length > 0) {
      finite.sort((p, q) => p - q);
      const med = finite[Math.floor(finite.length / 2)];
      const dev = finite.map(x => Math.abs(x - med)).sort((p, q) => p - q);
      const spread = dev[Math.floor(dev.length * 0.9)] || SPAN_MIN;
      const span = Math.max(spread * SPAN_MUL, SPAN_MIN);
      xmin = med - span;
      xmax = med + span;
    } else {
      xmin = -10; xmax = 10;
    }
  } else {
    const scale = Math.max(1, Math.cbrt(Math.abs(a6) + 1),
                              Math.sqrt(Math.abs(a4) + 1),
                              Math.abs(a2) + 1);
    xmin = -3 * scale; xmax = 3 * scale;
  }

  const yEst = Math.max(
    yMagAt(xmin, a1, a2, a3, a4, a6),
    yMagAt(xmax, a1, a2, a3, a4, a6),
    yMagAt((xmin + xmax) / 2, a1, a2, a3, a4, a6),
    1
  );
  const ymag = Math.max(yEst * Y_MUL, (xmax - xmin) * Y_FALLBACK_MUL);

  baseDomain.x = [xmin, xmax];
  baseDomain.y = [-ymag, ymag];

  const absX = xs.map(Math.abs).filter(v => v > 0).sort((p, q) => p - q);
  const absY = ys.map(Math.abs).filter(v => v > 0).sort((p, q) => p - q);
  symlogConst.x = absX.length ? Math.max(absX[Math.floor(absX.length / 2)], 1) : 1;
  symlogConst.y = absY.length ? Math.max(absY[Math.floor(absY.length / 2)], 1) : 1;
}

function yMagAt(x, a1, a2, a3, a4, a6) {
  const rhs = ((x + a2) * x + a4) * x + a6;
  const lin = a1 * x + a3;
  const D = lin * lin + 4 * rhs;
  if (!Number.isFinite(D) || D < 0) return 0;
  const r = Math.sqrt(D);
  return Math.max(Math.abs((-lin + r) / 2), Math.abs((-lin - r) / 2));
}

function rebuildScales() {
  const rect = document.getElementById('plot').getBoundingClientRect();
  const W = Math.max(rect.width, 400);
  const H = Math.max(rect.height, 300);
  const curve = state.get('currentCurve');
  if (!curve) return;
  const mode = state.get('scale');

  const pts = curve.points || [];
  const allX = pts.map(p => parseRational(p.xy[0])).filter(Number.isFinite);
  const allY = pts.map(p => parseRational(p.xy[1])).filter(Number.isFinite);

  let xDom, yDom;
  if (mode.x === 'symlog' && allX.length > 0) {
    const mx = Math.max(...allX.map(Math.abs),
                        Math.abs(baseDomain.x[0]), Math.abs(baseDomain.x[1]));
    xDom = [-mx * 1.2, mx * 1.2];
  } else {
    xDom = baseDomain.x.slice();
  }
  if (mode.y === 'symlog' && allY.length > 0) {
    const my = Math.max(...allY.map(Math.abs),
                        Math.abs(baseDomain.y[0]), Math.abs(baseDomain.y[1]));
    yDom = [-my * 1.2, my * 1.2];
  } else {
    yDom = baseDomain.y.slice();
  }

  xScale = buildScale(mode.x, xDom, [45, W - 15], symlogConst.x);
  yScale = buildScale(mode.y, yDom, [H - 30, 15], symlogConst.y);
}

function buildScale(mode, domain, range, c) {
  if (mode === 'symlog') {
    return d3.scaleSymlog().domain(domain).range(range).constant(c);
  }
  return d3.scaleLinear().domain(domain).range(range);
}

export function renderPlot() {
  const curve = state.get('currentCurve');
  if (!curve || !curve.a_invariants) {
    gAxes.selectAll('*').remove();
    gBands.selectAll('*').remove();
    gCurve.selectAll('*').remove();
    gPoints.selectAll('*').remove();
    return;
  }

  if (!xScale || !yScale || state.get('_needRebuild')) {
    computeBaseDomain(curve);
    rebuildScales();
  }

  const rect = document.getElementById('plot').getBoundingClientRect();
  const W = Math.max(rect.width, 400);
  const H = Math.max(rect.height, 300);
  svg.attr('width', W).attr('height', H);

  const zoom = state.get('zoom') || d3.zoomIdentity;
  const xs = zoom.rescaleX(xScale);
  const ys = zoom.rescaleY(yScale);

  // --- 軸 ---
  gAxes.selectAll('*').remove();
  const y0 = ys(0), x0 = xs(0);
  if (Number.isFinite(y0)) {
    gAxes.append('line')
      .attr('x1', 0).attr('x2', W)
      .attr('y1', y0).attr('y2', y0)
      .attr('stroke', '#bbb');
  }
  if (Number.isFinite(x0)) {
    gAxes.append('line')
      .attr('x1', x0).attr('x2', x0)
      .attr('y1', 0).attr('y2', H)
      .attr('stroke', '#bbb');
  }
  gAxes.append('g')
    .attr('transform', `translate(0,${H - 25})`)
    .call(d3.axisBottom(xs).ticks(8));
  gAxes.append('g')
    .attr('transform', 'translate(40,0)')
    .call(d3.axisLeft(ys).ticks(8));

  // --- y1 帯 ---
  drawBands(curve, xs, W, H);

  // --- 曲線 ---
  const a = curve.a_invariants.map(Number);
  const [a1, a2, a3, a4, a6] = a;
  const [xL, xR] = xs.domain();
  const n = Math.max(X_SAMPLES_MIN, Math.floor(W * 2.5));
  const upper = [], lower = [];
  const pxL = xs(xL), pxR = xs(xR);
  for (let i = 0; i <= n; i++) {
    const px = pxL + (pxR - pxL) * i / n;
    const x = xs.invert(px);
    const rhs = ((x + a2) * x + a4) * x + a6;
    const lin = a1 * x + a3;
    const D = lin * lin + 4 * rhs;
    if (Number.isFinite(D) && D >= 0) {
      const r = Math.sqrt(D);
      upper.push([x, (-lin + r) / 2]);
      lower.push([x, (-lin - r) / 2]);
    } else {
      upper.push(null);
      lower.push(null);
    }
  }
  const line = d3.line()
    .defined(d => d !== null && Number.isFinite(xs(d[0])) && Number.isFinite(ys(d[1])))
    .x(d => xs(d[0]))
    .y(d => ys(d[1]));

  gCurve.selectAll('path').remove();
  gCurve.append('path').attr('d', line(upper))
    .attr('fill', 'none').attr('stroke', '#222').attr('stroke-width', 1.5);
  gCurve.append('path').attr('d', line(lower))
    .attr('fill', 'none').attr('stroke', '#222').attr('stroke-width', 1.5);

  // --- 点 ---
  const points = curve.points || [];
  const sel = gPoints.selectAll('g.ptg').data(points);
  sel.exit().remove();
  const enter = sel.enter().append('g').attr('class', 'ptg');
  enter.append('path');
  const all = enter.merge(sel);

  all.attr('transform', (d) => {
    const x = parseRational(d.xy[0]);
    const y = parseRational(d.xy[1]);
    const px = xs(x), py = ys(y);
    if (!Number.isFinite(px) || !Number.isFinite(py)) {
      return 'translate(-9999,-9999)';
    }
    return `translate(${px},${py})`;
  });

  all.select('path')
    .attr('d', (d) => markerPath(d))
    .attr('fill', (d) => COLORS[d.role] || '#888')
    .attr('stroke', '#000')
    .attr('stroke-width', (d) => d.all_same_sign ? 3 : 0.5)
    .attr('opacity', (d) => d.role === 'combination' ? 0.7 : 0.95)
    .style('cursor', 'pointer');

  all.on('click', (event, d) => {
    event.stopPropagation();
    state.set({ selectedPoint: d });
  });
}

// y1 軸 (= プロットの x 軸) に沿った縦帯を描く
function drawBands(curve, xs, W, H) {
  gBands.selectAll('*').remove();
  if (!state.get('showBands')) return;
  const bands = curve.y1_bands || [];
  if (bands.length === 0) return;

  const [xDomL, xDomR] = xs.domain();

  for (const band of bands) {
    if (band.y1_min === null || band.y1_max === null) continue;
    const j = band.branch | 0;
    const color = BAND_COLORS[j % BAND_COLORS.length];

    // ドメインへクリップ
    let lo = Math.max(band.y1_min, xDomL);
    let hi = Math.min(band.y1_max, xDomR);
    if (!(hi > lo)) continue;

    let pxL = xs(lo);
    let pxR = xs(hi);
    if (!Number.isFinite(pxL) || !Number.isFinite(pxR)) continue;
    if (pxR < pxL) [pxL, pxR] = [pxR, pxL];

    const g = gBands.append('g').attr('class', `band b${j}`);
    g.append('rect')
      .attr('x', pxL)
      .attr('y', 0)
      .attr('width', Math.max(pxR - pxL, 0.5))
      .attr('height', H)
      .attr('fill', color)
      .attr('fill-opacity', BAND_OPACITY)
      .attr('stroke', color)
      .attr('stroke-opacity', 0.4)
      .attr('stroke-width', 0.5)
      .append('title')
      .text(`branch j=${j}: y1 ∈ [${band.y1_min.toExponential(4)}, ${band.y1_max.toExponential(4)}]`);
  }
}

function markerPath(d) {
  const sym = SVG_SYMBOL[d.role] || SVG_SYMBOL.combination;
  return d3.symbol().type(sym.type).size(sym.size)();
}

subscribe('scale', () => {
  if (state.get('currentCurve')) {
    rebuildScales();
  }
});

subscribe('currentCurve', (curve) => {
  if (!curve) return;
  computeBaseDomain(curve);
  rebuildScales();
  svg.call(zoomBehavior.transform, d3.zoomIdentity);
  state.set({ zoom: d3.zoomIdentity });
});

subscribe('showBands', () => {
  renderPlot();
});