// File: app.js
// 楕円曲線ビューア v2 — エントリポイント

import { state, subscribe } from './js/state.js';
import { loadIndex, loadCurve } from './js/loader.js';
import { renderHeatmap } from './js/heatmap.js';
import { initPlot, renderPlot } from './js/plot.js';
import { renderCurveInfo, renderPointInfo } from './js/panel.js';
import { initControls } from './js/controls.js';

window.addEventListener('error', (e) => {
  console.error(e);
  showFatal(`JS Error: ${e.message} @ ${e.filename}:${e.lineno}`);
});
window.addEventListener('unhandledrejection', (e) => {
  console.error(e);
  showFatal(`Promise rejected: ${e.reason}`);
});

function showFatal(msg) {
  const el = document.getElementById('curve-info');
  if (el) el.innerHTML = `<span style="color:red">${msg}</span>`;
}

async function bootstrap() {
  let index;
  try {
    index = await loadIndex();
  } catch (e) {
    showFatal(`index.json 読み込み失敗: ${e}`);
    return;
  }
  state.set({ index });

  renderHeatmap();
  initPlot();
  initControls();

  // セル選択時に曲線を読み込む
  subscribe('selection', async (sel) => {
    if (!sel) return;
    try {
      const curve = await loadCurve(sel.k, sel.n);
      state.set({ currentCurve: curve, selectedPoint: null });
    } catch (e) {
      showFatal(`曲線データ読み込み失敗: ${e}`);
    }
  });

  subscribe('currentCurve', () => {
    renderCurveInfo();
    renderPlot();
  });
  subscribe('selectedPoint', () => {
    renderPointInfo();
  });
  subscribe('scale', () => renderPlot());
  subscribe('zoom', () => renderPlot());
}

bootstrap();