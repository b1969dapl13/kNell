// File: js/state.js
const _state = {
  index: null,
  currentCurve: null,
  selection: null,
  selectedPoint: null,
  scale: { x: 'linear', y: 'linear' },
  zoom: null,           // 初期化は initPlot 内で d3.zoomIdentity を入れる
  hover: null,
  showBands: true,      // y1 帯の表示 ON/OFF
};

const _subs = new Map();

export const state = {
  get(key) { return _state[key]; },
  getAll() { return { ..._state }; },
  set(patch) {
    const changed = [];
    for (const k of Object.keys(patch)) {
      if (_state[k] !== patch[k]) {
        _state[k] = patch[k];
        changed.push(k);
      }
    }
    for (const k of changed) {
      const handlers = _subs.get(k) || [];
      handlers.forEach((h) => {
        try { h(_state[k]); } catch (e) { console.error(e); }
      });
    }
  },
};

export function subscribe(key, handler) {
  if (!_subs.has(key)) _subs.set(key, []);
  _subs.get(key).push(handler);
}