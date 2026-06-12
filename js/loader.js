// File: js/loader.js
const _curveCache = new Map();

function signed5(v) {
  return (v >= 0 ? '+' : '-') + String(Math.abs(v)).padStart(4, '0');
}

export async function loadIndex() {
  const res = await fetch('data/index.json');
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  return await res.json();
}

export async function loadCurve(k, n) {
  const key = `${k},${n}`;
  if (_curveCache.has(key)) return _curveCache.get(key);
  const fname = `curve_k${signed5(k)}_N${signed5(n)}.json`;
  const res = await fetch(`data/${fname}`);
  if (!res.ok) throw new Error(`HTTP ${res.status} for ${fname}`);
  const data = await res.json();
  _curveCache.set(key, data);
  return data;
}