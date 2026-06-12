// File: js/rational.js
export function parseRational(s) {
  if (typeof s === 'number') return s;
  if (typeof s !== 'string') return NaN;
  const t = s.trim();
  if (t.includes('/')) {
    const [num, den] = t.split('/');
    const n = parseFloat(num);
    const d = parseFloat(den);
    if (!isFinite(n) || !isFinite(d) || d === 0) return NaN;
    return n / d;
  }
  return parseFloat(t);
}