// File: js/panel.js
import { state } from './state.js';

export function renderCurveInfo() {
  const curve = state.get('currentCurve');
  const el = document.getElementById('curve-info');
  if (!curve) {
    el.innerHTML = '<i>(k, N) を選択してください</i>';
    return;
  }
  if (curve.status !== 'ok') {
    el.innerHTML = `<b>(k=${curve.k}, N=${curve.n})</b> status=${curve.status}`
      + (curve.note ? `<br><span class="info-label">note:</span> ${escape(curve.note)}` : '');
    return;
  }

  const torsionStr = formatTorsionStructure(curve.torsion_structure);
  const freeGenLines = (curve.free_generators || []).map(
    (g) => `  ${g.label} = (${g.xy[0]}, ${g.xy[1]})`
  ).join('\n');
  const torsionGenLines = (curve.torsion_generators || []).map(
    (g) => `  ${g.label} (order ${g.order}) = (${g.xy[0]}, ${g.xy[1]})`
  ).join('\n');

  const rankStr = curve.rank !== null && curve.rank !== undefined
    ? String(curve.rank)
    : (curve.rank_bounds ? `[${curve.rank_bounds[0]}, ${curve.rank_bounds[1]}]` : '?');

  el.textContent =
`(k=${curve.k}, N=${curve.n})  status=ok
a-invariants: [${curve.a_invariants.join(', ')}]
j-invariant : ${curve.j_invariant}
discriminant: ${curve.discriminant}
conductor   : ${curve.conductor}
rank        : ${rankStr}
torsion     : ${torsionStr} (order ${curve.torsion_order})
alpha       : ${curve.alpha}
beta        : ${curve.beta}
M (sweep)   : ${curve.M}
points      : ${curve.points_count}

free generators (P_i):
${freeGenLines || '  (none)'}

torsion generators (T_j):
${torsionGenLines || '  (none)'}

note: ${curve.proof_note || ''}`;
}

export function renderPointInfo() {
  const pt = state.get('selectedPoint');
  const el = document.getElementById('pt-info');
  if (!pt) {
    el.innerHTML = '<i>点をクリックすると詳細を表示</i>';
    return;
  }
  const abcStr = pt.abc ? `(${pt.abc.join(' : ')})` : '(N/A)';
  const signStr = pt.sign_pattern
    ? pt.sign_pattern.map(s => s > 0 ? '+' : s < 0 ? '-' : '0').join('')
    : '?';
  el.textContent =
`label : ${pt.label}
role  : ${pt.role}
xy    : (${pt.xy[0]}, ${pt.xy[1]})
abc   : ${abcStr}
sign  : ${signStr}${pt.all_same_sign ? '  <-- all same sign' : ''}
free  : [${(pt.free_coeffs || []).join(', ')}]
torsion: [${(pt.torsion_coeffs || []).join(', ') || '-'}]${pt.order ? `\norder : ${pt.order}` : ''}`;
}

function formatTorsionStructure(arr) {
  if (!arr || arr.length === 0) return 'trivial';
  return arr.map(o => `Z/${o}`).join(' + ');
}

function escape(s) {
  return String(s).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
}