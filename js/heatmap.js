// File: js/heatmap.js
import { state } from './state.js';

const CELL_SIZE = 16;
const KAXIS_HEIGHT = 14;
const NAXIS_WIDTH  = 22;

export function renderHeatmap() {
  const index = state.get('index');
  if (!index) return;

  const [kmin, kmax] = index.k_range;
  const [nmin, nmax] = index.n_range;
  const container = document.getElementById('grid-container');
  container.innerHTML = '';

  const table = el('div', 'grid-table');
  container.appendChild(table);

  // 左上隅
  table.appendChild(el('div', 'grid-corner'));

  // k 軸 (上部)
  const kHeader = el('div', 'grid-header-k');
  kHeader.style.gridTemplateColumns = `repeat(${kmax - kmin + 1}, ${CELL_SIZE}px)`;
  for (let k = kmin; k <= kmax; k++) {
    const a = el('div', 'axisnum col' + (k === 0 ? ' zero-k' : ''));
    a.textContent = k;
    kHeader.appendChild(a);
  }
  table.appendChild(kHeader);

  // N 軸 (左)
  const nHeader = el('div', 'grid-header-n');
  nHeader.style.gridTemplateRows = `repeat(${nmax - nmin + 1}, ${CELL_SIZE}px)`;
  for (let n = nmax; n >= nmin; n--) {
    const a = el('div', 'axisnum row' + (n === 0 ? ' zero-n' : ''));
    a.textContent = n;
    nHeader.appendChild(a);
  }
  table.appendChild(nHeader);

  // 本体
  const wrapper = el('div', 'grid-body-wrapper');
  const body = el('div', 'grid');
  body.style.gridTemplateColumns = `repeat(${kmax - kmin + 1}, ${CELL_SIZE}px)`;
  wrapper.appendChild(body);
  table.appendChild(wrapper);

  const hoverLabel = document.getElementById('hover-label');
  let selectedDiv = null;

  for (let n = nmax; n >= nmin; n--) {
    for (let k = kmin; k <= kmax; k++) {
      const key = `${k},${n}`;
      const cellInfo = index.cells[key];  // {status, rank, has_same_sign_point} or undefined
      const status = cellInfo ? cellInfo.status : null;

      const cell = el('div', 'cell');
      if (status) cell.classList.add(status);
      if (k === 0) cell.classList.add('zero-axis-col');
      if (n === 0) cell.classList.add('zero-axis-row');
      if (cellInfo && cellInfo.has_same_sign_point) {
        cell.classList.add('has-same-sign');
      }

      cell.addEventListener('mouseenter', () => {
        const tag = status || 'DB に無し';
        hoverLabel.textContent = `k=${k}, N=${n} [${tag}]`;
      });

      if (status === 'ok') {
        cell.addEventListener('click', () => {
          if (selectedDiv) selectedDiv.classList.remove('selected');
          cell.classList.add('selected');
          selectedDiv = cell;
          state.set({ selection: { k, n } });
        });
      }
      body.appendChild(cell);
    }
  }
}

function el(tag, className) {
  const e = document.createElement(tag);
  if (className) e.className = className;
  return e;
}