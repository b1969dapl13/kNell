// File: js/controls.js
import { state } from './state.js';
import { resetView } from './plot.js';

export function initControls() {
  document.getElementById('reset-view').addEventListener('click', resetView);

  document.querySelectorAll('.scale-ctrl button[data-axis]').forEach((btn) => {
    btn.addEventListener('click', () => {
      const axis = btn.dataset.axis;
      const mode = btn.dataset.mode;
      document.querySelectorAll(`.scale-ctrl button[data-axis="${axis}"]`)
        .forEach((b) => b.classList.toggle('active', b === btn));
      const current = state.get('scale');
      state.set({ scale: { ...current, [axis]: mode } });
    });
  });
}