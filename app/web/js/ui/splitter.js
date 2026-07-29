// splitter.js — drag-resizable panes. Keyboard-operable (arrow keys), which
// matters because a pane you cannot resize without a mouse is a pane a
// keyboard user cannot read.

export function makeSplitter(el, {
  axis = 'x', varName, min = 160, max = 720, invert = false, onChange, ariaLabel,
}) {
  const root = document.documentElement;
  const read = () => parseFloat(getComputedStyle(root).getPropertyValue(varName)) || min;
  const write = (v) => {
    const clamped = Math.max(min, Math.min(max, v));
    root.style.setProperty(varName, clamped + 'px');
    if (onChange) onChange(clamped);
    return clamped;
  };

  el.setAttribute('role', 'separator');
  el.setAttribute('tabindex', '0');
  el.setAttribute('aria-orientation', axis === 'x' ? 'vertical' : 'horizontal');
  el.setAttribute('aria-valuemin', String(min));
  el.setAttribute('aria-valuemax', String(max));
  // role="separator" carries no accessible name on its own — without one a
  // screen reader announces "separator, 214" and nothing about what it
  // resizes or that arrow keys act on it at all.
  if (ariaLabel) el.setAttribute('aria-label', ariaLabel);

  const sync = () => el.setAttribute('aria-valuenow', String(Math.round(read())));
  sync();

  let start = 0, base = 0, dragging = false;

  const move = (e) => {
    if (!dragging) return;
    const p = axis === 'x' ? e.clientX : e.clientY;
    const d = (p - start) * (invert ? -1 : 1);
    write(base + d);
    sync();
    e.preventDefault();
  };
  const up = () => {
    dragging = false;
    el.classList.remove('dragging');
    document.body.style.cursor = '';
    document.body.style.userSelect = '';
    window.removeEventListener('pointermove', move);
    window.removeEventListener('pointerup', up);
  };

  el.addEventListener('pointerdown', (e) => {
    dragging = true;
    start = axis === 'x' ? e.clientX : e.clientY;
    base = read();
    el.classList.add('dragging');
    document.body.style.cursor = axis === 'x' ? 'col-resize' : 'row-resize';
    document.body.style.userSelect = 'none';
    window.addEventListener('pointermove', move);
    window.addEventListener('pointerup', up, { once: true });
    e.preventDefault();
  });

  el.addEventListener('keydown', (e) => {
    const dec = axis === 'x' ? 'ArrowLeft' : 'ArrowUp';
    const inc = axis === 'x' ? 'ArrowRight' : 'ArrowDown';
    const step = e.shiftKey ? 40 : 8;
    if (e.key === dec) { write(read() - step * (invert ? -1 : 1)); sync(); e.preventDefault(); }
    else if (e.key === inc) { write(read() + step * (invert ? -1 : 1)); sync(); e.preventDefault(); }
    else if (e.key === 'Home') { write(min); sync(); e.preventDefault(); }
    else if (e.key === 'End') { write(max); sync(); e.preventDefault(); }
  });

  return { read, write: (v) => { write(v); sync(); } };
}
