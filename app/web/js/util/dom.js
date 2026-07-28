// dom.js — a ~60-line hyperscript. No framework, no build step.

/**
 * h('div.cls#id', {attrs}, ...children) -> HTMLElement
 * Attributes: strings/numbers set as attributes; `on*` set as listeners;
 * `style` accepts an object; `dataset` accepts an object; booleans toggle.
 */
export function h(sel, props, ...children) {
  const m = /^([a-zA-Z0-9-]+)?(#[^.\s]+)?((?:\.[^.#\s]+)*)$/.exec(sel);
  if (!m) throw new Error('bad selector: ' + sel);
  const el = document.createElement(m[1] || 'div');
  if (m[2]) el.id = m[2].slice(1);
  if (m[3]) el.className = m[3].slice(1).split('.').join(' ');

  if (props && (props.nodeType || Array.isArray(props) || typeof props !== 'object')) {
    children.unshift(props);
    props = null;
  }
  if (props) {
    for (const k in props) {
      const v = props[k];
      if (v === null || v === undefined || v === false) continue;
      if (k.startsWith('on') && typeof v === 'function') {
        el.addEventListener(k.slice(2).toLowerCase(), v);
      } else if (k === 'style' && typeof v === 'object') {
        Object.assign(el.style, v);
      } else if (k === 'dataset' && typeof v === 'object') {
        Object.assign(el.dataset, v);
      } else if (k === 'class' || k === 'className') {
        el.className = (el.className ? el.className + ' ' : '') + v;
      } else if (k === 'html') {
        el.innerHTML = v;
      } else if (v === true) {
        el.setAttribute(k, '');
      } else {
        el.setAttribute(k, v);
      }
    }
  }
  append(el, children);
  return el;
}

export function append(el, children) {
  for (const c of children) {
    if (c === null || c === undefined || c === false || c === '') continue;
    if (Array.isArray(c)) append(el, c);
    else el.appendChild(c.nodeType ? c : document.createTextNode(String(c)));
  }
  return el;
}

export const $ = (sel, root = document) => root.querySelector(sel);
export const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));

export function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); return el; }
export function fill(el, ...children) { return append(clear(el), children); }

/** Toggle a class and return the element. */
export function cls(el, name, on) { el.classList.toggle(name, !!on); return el; }

/** Debounce with a leading-edge option. */
export function debounce(fn, ms) {
  let t = 0;
  const d = (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); };
  d.cancel = () => clearTimeout(t);
  d.flush = (...a) => { clearTimeout(t); fn(...a); };
  return d;
}

/** rAF-coalesced callback. */
export function raf(fn) {
  let pending = 0;
  return (...a) => {
    if (pending) return;
    pending = requestAnimationFrame(() => { pending = 0; fn(...a); });
  };
}
