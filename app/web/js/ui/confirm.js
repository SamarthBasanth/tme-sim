// confirm.js — one arm-then-confirm affordance for every destructive control.
//
// Three buttons in this Studio deleted something on a single click with no way
// back: `remove` on one of the operator's own compiled compounds (which holds
// unpublished SMILES that exist nowhere else), `×` on a project, and `×` on a
// specimen. A stray click, a mis-aimed pointer on a 24px control, or a Tab +
// Space with the wrong row focused, and the thing was gone.
//
// A native window.confirm() would do the job and is what the protocol editor
// uses for its buffer-replace, but it is a modal that steals focus for a
// two-word question and gives no room to say WHAT is about to be lost. This is
// the in-place form: the button itself becomes the question.
//
//   click 1   the control turns danger-coloured and states the verb and the
//             target ("delete ACME-101?"). Nothing has happened yet.
//   click 2   the action runs.
//   Escape, blur, a click elsewhere, or 6 s of silence -> disarmed, untouched.
//
// Focus never moves, so the keyboard path is Space-Space, and the armed label
// is the accessible name, so a screen reader hears the question before the
// second press. `aria-live` is deliberately NOT used: the armed state is
// announced by the focused control's own name changing, which is what a
// screen-reader user is already listening to.

const DISARM_MS = 6000;

/**
 * Make an existing button two-step. Idempotent per element.
 *
 * @param {HTMLButtonElement} btn
 * @param {object} o
 * @param {string} o.confirmLabel  what the button says once armed — phrase it
 *   as the question, e.g. `delete ACME-101?`
 * @param {() => any} o.onConfirm  run on the second press
 * @param {string} [o.armedTitle]  tooltip while armed; defaults to a sentence
 *   saying how to back out
 */
export function confirmable(btn, { confirmLabel, onConfirm, armedTitle }) {
  if (btn.dataset.confirmable === '1') return btn;
  btn.dataset.confirmable = '1';

  const restLabel = btn.textContent;
  const restTitle = btn.getAttribute('title') || '';
  const restAria = btn.getAttribute('aria-label') || '';
  const restClass = btn.className;
  let armed = false;
  let timer = 0;

  function disarm() {
    if (!armed) return;
    armed = false;
    clearTimeout(timer);
    btn.textContent = restLabel;
    btn.className = restClass;
    btn.setAttribute('title', restTitle);
    if (restAria) btn.setAttribute('aria-label', restAria);
    else btn.removeAttribute('aria-label');
    document.removeEventListener('keydown', onKey, true);
    document.removeEventListener('pointerdown', onAway, true);
    btn.removeEventListener('blur', disarm);
  }

  function onKey(e) {
    if (e.key === 'Escape') { e.stopPropagation(); disarm(); }
  }
  function onAway(e) { if (e.target !== btn) disarm(); }

  function arm() {
    armed = true;
    btn.textContent = confirmLabel;
    btn.className = restClass + ' armed';
    btn.setAttribute('title', armedTitle
      || 'Press again to confirm. Escape, or a click anywhere else, cancels.');
    // The armed text IS the question, so it must also be the accessible name;
    // a stale aria-label would leave a screen reader hearing "remove" on the
    // press that actually deletes.
    btn.setAttribute('aria-label', confirmLabel);
    document.addEventListener('keydown', onKey, true);
    document.addEventListener('pointerdown', onAway, true);
    btn.addEventListener('blur', disarm);
    timer = setTimeout(disarm, DISARM_MS);
  }

  btn.addEventListener('click', (e) => {
    e.stopPropagation();
    e.preventDefault();
    if (!armed) { arm(); return; }
    disarm();
    onConfirm();
  });

  return btn;
}
