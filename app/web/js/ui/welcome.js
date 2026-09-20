// welcome.js — the first sixty seconds.
//
// WHAT THIS FILE IS FOR
// ---------------------
// Someone arrives knowing nothing. Before they touch anything they need three
// answers, in this order, and the FIRST STEP has to carry all three:
//
//   1. what does this do
//   2. what does it NOT do
//   3. what do I press
//
// Everything after step one is optional depth. That ordering is deliberate:
// until 2026-09-20 this tour opened on a paragraph of self-description and put
// the concrete first action on step SEVEN, so a reader who bailed after the
// first card left with a claim and no instruction — which is the one shape of
// onboarding that is worse than none.
//
// WHY THE LIMITS ARE IN THE FIRST CARD RATHER THAN THE LAST
// --------------------------------------------------------
// Onboarding copy is where overclaiming happens, because explaining is the act
// of compressing. A limit that arrives on card six arrives after the reader has
// already formed a model of the product from card one. So the first card states
// what is not defensible in the same breath, at the same type size, as what is —
// and then states the reading that DOES survive, with the reason it survives,
// because a product that lists only what it cannot do has told the reader
// nothing they can use. The limits are quoted from the project's own standing
// list (`ui/provenance.js::STANDING_LIMITS`), not re-derived here.
//
// WHY IT IS A POPOVER AND NOT A MODAL
// -----------------------------------
// It anchors to the real element each step is about (the SECTION dropdown, the
// Submit run button, the toolbar) and leaves the rest of the Studio live
// underneath, so "look at the dropdown" has a real dropdown to look at. Distinct
// from apikey.js's gate, which legitimately DOES block everything (nothing works
// without a key); this is advisory only and never traps focus.
//
// THREE PROPERTIES IT IS BUILT TO
// -------------------------------
//   SKIPPABLE   — Esc, the ×, and a Skip button on every step. Skipping is not
//                 punished: every claim in here is also reachable in place —
//                 the stage's empty-state notice (ui/overlays.js) repeats the
//                 three answers above, each tool panel states what its tool
//                 tells you (ui/tools.js), and the `tour` button replays this.
//   RESUMABLE   — the step index is persisted, so leaving at card four and
//                 coming back (or reloading) resumes at card four rather than
//                 restarting. `restart` is offered whenever it resumes.
//   NON-TRAPPING— focus is moved onto the card ONCE, at open, and RESTORED to
//                 whatever had it when the tour closes. Tab leaves the card
//                 freely in both directions; nothing here is aria-modal.
//
// Shown automatically once per browser (localStorage flag) and replayable any
// time from the "tour" button in the top bar (main.js's wiring, #btn-tour).

import { h, raf } from '../util/dom.js';

const SEEN_KEY = 'tmesim_welcome_seen_v1';
//: Separate from SEEN_KEY on purpose. "have you been offered this" and "where
//: did you stop" are different facts: a user who skipped on card one should not
//: be nagged again (SEEN), but should still land on card one when they choose
//: to open it (STEP). Bundling them into one flag is how a tour ends up either
//: re-nagging or un-resumable.
const STEP_KEY = 'tmesim_welcome_step_v1';
const MARGIN = 14;
//: The first card carries four labelled blocks and is the only one that has to
//: be read as prose rather than glanced at. 340px is a column of ~6 words.
const WIDE_PX = 560;

export function shouldShowWelcome() {
  try { return !localStorage.getItem(SEEN_KEY); }
  catch { return false; } // no localStorage (private mode) -- don't nag every load
}

function markSeen() {
  try { localStorage.setItem(SEEN_KEY, '1'); } catch { /* fine, just re-shows next time */ }
}

function readStep(max) {
  try {
    const raw = localStorage.getItem(STEP_KEY);
    if (raw === null) return 0;
    const n = parseInt(raw, 10);
    if (!Number.isFinite(n)) return 0;
    return Math.min(Math.max(n, 0), max);
  } catch { return 0; }
}
function writeStep(n) { try { localStorage.setItem(STEP_KEY, String(n)); } catch { /* not resumable, still usable */ } }
function clearStep() { try { localStorage.removeItem(STEP_KEY); } catch { /* ignore */ } }

/**
 * Did this session arrive from the workspace home at `/`?
 *
 * WORKSPACE_CONTRACT.md section 8: the Studio reads `?project=<id>` when it is
 * opened from a project card, and `main.js::showProjectContext` is what puts
 * the `← workspace` link and the project name in the top bar when it does. The
 * presence of that parameter is also the only evidence available in the browser
 * that the landing page EXISTS on this deployment — it is the home page that
 * writes it.
 *
 * So the project-aware copy is gated on it rather than on a hard-coded belief
 * about what `/` serves: on a deployment where `/` is still the Studio the
 * parameter is never present, and the tour never names a link or a page that is
 * not on the screen. A label must not describe behaviour the code does not have
 * (UI_ARCHITECTURE.md section 5), and that rule does not stop applying because
 * the behaviour belongs to another module.
 */
function openedFromWorkspace() {
  try { return !!new URLSearchParams(window.location.search).get('project'); }
  catch { return false; }
}

/** Each step: a title, body content, and an optional CSS selector for the
 * real element this step is anchored to (`.tour-highlight` outline + the
 * popover positions itself next to it). No selector -> centered, for the
 * cards that are not about any one control. Selectors are best-effort:
 * if not found (different viewport, a panel not built), the step still shows,
 * centered, rather than anchoring to nothing. `wide` widens the card. */
const STEPS = [
  {
    title: 'What this is',
    selector: null,
    wide: true,
    body: () => [
      h('div.hr-lbl', 'What it does'),
      h('p',
        'Vivum reconstructs a real tumour tissue section from a measured ' +
        'spatial-transcriptomics assay — every cell on the stage is a cell ' +
        'somebody sequenced, at its own coordinates — and then solves drug ' +
        'transport through that geometry. The output is a ', h('b', 'map'),
        ': for one compound on one section, where the drug reaches and where ' +
        'it does not. The viewer, the mechanism breakdown and every metric are ' +
        'computed from that solve. Nothing on screen is an illustration.'),

      // The tour was once the one first-run surface in the Studio that stated
      // no limit at all — it introduced a computed tumour response to someone
      // who had never seen the product and handed them straight to Submit run.
      // What follows is the project's own standing limit
      // (ui/provenance.js::STANDING_LIMITS), and it belongs in the first thing
      // a new user reads, not only in the panel they reach afterwards.
      h('div.hr-lbl', 'What it does not do'),
      h('p',
        'It does not predict what a patient will do. These sections carry ',
        h('b', 'no treatment-response label'),
        ', so no number in this Studio has ever been scored against a patient ' +
        'outcome. It is not a diagnostic, not a treatment-planning tool, and ',
        h('b', 'NOT FOR CLINICAL USE'),
        '. Two readings in particular are off the table: the absolute level — ' +
        'what fraction of a tumour dies — and setting one compound beside ' +
        'another. Both were measured, both failed, and both were withdrawn.'),

      h('div.hr-lbl', 'What survives, and why that is the interesting part'),
      h('p',
        'Inside ', h('b', 'one'), ' map, for ', h('b', 'one'), ' compound: the ' +
        'exposure in one region relative to the best-perfused tissue in that ' +
        'same map. That reading survives for a reason rather than by ' +
        'concession — the two factors measured to destroy the absolute level ' +
        'rescale every point of the map together, so they cancel out of a ' +
        'ratio between two points. What they do not cancel is error in the ' +
        'transport law itself, which tilts the falloff instead of lifting it. ' +
        'That law was scored blind and unfitted against 71 compounds — measured ' +
        'by one method, in one cell line — and reproduces about a quarter of ' +
        'the spread those measurements show. Read it as the best available ' +
        'estimate of a shape, not as a measurement of one. Both halves of this ' +
        'are restated, per run, in the CAVEATS panel — and never in smaller ' +
        'type than the number they qualify.'),

      h('div.hr-lbl', 'What to press'),
      h('p',
        'Left rail: a section is already chosen and a compound is already ' +
        'filled in. Press ', h('b', 'Submit run'),
        '. Right rail, when it lands: read ', h('b', 'The finding'), ', then ',
        h('b', 'Caveats'), '. That is the entire loop, and everything else in ' +
        'here is easier to read once a real result is on the screen.'),
      // Only when the workspace handed this session a project. A reader who
      // arrived from a project card has one question the others do not have —
      // "does my run come back with me" — and it is worth answering on the
      // first card rather than the last.
      openedFromWorkspace()
        ? h('p',
          'You came in from a project, and the top bar names it. A run is not ' +
          'filed under it automatically: when one finishes you are offered a ' +
          'specimen to attach it to, and nothing is attached until you pick one.')
        : null,
      h('p.dim',
        `The ${STEPS.length - 1} cards after this one each point at the control they ` +
        'are about. Leave whenever you like — this remembers where you stopped, ' +
        'and ', h('b', 'tour'), ' in the top bar brings it back.'),
    ],
  },
  {
    title: '1 · The specimen you are looking at',
    selector: '#f-section',
    side: 'bottom',
    body: () => [
      /*
       * DO NOT PUT A SECTION COUNT BACK IN THIS SENTENCE.
       *
       * It said "38 real sections across 13 patients ... spanning
       * glioblastoma and breast cancer". That was true on the machine it was
       * written on and false on a fresh install, which is the only machine
       * whose reader is being addressed here. Two independent reasons:
       *
       *   (1) `app/runs/` is gitignored, so the ingest registry does not ship.
       *       A clean checkout has SEVEN built-in sections -- all glioblastoma
       *       (Xenium 319-panel, Xenium IO 480-panel, GSE324268, and four
       *       Visium). There is no built-in breast section. The other 31 and
       *       the breast tissue were locally ingested by the developer.
       *   (2) The count is now user-mutable: the "Your tissue" panel changes
       *       it, so any literal is stale the moment a field tester ingests.
       *
       * A first-run tour is the one surface a new user cannot check against
       * anything else, so a wrong number here is caught immediately and costs
       * the trust every other disclosure depends on. Describe the KIND of
       * thing in the list and let the list speak for its own length.
       */
      h('p',
        'This ', h('b', 'SECTION'), ' dropdown is not one fixed demo tumour. ' +
        'It lists the built-in specimens this project ships measured evidence ' +
        'for — real glioblastoma tissue across two assay families, Xenium and ' +
        'Visium — alongside any section registered on this server. Switching ' +
        'sections loads genuinely different measured tissue, not a re-skin of ' +
        'the same one.'),
      h('p',
        'A section\'s own caveats print directly underneath this control, at ' +
        'the moment you select it — not in a help file, and not after the run. ' +
        'What a given panel cannot support is a property of the specimen, so ' +
        'it is stated where the specimen is chosen.'),
    ],
  },
  {
    title: '2 · Or point it at your own tissue',
    selector: '#sect-tissue',
    side: 'right',
    body: () => [
      h('p',
        'The list is not the limit. ', h('b', 'Your tissue'),
        ' registers a 10x Xenium bundle, a Space Ranger (Visium) output or a ' +
        'bare MEX export that already sits in a ', h('b', 'directory on the ' +
        'machine running this server'),
        ': you name the path, that server process reads it off its own disk, ' +
        'and the browser sends no tissue anywhere. It then tells you what the ' +
        'specimen can and cannot support before anything is run — and a ' +
        'directory it refuses names the reason and what to do about it.'),
      h('p',
        'Two things it will say to you, and both are load-bearing: a section ' +
        'you register this way is ', h('b', 'not validated'), ' — the built-in ' +
        'specimens are the ones this project ships measured evidence for — and ' +
        'Vivum may be pointed at ', h('b', 'de-identified tissue only'), '.'),
    ],
  },
  {
    title: '3 · Set it up, then submit',
    selector: '#btn-submit-run',
    side: 'top',
    body: () => [
      h('p',
        'Pick a compound, adjust dose and schedule if you want, then ',
        h('b', 'Submit run'), '. ', h('b', 'Validate'),
        ' beside it checks the RunSpec without spending a solve, and the chip ' +
        'next to them says which of those two states you are in. At the ' +
        'default grid spacing a run is usually well under a minute.'),
      h('p',
        'You can leave and come back. A run in flight is picked up again when ' +
        'you reload, and the run history in the dock along the bottom reopens ' +
        'any finished run exactly as it was.'),
    ],
  },
  {
    title: '4 · Read the result in this order',
    selector: '#sect-finding',
    side: 'left',
    body: () => [
      // The tour previously stopped at "it populates", which leaves a new user
      // staring at six numbers with no idea which one is the answer. Naming the
      // reading order is the actual skill being taught, and it doubles as the
      // product's central claim: the map is the output, the level is not.
      h('p',
        'Start at ', h('b', 'The finding'), ', at the top of the right rail. ' +
        'It states what the run found in counts of cells, with the ' +
        'denominator visible beside every fraction, and it names the ' +
        'measurement that would refute it.'),
      h('p',
        'Then ', h('b', 'Mechanism'), ', below it: why each surviving tumour ' +
        'cell survived. Cells the drug never reached are counted separately ' +
        'from cells it reached and under-dosed, and from cells that resisted ' +
        'it. That separation — ', h('b', 'delivery or dose'),
        ' — is the question this tool exists to answer.'),
      h('p',
        h('b', 'Metrics'), ' sits below both on purpose. It is supporting ' +
        'detail, the absolute values among them are not defensible on their ' +
        'own, and the rail is ordered so that reading it top-down does not ' +
        'hand you one of those first.'),
    ],
  },
  {
    title: '5 · Caveats are not the small print',
    selector: '.caveats',
    side: 'left',
    body: () => [
      h('p',
        h('b', 'This panel'), ' states plainly what each number can and ' +
        'cannot be read as, for the run you have open. It is never truncated, ' +
        'never scrolled inside itself, and never set smaller than the number ' +
        'it qualifies. Read it before quoting anything.'),
      h('p',
        'A stored run keeps the caveats it was written with, word for word, ' +
        'even where this project has since retracted one of them — such a ' +
        'record is marked as superseded rather than quietly rewritten. An ' +
        'honest "we do not know" is treated here as a result, not as a gap to ' +
        'hide.'),
    ],
  },
  {
    title: '6 · Measure it yourself',
    selector: '#tb-tools',
    side: 'bottom',
    body: () => [
      h('p',
        'The map is meant to be interrogated, not looked at. Four pointer ' +
        'tools live in this group above the stage, and each one writes its ' +
        'numbers into the ', h('b', 'Measure'), ' group of the right rail with ' +
        'provenance on every value:'),
      h('p',
        h('b', 'pick'), ' pins one cell and fills ', h('b', 'Selected cell'),
        '. ', h('b', 'profile'), ' drags a line and plots the field along it — ' +
        'the gradient out from a vessel, read as numbers instead of guessed ' +
        'from a colour. ', h('b', 'region'), ' drags a rectangle and reports ' +
        'the composition and exposure inside it; two regions in the same map, ' +
        'for the same compound, is the comparison this product defends. ',
        h('b', 'distance'), ' measures a separation in micrometres.'),
      h('p.dim',
        'Keys ', h('kbd', 'V'), ' ', h('kbd', 'L'), ' ', h('kbd', 'R'), ' ',
        h('kbd', 'M'), ' select them, ', h('kbd', 'Esc'),
        ' drops the newest annotation of the active tool, and ', h('kbd', '?'),
        ' lists every key in the application.'),
    ],
  },
  // A tour that ends on a feature leaves the reader with nothing to do. This
  // card hands over one concrete first action, says where the work is kept, and
  // pre-empts the first real failure — which for a new operator is far more
  // likely than a successful exotic run.
  {
    title: '7 · Where the work is kept, and when it goes wrong',
    selector: null,
    side: 'center',
    wide: true,
    body: (ctx) => [
      h('div.hr-lbl', 'Where the work is kept'),
      openedFromWorkspace()
        ? h('p',
          'The Studio is where a run ', h('b', 'happens'), '; the workspace is ' +
          'where it is ', h('b', 'kept'), '. A run submitted here is not filed ' +
          'under your project automatically — when it finishes you are offered ' +
          'a specimen to attach it to, and nothing is attached until you pick ' +
          'one. ', h('b', '← workspace'), ' in the top bar goes back, and ' +
          'leaving changes nothing that belongs to the project: a project holds ' +
          'ids rather than copies, so deleting one never deletes its runs.')
        : h('p',
          'The Studio is where a run ', h('b', 'happens'), '; a project is ' +
          'where it is ', h('b', 'kept'), '. ', h('b', 'projects'),
          ' in the top bar groups runs by specimen and holds the ' +
          'specimen-to-specimen comparison. A project holds ids rather than ' +
          'copies — so detaching a specimen, or deleting the project itself, ' +
          'never deletes a run.'),

      h('div.hr-lbl', 'When it goes wrong'),
      h('p',
        'Failures here are meant to be readable. A run that fails says what ' +
        'went wrong in plain language, and the stage clears rather than ' +
        'leaving the previous run\'s numbers standing under a failed heading — ' +
        'if there is nothing to draw, you get told why.'),
      h('p',
        'If a failure was unexpected the message carries a short ',
        h('b', 'reference code'), '. That code matches one line in the ' +
        'server\'s own log, so quoting it to whoever runs this server points ' +
        'them straight at the cause. Copy it as-is.'),

      h('div.hr-lbl', 'If you would rather script it'),
      h('p',
        h('b', 'protocol'), ' in the top bar is a real code editor for ' +
        'defining parameter sweeps — working templates to start from, syntax ' +
        'highlighting, inline diagnostics, and a plan you confirm before ' +
        'anything is executed. ',
        ctx && ctx.onGoToProtocol
          ? h('button.btn.ghost.tiny', {
            type: 'button',
            onclick: () => { ctx.close(); ctx.onGoToProtocol(); },
          }, 'open the protocol editor')
          : null),

      h('div.hr-lbl', 'A good first run'),
      h('p',
        'Take the section already selected, leave every setting alone, and ' +
        'press ', h('b', 'Submit run'), '. It is the fastest way to see the ' +
        'whole loop, and everything you have just been shown is easier to read ' +
        'once there is a real result on the screen.'),
    ],
  },
];

//: The one tour on screen, if any. Module-level rather than per-call because
//: the invariant belongs to the tour and a caller cannot be trusted to hold it.
let active = null;

function clearHighlight() {
  document.querySelectorAll('.tour-highlight').forEach((el) => el.classList.remove('tour-highlight'));
}

/**
 * Outline the real element this step is about and make sure it can be SEEN.
 *
 * Scrolling and opening are not decoration. Four of the eight steps point into
 * a rail that scrolls and two of those targets live inside a `<details>` that
 * ships collapsed, so without this the card says "this panel" while drawing a
 * 2px outline on something off-screen or folded shut — which is worse than no
 * anchor at all, because the reader goes looking for it.
 *
 * `<details>` is left OPEN afterwards on purpose: it is exactly the state the
 * reader would have produced by clicking the summary themselves, and silently
 * re-collapsing a panel the tour just told them to read is the tour taking
 * something back.
 */
function highlight(selector) {
  clearHighlight();
  if (!selector) return null;
  const el = document.querySelector(selector);
  if (!el) return null;
  const det = typeof el.closest === 'function' ? el.closest('details') : null;
  if (det && !det.open) det.open = true;
  // Instant, never smooth: this runs on every step change and a smooth scroll
  // would still be moving when positionCard() measures the rect.
  if (typeof el.scrollIntoView === 'function') {
    try { el.scrollIntoView({ block: 'nearest', inline: 'nearest' }); }
    catch { el.scrollIntoView(); }
  }
  el.classList.add('tour-highlight');
  return el;
}

/** Position `card` (already in the DOM, so its real measured size is
 * available) next to `targetEl`'s current bounding rect, preferring `side`
 * but flipping or clamping to whatever actually fits in the viewport --
 * a fixed side would run the card off-screen on a small window or a target
 * near an edge. Centers the card (no anchor) when there is no target. */
function positionCard(card, targetEl, side) {
  const cw = card.offsetWidth;
  const ch = card.offsetHeight;
  const vw = window.innerWidth;
  const vh = window.innerHeight;

  if (!targetEl) {
    card.style.top = `${Math.max(MARGIN, (vh - ch) / 2)}px`;
    card.style.left = `${Math.max(MARGIN, (vw - cw) / 2)}px`;
    return;
  }

  const r = targetEl.getBoundingClientRect();
  const fits = {
    bottom: r.bottom + MARGIN + ch <= vh,
    top: r.top - MARGIN - ch >= 0,
    left: r.left - MARGIN - cw >= 0,
    right: r.right + MARGIN + cw <= vw,
  };
  // Preferred side first, then any other side that actually fits, then fall
  // back to the preferred side clamped into the viewport rather than hidden.
  const order = [side, 'bottom', 'top', 'right', 'left'].filter((s, i, a) => a.indexOf(s) === i);
  const chosen = order.find((s) => fits[s]) || side || 'bottom';

  let top, left;
  if (chosen === 'bottom') { top = r.bottom + MARGIN; left = r.left; }
  else if (chosen === 'top') { top = r.top - MARGIN - ch; left = r.left; }
  else if (chosen === 'left') { top = r.top; left = r.left - MARGIN - cw; }
  else { top = r.top; left = r.right + MARGIN; }

  top = Math.min(Math.max(top, MARGIN), vh - ch - MARGIN);
  left = Math.min(Math.max(left, MARGIN), vw - cw - MARGIN);
  card.style.top = `${top}px`;
  card.style.left = `${left}px`;
}

/**
 * @param {{onGoToProtocol?: () => void}} [opts]
 *   `onGoToProtocol` is invoked by the protocol action on the final card, after
 *   the tour has closed itself — the protocol workspace is a full-window sheet
 *   and a popover floating over it is not a layout anybody designed.
 * @returns {{node: HTMLElement, close: () => void, destroy: () => void}}
 */
export function createWelcomeOverlay({ onGoToProtocol } = {}) {
  /*
   * Tear down any tour that is already on screen before mounting this one.
   *
   * Every call used to append a fresh `.tour-popover` and nothing ever removed
   * the previous card, so pressing the `tour` button twice left two dialogs
   * stacked at identical coordinates, three times left three, and so on --
   * measured 1/2/3/4 popovers over four presses. Only the topmost receives the
   * click, so `Next` advanced one card while the ones underneath stayed on
   * step 1 for ever; dismissing the top one uncovered an identical dialog and
   * the tour read as refusing to close. It also put N elements with
   * role="dialog" and the same aria-labelledby into the accessibility tree.
   *
   * The node-sweep below was the fix, and it was only HALF a fix: removing an
   * element does not remove the `keydown` and `resize` listeners its closure
   * registered on `document` and `window`. Those survived every re-open, so a
   * session that opened the tour four times had four live ArrowRight handlers,
   * three of them stepping a detached card and re-running its highlight against
   * the live one. The instance is therefore DESTROYED, which unbinds, and the
   * node sweep is kept underneath it as a belt for anything that predates this.
   */
  if (active) active.destroy();
  document.querySelectorAll('.tour-popover').forEach((el) => el.remove());
  clearHighlight();

  let self = null;
  let i = readStep(STEPS.length - 1);
  // Whoever had the keyboard when the tour opened gets it back when the tour
  // closes. Without this, closing dumps focus on <body> and a keyboard user
  // restarts from the top of a 445-element tab order -- the same cost the
  // opening focus jump exists to avoid, just paid on the way out.
  const opener = document.activeElement;
  //: A resumed tour has to offer a way out of the resume, or a reader who
  //: wanted the beginning has no route to it short of clearing localStorage.
  const resumed = i > 0;

  function restoreFocus() {
    const el = opener;
    if (!el || el === document.body) return;
    if (!el.isConnected || typeof el.focus !== 'function') return;
    try { el.focus(); } catch { /* element refused focus; leave it alone */ }
  }

  function detach() {
    clearHighlight();
    window.removeEventListener('resize', onReflow);
    document.removeEventListener('scroll', onReflow, true);
    document.removeEventListener('keydown', onKey);
    card.remove();
    if (active === self) active = null;
  }

  /** Close without recording anything: used when a SECOND tour replaces this
   *  one, where marking "seen" and restoring focus both belong to the new
   *  instance rather than to this one. */
  function destroy() { detach(); }

  /** Dismissed part-way. The step is KEPT, so re-opening resumes here. */
  function close() { markSeen(); detach(); restoreFocus(); }

  /** Reached the end. The step is CLEARED, so a replay starts at the top. */
  function finish() { markSeen(); clearStep(); detach(); restoreFocus(); }

  function onKey(ev) {
    if (ev.key === 'Escape') { close(); return; }
    // Arrow keys advance the tour ONLY when the operator is not typing. This
    // listener is on `document`, and the tour is deliberately non-modal (see
    // the module comment) — so the rest of the Studio stays live underneath
    // it, including the RunSpec editor and every form field. Without this
    // guard, moving the caret with ArrowLeft inside the JSON editor stepped
    // the tour backwards, which reads as the tour fighting the user for the
    // keyboard. Escape is exempt: closing the tour is what Escape should do
    // from anywhere, and the fields that care already stopPropagation on it.
    const t = ev.target;
    const typing = t && (t.isContentEditable
      || /^(INPUT|TEXTAREA|SELECT)$/.test(t.tagName || ''));
    if (typing || ev.metaKey || ev.ctrlKey || ev.altKey) return;
    if (ev.key === 'ArrowRight') next();
    else if (ev.key === 'ArrowLeft') back();
  }

  // Scroll as well as resize: four of the eight cards anchor into a rail that
  // scrolls independently of the window, and a card that stays put while its
  // target slides away is pointing at the wrong control. rAF-coalesced because
  // a scroll fires far more often than the display presents frames.
  const onReflow = raf(() => {
    const step = STEPS[i];
    const el = step.selector ? document.querySelector(step.selector) : null;
    positionCard(card, el, step.side);
  });

  function render() {
    const step = STEPS[i];
    titleEl.textContent = step.title;
    bodyEl.innerHTML = '';
    for (const node of step.body({ onGoToProtocol, close: finish })) {
      if (node) bodyEl.appendChild(node);
    }
    // `--fs-*` and the 340px default both live in the stylesheet; this is a
    // per-step override of a WIDTH, not of a colour, and it is set through
    // element.style, which CSP does not govern (UI_ARCHITECTURE.md 4.6).
    card.style.maxWidth = step.wide
      ? `${Math.min(WIDE_PX, Math.max(280, window.innerWidth - 2 * MARGIN))}px`
      : '';

    dotsEl.innerHTML = '';
    STEPS.forEach((_, idx) => {
      dotsEl.appendChild(h(`span.tour-dot${idx === i ? '.active' : ''}`));
    });
    // Colour is never the only channel (UI_ARCHITECTURE.md 4.5). The dots
    // differ from each other by background-color alone, so the position is
    // also written out — which is what a screen reader and a CVD reader both
    // get instead of the dots.
    dotsEl.appendChild(h('span.tour-count.mono.dim', {
      style: { marginLeft: 'var(--sp-4)', fontSize: 'var(--fs-xs)' },
    }, `${i + 1} / ${STEPS.length}`));

    backBtn.disabled = i === 0;
    const last = i === STEPS.length - 1;
    nextBtn.textContent = last ? 'Done' : 'Next';
    nextBtn.title = last
      ? 'Close the tour. Nothing in it is lost — the `tour` button replays it.'
      : 'Next card (→)';
    // Kept on the LAST card too. A reader who wants out on the final step
    // should not have to find the × to get it, and "Skip" on a card with
    // nothing after it still means "I am done reading".
    skipBtn.textContent = last ? 'Close' : 'Skip tour';
    restartBtn.classList.toggle('hidden', !(resumed && i > 0));

    writeStep(i);
    positionCard(card, highlight(step.selector), step.side);
  }

  function next() {
    if (i === STEPS.length - 1) { finish(); return; }
    i = Math.min(i + 1, STEPS.length - 1);
    render();
  }
  function back() { i = Math.max(i - 1, 0); render(); }
  function restart() { i = 0; render(); }

  const titleEl = h('h3');
  // A live region, because the card does not move focus between steps (moving
  // it on every render makes Back/Next unusable from the keyboard). Without
  // this, pressing Next changes every word on screen and announces nothing.
  const bodyEl = h('div.tour-body', { 'aria-live': 'polite' });
  const dotsEl = h('div.tour-dots', { role: 'group', 'aria-label': 'Tour progress' });
  const closeX = h('button.btn.ghost.icon.sm', { type: 'button', 'aria-label': 'Close tour', onclick: close }, '×');
  const backBtn = h('button.btn.ghost', { type: 'button', onclick: back, title: 'Previous card (←)' }, 'Back');
  const nextBtn = h('button.btn.primary', { type: 'button', onclick: next }, 'Next');
  const skipBtn = h('button.btn.ghost', { type: 'button', onclick: close }, 'Skip tour');
  const restartBtn = h('button.btn.ghost.tiny.hidden', {
    type: 'button', onclick: restart,
    title: 'Go back to the first card',
  }, 'restart');

  // `role="dialog"` with aria-modal FALSE, which is the honest description:
  // this popover deliberately does not block the app behind it (that is the
  // whole point of anchoring it to real controls), so claiming modality would
  // tell a screen reader the rest of the page is inert when it is not.
  // `aria-labelledby` points at the step title, so each step announces itself.
  titleEl.id = 'tour-title';
  const card = h('div.tour-popover.notice', {
    role: 'dialog', 'aria-modal': 'false', 'aria-labelledby': 'tour-title',
  },
  h('div.tour-head', titleEl, closeX),
    bodyEl,
    h('div.frow', { style: { alignItems: 'center', gap: 'var(--sp-4)' } }, dotsEl, restartBtn),
    h('div.btn-row', { style: { marginTop: 'var(--sp-5)', justifyContent: 'space-between' } },
      h('div.btn-row', backBtn, nextBtn),
      skipBtn,
    ),
  );

  document.body.appendChild(card);
  document.addEventListener('keydown', onKey);
  window.addEventListener('resize', onReflow);
  // Capture phase: a rail's scroll event does not bubble to window, so a
  // bubble-phase listener on document never sees it.
  document.addEventListener('scroll', onReflow, true);
  render();
  // Land the keyboard on the tour once, at creation — not on every render,
  // which would yank focus back on each step and make Back/Next unusable from
  // the keyboard. Without this a keyboard user arriving at a first run had to
  // tab past the entire application (445 focusable elements were counted on
  // this workstation) before reaching `Next` or `Skip tour`. Focus is handed
  // back to `opener` on close, so this is a loan and not a seizure.
  queueMicrotask(() => nextBtn.focus());

  self = { node: card, close, destroy };
  active = self;
  return self;
}
