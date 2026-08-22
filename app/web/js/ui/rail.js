// rail.js — the grouping and navigation primitive shared by both rails.
//
// Both rails grew as flat stacks of <details> sections: ten of them on the
// right, 2350 px of content in a 1000 px column, with no landmark between the
// top and the bottom. Every section was legible on its own and the column as a
// whole was not navigable.
//
// This adds two things and removes nothing:
//
//   railGroup(id, title, note, ...sections)
//     a labelled band with a STICKY header, so "where am I" is always answered
//     on screen without scrolling back up;
//
//   railNav(groups, scroller)
//     a strip of group names at the top of the rail that scrolls to one.
//
// The navigator is deliberately NOT a menu: nothing lives behind it that is not
// also reachable by scrolling, and no section is closed to make room for it. It
// is a table of contents for a document that is entirely present.
//
// See docs/UI_ARCHITECTURE.md sections 2.2 and 2.4.

import { h } from '../util/dom.js';

/**
 * A labelled band inside a rail.
 * @param {string} id     stable id fragment; the DOM id becomes `rg-<id>`
 * @param {string} title  the group name, rendered upper-case by CSS
 * @param {string} note   a short right-aligned annotation, or ''
 * @param {...Node} sections the <details> sections in the group
 */
export function railGroup(id, title, note, ...sections) {
  return h('section.rail-group', { id: 'rg-' + id, 'aria-label': title },
    // <h2>, not <span>. See sectTitle() below for why the rail had no heading
    // structure at all and what that cost a screen-reader user.
    h('div.rg-head', h('h2.rg-h', title), note ? h('span.rg-note', note) : null),
    ...sections);
}

/**
 * A section title, as an actual heading.
 *
 * MEASURED 2026-08-22, on the READ rail with a finished run loaded: zero
 * elements matching `h1,h2,h3,h4,h5,h6`. Every group name was a `<span>` and
 * every section title a bare text node inside a `<summary>`, over roughly
 * 30,000 characters of disclosure in one scroller. A screen-reader user had no
 * heading list, so the only way to find "Delivery or dose" was to read from the
 * top — and the whole argument for this rail is that the limits sit next to the
 * numbers they qualify, which is a spatial claim that only pays off if you can
 * get to the place.
 *
 * `<h3>` INSIDE the `<summary>`, never `role="heading"` ON it: the summary's
 * own role is what makes the section expandable, and overwriting it would trade
 * one affordance for the other. The HTML content model for `<summary>` is
 * "phrasing content, optionally intermixed with heading content", so this is
 * the sanctioned shape rather than a trick.
 *
 * Levels: group header h2, section title h3, in-panel label h4 (see
 * js/ui/reach.js::microLabel). The page's h1 is the app title in index.html.
 */
export function sectTitle(text) {
  return h('h3.sect-h', text);
}

/**
 * The group navigator. Returns the strip; call `sync()` to repaint which group
 * is current after a scroll.
 *
 * @param {{id:string,label:string}[]} groups
 * @param {HTMLElement} scroller the scrolling rail body the groups live in
 */
export function railNav(groups, scroller) {
  const btns = groups.map((g) => h('button', {
    type: 'button',
    title: `Jump to ${g.label}`,
    'aria-current': 'false',
    onclick: () => {
      const el = scroller.querySelector('#rg-' + g.id);
      if (!el) return;
      // scrollTop rather than scrollIntoView: the latter also scrolls the
      // nearest scrollable ancestor, which here is the shell, and the shell
      // must never scroll (html/body are overflow:hidden by contract).
      //
      // Bounding rects, not offsetTop: offsetTop is measured from the nearest
      // POSITIONED ancestor, which in this shell is `.shell` (position:fixed)
      // rather than the scrolling rail body, so the arithmetic was off by the
      // rail's own origin and the navigator lit the wrong group.
      scroller.scrollTop += el.getBoundingClientRect().top
        - scroller.getBoundingClientRect().top;
      sync();
    },
  }, g.label));

  // The strip names its own verb. Without the label it is a row of group names
  // with the current one underlined in the accent — which is, to the letter,
  // the selected-tab idiom, and it was read as one: it says "RESULT is open,
  // VIEW and MEASURE are behind it" about a rail where all three are on the
  // same scroll and nothing is behind anything. Worse, the sticky group header
  // 24 px below it is the OTHER answer to "where am I", and the two can differ
  // on screen at the same instant (the navigator marks COMMIT at the bottom of
  // the rail while DEFINE's header is still stuck to the top of it) — two
  // orientation devices, adjacent, disagreeing, one of them lying about how
  // much of the panel exists. Two words fix it: this is where you can GO, and
  // the sticky header below is where you ARE.
  const el = h('nav.rail-nav', { 'aria-label': 'Section groups' },
    h('span.rn-lbl', { 'aria-hidden': 'true' }, 'go to'), btns);

  /** Mark the last group whose header has reached the top of the rail body. */
  function sync() {
    // Before the rail is appended to the document every rect is 0x0 at 0,0, so
    // every group would satisfy the test below and the LAST one would light up
    // on a rail that has not been scrolled. Nothing to say yet — say nothing.
    if (!scroller.isConnected || !scroller.clientHeight) return;

    const top = scroller.getBoundingClientRect().top;
    let cur = 0;
    groups.forEach((g, i) => {
      const n = scroller.querySelector('#rg-' + g.id);
      // +4 tolerance so a sub-pixel scroll position does not flicker between
      // two groups; the sticky header itself sits exactly at `top`.
      if (n && n.getBoundingClientRect().top <= top + 4) cur = i;
    });
    // At maximum scroll the trailing groups can never reach the top — a short
    // last group simply has nothing below it to push it up. Without this the
    // navigator says "Result" while the operator is looking at Measure, and
    // clicking Measure appears to do nothing.
    //
    // Guarded on the rail ACTUALLY overflowing: a rail whose content fits is
    // trivially "at the bottom", and marking its last group current would be
    // wrong from the first frame — which is exactly what happened on the left
    // rail before the RunSpec editor had any content in it.
    const overflows = scroller.scrollHeight > scroller.clientHeight + 2;
    if (overflows && scroller.scrollTop + scroller.clientHeight >= scroller.scrollHeight - 2) {
      cur = groups.length - 1;
    }
    btns.forEach((b, i) => b.setAttribute('aria-current', i === cur ? 'true' : 'false'));
  }

  let queued = 0;
  scroller.addEventListener('scroll', () => {
    if (queued) return;
    queued = requestAnimationFrame(() => { queued = 0; sync(); });
  }, { passive: true });

  // The rail is built detached and mounted by main.js, so the first meaningful
  // measurement is a frame away. A ResizeObserver on the scroller also catches
  // the rail being resized by its splitter, which changes which group is at the
  // top without ever firing a scroll event.
  requestAnimationFrame(sync);
  if (typeof ResizeObserver === 'function') {
    // The scroller AND every group: the rail is resized by its splitter (which
    // fires no scroll event), and the groups grow as panels fill with a run's
    // data (which also fires no scroll event and does not resize the scroller).
    const ro = new ResizeObserver(sync);
    ro.observe(scroller);
    for (const g of groups) {
      const n = scroller.querySelector('#rg-' + g.id);
      if (n) ro.observe(n);
    }
  }

  return { el, sync };
}
