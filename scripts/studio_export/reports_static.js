// reports_static.js — the STATIC EXPORT's stand-in for js/ui/reports.js.
//
// A therapy report is not data, it is a JOB: POST /api/reports queues two runs
// per compound (the real arm and a resistance-off counterfactual), and the finished
// document is HTML rendered and served by the backend at
// /api/reports/{id}/html. There is nothing to bake — no report exists until the
// engine has run one, and the viewer is an iframe pointed at a server route.
//
// So the export ships the panel's frame with its explanation intact and NO
// controls: no "new report" button that would 501, no refresh that would fail,
// no iframe whose src cannot resolve. What is left says what the feature is and
// what it takes to use it.
//
// Same constructor signature and same returned interface as the real module,
// so js/main.js is bundled unmodified.
//
// DRIFT WARNING, because this file has already drifted once. The lead paragraph
// below is a COPY of js/ui/reports.js's, and a copy of a claim is a claim that
// can go stale on its own. It did: when compound ordering was withdrawn on
// 2026-08-21 the live module was corrected and this stand-in was not, so every
// export built after that date shipped "One specimen in, one ranked drug panel
// out" and "the ranking is the output" as CURRENT product copy. That is worse
// than the stored caveat the claim-set banner marks — a stored caveat is a
// dated record and this is the page speaking in the present tense, with no
// banner possible because nothing here is a record.
//
// tests/test_claim_set.py's retired-sentence scan does not cover this file:
// it reads .py only, and it excludes scripts/ on the reasoning that scripts do
// not render to a user. These modules do. If the live panel's wording changes
// again, change it here in the same commit.

import { h } from '../util/dom.js';

export function createReportsPanel() {
  const root = h('div.panel.rp',
    h('div.panel-head',
      h('h2', 'Therapy reports'),
      h('span.hint', 'not in this export'),
      h('span.grow')),
    h('div.panel-body.rp-body',
      h('div.rp-lead',
        h('b', 'One specimen and one compound in, one discrimination out: '
          + 'was it delivery, or was it dose?'),
        ' The compound is simulated twice on the same cached reconstruction — the real '
        + 'run, and a resistance-off counterfactual — so the O',
        h('sub', '2'),
        '-resistance term of the attribution is measured, not asserted. The document '
        + 'carries its own limits: ',
        h('b', 'the map is the output; the absolute level is not, and neither '
          + 'is a comparison against another compound.')),
      h('div.xnote',
        h('b', 'A report is a job, not a file. '),
        'It queues two engine runs per compound and is then rendered and served by the '
        + 'backend at ', h('code.badge-code', 'GET /api/reports/{id}/html'),
        '. Neither the queue nor that URL exists in a static export, so the controls '
        + 'have been removed rather than left here to fail. Start the local server '
        + 'and the panel comes back with them: ',
        h('code.badge-code', './.venv/bin/python -m app.server'), '.')),
  );

  return {
    el: root,
    dialog: null,
    setDrugs() {},
    setSection() {},
    async refresh() {},
    busy() { return false; },
    get selected() { return null; },
  };
}
