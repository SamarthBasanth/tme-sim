// reports_static.js — the STATIC EXPORT's stand-in for js/ui/reports.js.
//
// A therapy report is not data, it is a JOB: POST /api/reports queues two runs
// per drug (the real arm and a resistance-off counterfactual), and the finished
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

import { h } from '../util/dom.js';

export function createReportsPanel() {
  const root = h('div.panel.rp',
    h('div.panel-head',
      h('h2', 'Therapy reports'),
      h('span.hint', 'not in this export'),
      h('span.grow')),
    h('div.panel-body.rp-body',
      h('div.rp-lead',
        h('b', 'One specimen in, one ranked drug panel out.'),
        ' Each drug is simulated twice on the same cached reconstruction — the real '
        + 'run, and a resistance-off counterfactual — so the O',
        h('sub', '2'),
        '-resistance term of the attribution is measured, not asserted. The document '
        + 'carries its own limits: ',
        h('b', 'the ranking is the output; the absolute kill is not.')),
      h('div.xnote',
        h('b', 'A report is a job, not a file. '),
        'It queues two engine runs per drug and is then rendered and served by the '
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
