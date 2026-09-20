// recentRuns.js — the cross-project recent-runs strip on the workspace home.
//
// WHAT IT IS FOR. The grid above it answers "what work do I keep"; this strip
// answers "what did I just do". It is the only place on the home page that
// crosses project boundaries, which is why every entry has to say WHERE it
// belongs — a run id with no project and no specimen beside it is an orphan,
// and on a page whose whole subject is organisation that is the one thing it
// must not be.
//
// FAILED RUNS ARE NOT HIDDEN, AND THAT IS THE POINT. A recent-activity strip
// that quietly drops failures reads as a product that always works. Every
// failed run appears in place, in time order, carrying the word "failed", and
// its reason is reachable in one keypress from the row. The reason is not in
// RunCard — the contract's flattened home-page view is {run_id, project_id,
// project_name, specimen_id, specimen_label, section, drug, status,
// created_at, finished_at, runtime_s} and deliberately carries no error text —
// so the row fetches it on demand from GET /api/runs/{id}, whose result_view
// carries `error` (app/server/jobs.py). The strip never invents one: in stub
// mode it refuses to fetch at all and says so, because js/mock.js would answer
// "why did this run fail" with a seeded PRNG.
//
// NO ORDERING BEYOND TIME. The strip is newest-first and says so on screen.
// It does not sort by runtime, by any metric, or by compound, and it shows no
// metric at all — an activity list that starts ranking outcomes is a
// cross-compound ordering wearing a different hat, and that is withdrawn.

import { h } from '../util/dom.js';
import { dur, ts, int } from '../util/fmt.js';
import { relTime, epochMs } from './projectCard.js';
import { getRun, api as apiState } from '../api.js';

/** Statuses app/server/jobs.py actually produces: TERMINAL is
 *  ("done","failed","cancelled") and a job is "queued" or "running" before
 *  that. components.css ships .chip.st-* for exactly these five. Anything
 *  else is printed verbatim on a neutral chip rather than coloured as a
 *  status it might not be. */
const KNOWN_STATUS = ['queued', 'running', 'done', 'failed', 'cancelled'];

/**
 * Where a run lives in the Studio.
 *
 * `?project=` is the contract (WORKSPACE_CONTRACT.md §8). `&run=` is a
 * CONVENTION this strip introduces and the Studio's boot must honour for the
 * deep link to land on the run rather than merely inside the project; it is
 * defined here, once, and can be replaced wholesale by passing
 * `handlers.hrefForRun`.
 */
export function studioHrefForRun(run) {
  const q = [];
  if (run && run.project_id) q.push('project=' + encodeURIComponent(String(run.project_id)));
  if (run && run.run_id) q.push('run=' + encodeURIComponent(String(run.run_id)));
  return '/studio' + (q.length ? '?' + q.join('&') : '');
}

const safeId = (s) => String(s == null ? 'unknown' : s).replace(/[^A-Za-z0-9_-]/g, '-');

/* --------------------------------------------------------------- a row */

function whenLine(run) {
  const fin = epochMs(run.finished_at);
  const made = epochMs(run.created_at);
  if (fin !== null) {
    return { text: 'finished ' + relTime(fin), abs: ts(fin / 1000), label: 'finished_at' };
  }
  if (made !== null) {
    return { text: 'submitted ' + relTime(made), abs: ts(made / 1000), label: 'created_at' };
  }
  return { text: 'no timestamp reported', abs: null, label: '' };
}

function reasonBlock(run, handlers, uid) {
  const id = 'rr-reason-' + uid;
  const body = h('div.rr-reason', { id, hidden: true });
  let loaded = false;

  const btn = h('button.btn.ghost.sm.rr-why', {
    type: 'button',
    'aria-expanded': 'false',
    'aria-controls': id,
    title: 'What the server said about this failure. Fetched from '
      + 'GET /api/runs/' + String(run.run_id) + ' when you ask for it.',
  }, 'why it failed');

  function show(...kids) {
    body.textContent = '';
    for (const k of kids) if (k) body.appendChild(k);
  }

  function load() {
    if (loaded) return;
    loaded = true;
    if (apiState && apiState.stub) {
      // A fabricated reason is worse than no reason. Say which one this is.
      show(h('p.rr-reason-text',
        'This session is showing stub data, so the real failure reason is not '
        + 'reachable from here. It is on the server that ran it: read '
        + 'GET /api/runs/' + String(run.run_id) + ', or the server’s stderr.'));
      return;
    }
    show(h('p.rr-reason-pending', 'Asking the server…'));
    Promise.resolve()
      .then(() => (typeof handlers.reasonForRun === 'function'
        ? handlers.reasonForRun(run)
        : getRun(run.run_id).then((r) => (r && r.error) || null)))
      .then((said) => {
        const text = said === null || said === undefined ? '' : String(said).trim();
        if (text) {
          show(h('p.rr-reason-text', text),
            h('p.rr-reason-src', 'What the server said, verbatim, for run '
              + String(run.run_id) + '. Quote it to whoever runs this server.'));
        } else {
          show(h('p.rr-reason-text',
            'The server recorded no reason for this failure. Run '
            + String(run.run_id) + ' is failed with an empty `error` field; '
            + 'the server’s stderr is the next place to look.'));
        }
      })
      .catch((err) => {
        loaded = false;   // a failed lookup must be retryable
        show(h('p.rr-reason-text',
          'Could not reach the server for this run’s reason: '
          + (err && err.message ? String(err.message) : 'the request failed')
          + '. Press again to retry.'));
      });
  }

  btn.addEventListener('click', () => {
    const opening = body.hidden;
    body.hidden = !opening;
    btn.setAttribute('aria-expanded', opening ? 'true' : 'false');
    if (opening) load();
  });

  return { btn, body };
}

function runRow(run, handlers) {
  const r = (run && typeof run === 'object') ? run : {};
  const rid = r.run_id == null ? '' : String(r.run_id);
  const uid = safeId(rid || Math.random().toString(36).slice(2));
  const status = typeof r.status === 'string' && r.status.trim() ? r.status.trim() : null;
  const failed = status === 'failed';
  const chipCls = status && KNOWN_STATUS.indexOf(status) >= 0 ? '.st-' + status : '';

  const href = typeof handlers.hrefForRun === 'function'
    ? handlers.hrefForRun(r) : studioHrefForRun(r);

  const idEl = rid
    ? h('a.rr-run.mono', {
      href,
      title: 'Open run ' + rid + ' in the Studio (' + href + ')',
      onclick: (ev) => {
        if (ev.defaultPrevented) return;
        if (ev.button !== 0 || ev.metaKey || ev.ctrlKey || ev.shiftKey || ev.altKey) return;
        if (typeof handlers.onOpenRun !== 'function') return;  // native navigation
        ev.preventDefault();
        try { handlers.onOpenRun(r, ev); } catch (err) {
          console.error('recentRuns: onOpenRun threw', err);
        }
      },
    }, rid)
    : h('span.rr-run.mono', 'run id not reported');

  const where = h('span.rr-where',
    h('span.rr-lbl.lbl', 'in'),
    h('span.rr-proj', {
      title: r.project_id ? 'project ' + String(r.project_id) : 'no project id reported',
    }, r.project_name ? String(r.project_name)
      : (r.project_id ? String(r.project_id) : 'no project reported')),
    h('span.rr-sep', '·'),
    h('span.rr-lbl.lbl', 'specimen'),
    h('span.rr-spec', {
      title: r.specimen_id ? 'specimen ' + String(r.specimen_id) : 'no specimen id reported',
    }, r.specimen_label ? String(r.specimen_label)
      : (r.specimen_id ? String(r.specimen_id) : 'no specimen reported')));

  const what = h('span.rr-what',
    h('span.rr-lbl.lbl', 'section'),
    h('span.rr-sect.mono', r.section ? String(r.section) : 'not reported'),
    h('span.rr-sep', '·'),
    h('span.rr-lbl.lbl', 'compound'),
    h('span.rr-drug.mono', r.drug ? String(r.drug) : 'none'));

  const w = whenLine(r);
  // Number(null) is 0, and a fabricated "ran 0.0s" on a run whose runtime the
  // server did not report is exactly the kind of zero util/fmt.js exists to
  // refuse. Absent stays absent.
  const rt = (r.runtime_s === null || r.runtime_s === undefined || r.runtime_s === ''
    || !Number.isFinite(Number(r.runtime_s))) ? null : Number(r.runtime_s);
  const when = h('span.rr-when', { title: w.abs ? w.label + ' = ' + w.abs : null },
    w.text,
    rt !== null ? h('span.rr-runtime.mono', { title: 'runtime_s' }, 'ran ' + dur(rt)) : null);

  const row = h('li.rr-row' + (failed ? '.rr-row-failed' : ''), {
    dataset: { runId: rid, status: status || '' },
  },
    h('span.rr-status',
      h('span.chip' + chipCls, {
        title: status && chipCls ? 'job status: ' + status
          : 'The server reported ' + (status
            ? 'an unrecognised status (\u201c' + status + '\u201d)'
            : 'no status') + ' for this run. It is printed as sent rather '
            + 'than coloured as a status it might not be.',
      }, status || 'status not reported')),
    idEl, where, what, when);

  if (failed) {
    const { btn, body } = reasonBlock(r, handlers, uid);
    row.appendChild(btn);
    row.appendChild(body);
  }
  return row;
}

/* ------------------------------------------------------------ the strip */

/**
 * @param {Array} runs  the contract's `recent_runs`: RunCard[], across ALL
 *   projects, newest first, <= 12. Rendered in the order given; this module
 *   does not re-order and does not filter.
 * @param {object} [handlers]
 *   onOpenRun(run, ev)    optional. Absent -> the link navigates natively.
 *   hrefForRun(run)       optional. Replaces the /studio?project=&run= URL.
 *   reasonForRun(run)     optional. Returns a string or a promise of one;
 *                         replaces the GET /api/runs/{id} lookup.
 * @returns {HTMLElement} a <section>. Never null.
 */
export function createRecentRuns(runs, handlers = {}) {
  const hs = handlers || {};
  const list = Array.isArray(runs) ? runs : null;
  const titleId = 'rr-title-h';   // an id, kept distinct from the .rr-title class

  const head = h('header.rr-head',
    h('h2.rr-title', { id: titleId }, 'Recent runs'),
    h('span.rr-sub', 'across every project, newest first'));

  const section = h('section.rr', { 'aria-labelledby': titleId }, head);

  if (list === null) {
    section.appendChild(h('p.rr-empty',
      'This server did not report recent runs. The runs themselves are '
      + 'unaffected — they live in app/runs/ — but the home page '
      + 'cannot list them. Open a project to see its own runs in the Studio.'));
    return section;
  }
  if (!list.length) {
    section.appendChild(h('p.rr-empty',
      'No runs are attached to a project yet.',
      h('span.rr-empty-next',
        'This strip lists runs that belong to a project, so a run appears '
        + 'here once you attach it to a specimen — the Studio offers that '
        + 'when a run finishes, and nothing is attached until you choose. '
        + 'Runs you never attach are not lost: they stay in the Studio\u2019s '
        + 'own run history, which lists every run on this server.')));
    return section;
  }

  const nFailed = list.filter((r) => r && r.status === 'failed').length;
  head.appendChild(h('span.rr-count.mono', {
    title: 'How many runs this strip is showing. The server sends at most 12, '
      + 'and only runs ATTACHED to a project — so this is neither the total '
      + 'in the workspace nor the total on this server. The Studio\u2019s run '
      + 'history is the complete list.',
  }, int(list.length) + ' of the most recent runs shown'));
  if (nFailed) {
    head.appendChild(h('p.rr-failed-count',
      int(nFailed) + ' of the ' + int(list.length)
      + ' runs shown failed. Each one is in the list below, marked failed, '
      + 'with its reason one press away.'));
  }

  section.appendChild(h('ol.rr-list', list.map((r) => runRow(r, hs))));
  return section;
}
