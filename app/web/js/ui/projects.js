// projects.js — the project / specimen browser and the cross-specimen
// comparison view. Self-contained: it injects its own toolbar button and its
// own <dialog>, so it adds nothing to the three-rail workstation layout.
//
// WHY THE COMPARISON IS LAID OUT THIS WAY
// ---------------------------------------
// WHAT CHANGED: this view used to open on "2. RANKING — THE PRIMARY RESULT", a
// per-compound column table, a Kendall tau-b pairwise agreement table and a
// RANKING AGREES / DISAGREES verdict. All of it is gone. Ordering compounds by
// predicted effect is withdrawn: per-compound plasma exposure spans 113.1x
// across the compounds the GBM panel ranks and the published ordering breaks at
// a spread of only 1.53x (results/pk_ranking_robustness.md), and the breast
// top-1 pick changes under this repo's own structural capacity rule
// (results/capacity_ranking_robustness.md). Comparing runs did NOT go with it —
// the same compound on two specimens is a different question and a legitimate
// one, so that is the only comparison the page now draws.
//
// The order on screen is the order the evidence supports, and it is not
// negotiable by the reader:
//
//   1. SPECIMENS      — what each panel can and cannot support. A number is
//                       never shown above its own capability verdict. With (3),
//                       this is the primary result of a comparison now: the
//                       capability contrast and the tissue contrast.
//   2. SAME COMPOUND, — one row per compound both specimens ran at an identical
//      TWO SPECIMENS    RunSpec: that compound's own value on each side. Shown,
//                       but blocked wherever the backend says
//                       `comparable: false`. Absolute kill moves 0.3252
//                       absolute (1.77x) with the state route ALONE, and a panel
//                       with no hypoxia-response transcript cannot pin it down
//                       under any route (cell_state_impact.md §5). Two DIFFERENT
//                       compounds are never placed on one axis here.
//   3. COMPOSITION    — where a user actually sees that patient B is more
//                       vascularised than patient A, with the two rows that
//                       moved in OPPOSITE directions on the real pair
//                       (endothelial fraction 0.67x, areal density 4.19x —
//                       heldout_patient.md §1.1) sitting adjacent on purpose.
//
// Every verdict, reason and caveat string comes from the backend verbatim. This
// module decides layout and emphasis; it never decides comparability.
//
// WHAT THIS DIALOG IS FOR NOW (docs/WORKSPACE_CONTRACT.md §8)
// ----------------------------------------------------------
// It used to be the ONLY door to the project store, so it carried the whole of
// it: a `new project name` box, a per-row delete, and a project list that was
// the only place a project could be seen. There is a workspace home now, and
// two doors to the same lifecycle is two places for the same disclosure to
// drift out of — specifically §2's, that DELETING A PROJECT DOES NOT DELETE
// ITS RUNS, which has to be stated at the point of the click and is therefore
// stated where the click is.
//
// So the lifecycle hands off and this dialog keeps the half the home page does
// not do:
//
//   handed off to `/`   create · rename · duplicate · delete · browse
//   kept here           specimens, the run associations on each, and the
//                       cross-specimen comparison — the things you do INSIDE
//                       one project, beside the run you just produced
//
// The project list stays, because the comparison needs something to be a
// comparison OF and switching it must not cost a page load. What each row
// gained instead of a delete button is a link that reloads the Studio scoped
// to that project (`?project=<id>`), which is the one thing only this page can
// offer: it is already the Studio.

import { h, clear, fill, $ } from '../util/dom.js';
import { chip } from './provenance.js';
import { auto, num, int, pct, sci, ts, NULL_GLYPH } from '../util/fmt.js';
import * as api from '../api.js';
import { toast, notifyError } from './toast.js';
import { confirmable } from './confirm.js';

const METRIC_LABEL = {
  kill_fraction: 'Kill fraction',
  mean_tumor_auc: 'Mean tumour AUC',
  therapeutic_coverage_1x: 'Coverage @ 1x IC50',
  penetration_length_um: 'Penetration length',
};
const METRIC_UNIT = {
  kill_fraction: '%', mean_tumor_auc: 'uM.s',
  therapeutic_coverage_1x: '%', penetration_length_um: 'um',
};

function fmtBy(fmt, v) {
  if (v === null || v === undefined) return NULL_GLYPH;
  if (Array.isArray(v)) return v.length ? v.join(', ') : NULL_GLYPH;
  if (typeof v !== 'number' || Number.isNaN(v)) return String(v);
  switch (fmt) {
    case 'pct': return pct(v, 2);
    case 'int': return int(v);
    case 'sci': return sci(v, 3);
    default: return auto(v);
  }
}
const metricFmt = (m) =>
  (m === 'kill_fraction' || m === 'therapeutic_coverage_1x') ? 'pct'
    : m === 'mean_tumor_auc' ? 'sci' : 'num';

function ratioText(r) {
  if (r === null || r === undefined || !Number.isFinite(r)) return NULL_GLYPH;
  return num(r, r >= 100 ? 0 : 2) + 'x';
}

export function createProjectsPanel({ store, projectId = null }) {
  let projects = [];
  let project = null;         // the full document for the selected project
  let comparison = null;
  let metric = 'kill_fraction';
  let busy = false;
  // The project THE STUDIO ITSELF is scoped to, from `?project=` — null in the
  // ordinary no-project Studio, which is every Studio that was ever opened
  // before the workspace home existed. It is a label and a default selection
  // here, never a restriction: the dialog still shows every project, because a
  // comparison you cannot reach is not a comparison.
  let contextId = projectId || null;

  /* ---------------------------------------------------------------- shell */

  const railEl = h('div.proj-rail');
  const mainEl = h('div.proj-main');
  const titleEl = h('span.mono.strong.proj-title', { title: '' }, NULL_GLYPH);
  // Was "rank by". It no longer selects anything about an ordering — there is
  // none — so it is relabelled to what it actually does: pick which metric the
  // same-compound cross-specimen rows are featured on. Every row carries all
  // four whatever is chosen, so this changes a column, never a conclusion.
  const metricSel = h('select.inp.sm', {
    title: 'Which metric the same-compound cross-specimen rows are read on. '
      + 'It selects a column, not a conclusion: every row carries all four '
      + 'metrics regardless, and no metric produces an ordering of compounds.',
    onchange: (e) => { metric = e.target.value; refreshComparison(); },
  }, Object.keys(METRIC_LABEL).map((k) =>
    h('option', { value: k, selected: k === metric }, METRIC_LABEL[k])));

  /** A link out to the workspace home, where a project's lifecycle lives. */
  function homeLink(label, title) {
    return h('a.btn.ghost.sm', {
      href: '/',
      // `.btn` was written for <button>; an <a> would otherwise arrive
      // underlined and in the UA's link blue, the one colour tokens.css does
      // not own.
      style: { textDecoration: 'none' },
      title,
    }, label);
  }

  /**
   * A link that reloads the Studio scoped to a project.
   *
   * The href is built from the CURRENT path on purpose. This page is served at
   * `/` today and at `/studio` once the workspace contract's routes land
   * (§1), and a hard-coded path would be wrong on exactly one of those two
   * days. `?project=` is the whole of the Studio's project state, so the same
   * link is also what a bookmark of this view would be.
   */
  function studioLink(id) {
    return h('a.btn.ghost.sm', {
      href: `${window.location.pathname}?project=${encodeURIComponent(id)}`,
      style: { textDecoration: 'none' },
      title: `Reload the Studio scoped to ${id}. A run you finish there is `
        + 'OFFERED for filing under one of its specimens when it completes — '
        + 'offered, never filed for you.',
      'aria-label': `Open project ${id} in the Studio`,
      // The row behind this link is itself a click target that SELECTS the
      // project. Without this, following the link would also fire a select
      // and a network round-trip on a page that is being navigated away from.
      onclick: (e) => e.stopPropagation(),
    }, 'open in Studio');
  }

  const dlg = h('dialog#proj-dlg',
    h('div.panel-head',
      h('h2', 'Projects'),
      // Was "multi-patient studies". The grouping key is the SPECIMEN — a
      // string the operator chose, which must not be a patient identifier —
      // and the rail two hundred lines below already refuses to name the unit
      // "patient" for exactly that reason. The header was the last place in
      // this file still doing it.
      h('span.hint', 'specimens and the cross-specimen comparison'),
      h('span.vrule'),
      h('span.lbl', 'project'), titleEl,
      h('span.grow'),
      h('span.lbl', 'read on'), metricSel,
      h('button.btn.ghost.sm', { onclick: () => reload(), title: 'GET /api/projects' }, 'refresh'),
      homeLink('workspace home \u2192',
        'Create, rename, duplicate and delete projects on the workspace home. '
        + 'This dialog no longer does any of those: one place for a lifecycle '
        + 'is one place for its disclosures.'),
      h('button.btn.ghost.sm', { id: 'proj-close', onclick: () => dlg.close() }, 'close')),
    h('div.proj-body', railEl, h('div.proj-div'), mainEl));

  document.body.appendChild(dlg);
  dlg.addEventListener('close', () => { /* keep state for the next open */ });

  /* ------------------------------------------------------------- toolbar */

  const btn = h('button.btn.ghost.sm', {
    id: 'btn-projects',
    // DELIBERATELY UNCHANGED, tooltip included. The top bar is the Studio's
    // persistent chrome — the part you see without opening anything — and
    // docs/WORKSPACE_CONTRACT.md §8 says a user who never visits the
    // workspace home notices no change. Re-pointing this dialog (see the note
    // at the top of this file) necessarily changes what is INSIDE it; it does
    // not get to change the shell on the way. So the hand-off is stated in the
    // dialog's own head row and rail, where the reader who opened it is
    // already looking, and this string stays byte-for-byte what it was.
    title: 'Projects: group runs by patient and compare across specimens',
    onclick: open,
  }, 'projects');
  const slot = $('.topbar-group.last');
  if (slot) slot.insertBefore(btn, slot.firstChild);

  /* ----------------------------------------------------------- data flow */

  function offline() {
    return !!(api.api && api.api.stub);
  }

  async function open() {
    if (!dlg.open) dlg.showModal();
    await reload();
  }

  async function reload() {
    if (offline()) {
      renderOffline();
      return;
    }
    try {
      projects = await api.listProjects();
    } catch (err) {
      notifyError('GET /api/projects failed', err);
      projects = [];
    }
    if (project && !projects.some((p) => p.id === project.id)) project = null;
    // Open on the project THE STUDIO IS IN, when it is in one. The dialog and
    // the top bar naming different projects is precisely the "unsure which
    // project this run lands in" state the contract's §8 exists to prevent,
    // and the newest-first first row is not a reason to create it.
    const preferred = contextId && projects.some((p) => p.id === contextId)
      ? contextId : (projects.length ? projects[0].id : null);
    if (!project && preferred) await select(preferred);
    else if (project) await select(project.id);
    else { comparison = null; render(); }
  }

  async function select(id) {
    try {
      project = await api.getProject(id);
    } catch (err) {
      notifyError(`GET /api/projects/${id} failed`, err);
      project = null;
    }
    render();
    await refreshComparison();
  }

  async function refreshComparison() {
    if (!project) { comparison = null; render(); return; }
    try {
      comparison = await api.compareProject(project.id, metric);
    } catch (err) {
      comparison = null;
      notifyError(`GET /api/projects/${project.id}/compare failed`, err);
    }
    render();
  }

  async function guard(label, fn) {
    if (busy) return null;
    busy = true;
    try {
      return await fn();
    } catch (err) {
      notifyError(label + ' failed', err);
      return null;
    } finally {
      busy = false;
    }
  }

  /* -------------------------------------------------------------- render */

  function renderOffline() {
    clear(railEl);
    fill(mainEl, h('div.notice', { style: { maxWidth: '520px' } },
      h('h3', 'Projects need the live backend'),
      h('p', 'A project is state you built — specimens, and the run ids that ' +
        'belong to each patient. js/mock.js deliberately does not synthesise ' +
        'one, because a fabricated study on screen is worse than no study.'),
      h('p.mono.dim', { style: { marginTop: '8px' } },
        './.venv/bin/python -m app.server')));
    titleEl.textContent = NULL_GLYPH;
    titleEl.title = '';
  }

  function render() {
    renderRail();
    renderMain();
    const label = project ? `${project.id} · ${project.name}` : NULL_GLYPH;
    titleEl.textContent = label;
    // The head row truncates this with an ellipsis (a long project name plus
    // the read-on control and the refresh/close buttons do not all fit on
    // one --h-head line at the dialog's min width) — the full name is still
    // one hover away.
    titleEl.title = project ? label : '';
  }

  /* --- rail: projects, specimens, run association ---------------------- */

  function renderRail() {
    clear(railEl);

    /* the Studio's own scope, when it has one — stated before the list, so
       the list is read as "and here is everywhere else", not as a menu of
       equals one of which is silently special. */
    if (contextId) {
      const inList = projects.find((p) => p.id === contextId);
      railEl.appendChild(h('div.caveat.info', h('span',
        h('b', 'This Studio is scoped to '),
        h('span.mono.strong', contextId),
        inList && inList.name ? ` · ${inList.name}` : '',
        '. A run you finish in this tab is offered for filing under one of '
        + 'its specimens — offered, never filed for you. Selecting another '
        + 'project below changes what is COMPARED here; it does not change '
        + 'where this tab\'s runs are offered. Use "open in Studio" for that.')));
    }

    /* projects */
    railEl.appendChild(h('div.hr-lbl', `Projects · ${projects.length}`));
    // WHAT CHANGED: a `new project name` box and a `create` button stood here,
    // and every row carried a delete. Both are the workspace home's now — see
    // the note at the top of this file. The link is not a consolation prize
    // for a removed feature; it is the removal being honest about where the
    // feature went.
    railEl.appendChild(h('div.proj-new',
      homeLink('workspace home \u2192',
        'Projects are created, renamed, duplicated and deleted on the '
        + 'workspace home, where deleting one can say what it does and does '
        + 'not delete at the point of the click: a project owns run IDS, so '
        + 'deleting it never deletes a run.'),
      h('span.hint', { style: { minWidth: 0, overflow: 'hidden', textOverflow: 'ellipsis' } },
        'create · rename · duplicate · delete')));

    if (!projects.length) {
      railEl.appendChild(h('div.caveat.info', h('span',
        // "by PATIENT" until 2026-08-15. Vivum is scoped to de-identified
        // tissue (docs/DEPLOYMENT.md §0), and the grouping key is the specimen
        // anyway — a string the operator chose, which must not be a patient
        // identifier. A UI that names the unit "patient" invites one.
        'No projects yet. A project groups runs by SPECIMEN so a compound ' +
        'panel on one specimen can be compared against the same panel on ' +
        'another. Create the first one on the workspace home — the link above.')));
    }
    for (const p of projects) {
      const on = project && p.id === project.id;
      const here = contextId === p.id;
      railEl.appendChild(h('div.proj-item', {
        class: on ? 'on' : null,
        tabindex: '0',
        onclick: () => select(p.id),
        // `e.target !== e.currentTarget` because the row now contains its own
        // focusable control. Without the guard, Enter on the link inside would
        // ALSO be swallowed by preventDefault here and turn a navigation into
        // a selection — the same class of mis-aimed-click defect the delete
        // button that used to sit here was given a two-step confirm for.
        onkeydown: (e) => {
          if (e.target !== e.currentTarget) return;
          if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); select(p.id); }
        },
      },
        h('div.proj-item-top',
          h('span.mono.strong', p.id),
          h('span.grow', p.name || NULL_GLYPH),
          // Deleting a project moved to the workspace home. It is not that the
          // button was dangerous here — it had a two-step confirm — but that
          // §2's consequence ("deleting a project does not delete its runs")
          // has to be stated at the point of deletion, and it can only be
          // stated in one place without drifting.
          here
            ? h('span.chip.st-running', {
                title: 'This tab is already scoped to this project. A run that '
                  + 'finishes here is offered for filing under one of its '
                  + 'specimens.',
              }, 'this tab')
            : studioLink(p.id)),
        h('div.proj-item-sub.mono',
          `${int(p.n_specimens)} specimen${p.n_specimens === 1 ? '' : 's'} · `,
          `${int(p.n_runs)} run${p.n_runs === 1 ? '' : 's'} · `,
          ts(p.created_at))));
    }

    if (!project) return;

    /* specimens */
    const specimens = project.specimens || [];
    railEl.appendChild(h('div.hr-lbl', `Specimens · ${specimens.length}`));

    const sections = (store.get().sections || []);
    const used = new Set(specimens.map((s) => s.section));
    const free = sections.filter((s) => !used.has(s.id));
    const sectionSel = h('select.inp.sm', { 'aria-label': 'section to add as a specimen' },
      free.length
        ? free.map((s) => h('option', {
            value: s.id, disabled: !s.available,
            title: s.available ? s.panel : s.reason,
          }, `${s.id} — ${s.label}${s.available ? '' : ' (not cached)'}`))
        : [h('option', { value: '' }, 'every known section is already a specimen')]);
    railEl.appendChild(h('div.proj-new',
      sectionSel,
      h('button.btn.sm', {
        disabled: !free.length,
        onclick: async () => {
          const sid = sectionSel.value;
          if (!sid) return;
          const out = await guard('POST /specimens',
            () => api.addSpecimen(project.id, sid, null));
          if (!out) return;
          toast('ok', `Specimen ${out.id} added`, `${sid} — ${out.label}`, 3000);
          await select(project.id);
        },
      }, 'add specimen')));

    for (const sp of specimens) railEl.appendChild(specimenCard(sp));
  }

  function specimenCard(sp) {
    const cap = sp.capability || {};
    const runs = store.get().runs || [];
    const attached = sp.run_ids || [];
    const candidates = runs.filter(
      (r) => r.spec && r.spec.section === sp.section && !attached.includes(r.run_id));

    const runSel = h('select.inp.sm', { 'aria-label': `run to attach to ${sp.id}` },
      candidates.length
        ? candidates.map((r) => h('option', { value: r.run_id },
            `${r.run_id} · ${(r.spec && (r.spec.drug ?? 'none'))} · ${r.status}`))
        : [h('option', { value: '' }, `no unattached run on ${sp.section}`)]);

    // Capability is a property of the PANEL, not a provenance claim about a
    // number, so it does not wear a provenance chip.
    const capChip = cap.absolute_kill_supported
      ? h('span.chip.st-done', { title: cap.absolute_kill_reason }, 'state supported')
      : h('span.chip.st-failed', { title: cap.absolute_kill_reason }, 'no absolute kill');

    return h('div.proj-item.spec',
      h('div.proj-item-top',
        h('span.mono.strong', sp.id),
        h('span.grow', sp.label || NULL_GLYPH),
        capChip,
        confirmable(h('button.btn.ghost.sm.danger', {
          title: 'Remove this specimen. Its runs are NOT deleted.',
          'aria-label': `Remove specimen ${sp.id}`,
        }, '×'), {
          confirmLabel: `remove ${sp.id}?`,
          armedTitle: `Removes ${sp.id} from this study, with the run `
            + 'associations it carries. Its runs stay in app/runs/. Press again '
            + 'to confirm; Escape cancels.',
          onConfirm: async () => {
            const out = await guard('DELETE specimen',
              () => api.deleteSpecimen(project.id, sp.id));
            if (!out) return;
            toast('ok', `${sp.id} removed`,
              `${(out.runs_retained || []).length} run(s) retained in app/runs/.`, 4000);
            await select(project.id);
          },
        })),
      h('div.proj-item-sub.mono',
        `${sp.section} · ${cap.panel || NULL_GLYPH} · ${int(cap.n_cells)} cells`),
      cap.hypoxia_present && cap.hypoxia_present.length
        ? h('div.proj-item-sub.dim.mono', `hypoxia markers: ${cap.hypoxia_present.join(', ')}`)
        : h('div.proj-item-sub.mono', { style: { color: 'var(--warn)' } },
            `no hypoxia-response transcript on this panel (${(cap.hypoxia_missing || []).join(', ')})`),
      h('div.proj-runs',
        attached.length
          ? attached.map((rid) => h('span.runchip',
              h('span.mono', rid),
              h('button.btn.ghost.sm', {
                title: 'Detach. The run itself is untouched.',
                'aria-label': `Detach run ${rid}`,
                onclick: async () => {
                  const out = await guard('DELETE run association',
                    () => api.detachRun(project.id, sp.id, rid));
                  if (!out) return;
                  await select(project.id);
                },
              }, '×')))
          : h('span.dim', 'no runs attached')),
      h('div.proj-new',
        runSel,
        h('button.btn.sm', {
          disabled: !candidates.length,
          onclick: async () => {
            const rid = runSel.value;
            if (!rid) return;
            const out = await guard('POST run association',
              () => api.attachRun(project.id, sp.id, rid));
            if (!out) return;
            toast(out.attached ? 'ok' : 'warn',
              out.attached ? `${rid} attached to ${sp.id}` : `${rid} was already attached`,
              'Runs are associated, never copied — app/runs/ stays the only store.',
              3000);
            await select(project.id);
          },
        }, 'attach run')));
  }

  /* --- main: the comparison -------------------------------------------- */

  function renderMain() {
    clear(mainEl);
    if (!project) {
      mainEl.appendChild(h('div.notice',
        h('h3', 'No project selected'),
        h('p', 'Create a project on the workspace home, then come back here, ' +
          'add two sections as specimens, and run the SAME compound at the ' +
          'SAME spec on each. The comparison then puts that one compound side ' +
          'by side across the two specimens, and the two tissues side by side ' +
          'under it. It does not order compounds against each other.')));
      return;
    }
    if (!comparison) {
      mainEl.appendChild(h('div.caveat.info', h('span',
        'No comparison yet — GET /api/projects/{id}/compare has not returned.')));
      return;
    }
    const c = comparison;
    const sids = c.specimens.map((s) => s.id);
    const label = (sid) => {
      const s = c.specimens.find((x) => x.id === sid);
      return s ? s.label : sid;
    };

    /* 0. headline — what this comparison established. Never an ordering: the
       tag and its severity are the backend's, so the browser cannot upgrade a
       "not comparable" into a finding. */
    mainEl.appendChild(h('div.verdict',
      { class: c.headline_level === 'ok' ? 'ok'
        : c.headline_level === 'crit' ? 'crit' : 'info' },
      h('div.verdict-tag', c.headline_tag || 'NOT EVALUABLE'),
      h('div.verdict-body', c.headline)));
    if (c.compares) {
      mainEl.appendChild(h('div.caveat.crit', h('span',
        h('b', 'WHAT THIS COMPARES: '), c.compares)));
    }

    /* 1. specimens */
    mainEl.appendChild(h('div.hr-lbl', '1 · Specimens — what each panel can support'));
    mainEl.appendChild(h('div.tblwrap', h('table.tbl.zebra',
      h('thead', h('tr',
        h('th', 'specimen'), h('th', 'section'), h('th', 'panel'),
        h('th.right', 'cells'), h('th.right', 'runs'),
        h('th', 'hypoxia on panel'), h('th', 'absolute kill'))),
      h('tbody', c.specimens.map((s) => {
        const cap = s.capability || {};
        return h('tr',
          h('td', h('span.mono.strong', s.id), ' ', s.label),
          h('td.mono', s.section),
          h('td', { title: cap.accession || '' }, cap.panel || NULL_GLYPH),
          h('td.right.num', int(cap.n_cells)),
          h('td.right.num', `${int(s.n_runs_done)} / ${int(s.n_runs_attached)}`),
          h('td.mono', { style: { color: (cap.hypoxia_present || []).length ? null : 'var(--warn)' } },
            (cap.hypoxia_present || []).length
              ? cap.hypoxia_present.join(', ') : 'none'),
          h('td', s.absolute_kill_supported
            ? h('span.chip.st-done', { title: cap.absolute_kill_reason }, 'supported')
            : h('span.chip.st-failed', { title: cap.absolute_kill_reason }, 'no abs kill')));
      })))));
    for (const s of c.specimens) {
      for (const lab of s.labels || []) {
        mainEl.appendChild(h('div.caveat.crit',
          h('span', h('b', `${s.label} · ${s.section}: `), lab, ' — ',
            h('span.dim', s.absolute_kill_reason || ''))));
      }
      for (const nd of s.runs_not_done || []) {
        mainEl.appendChild(h('div.caveat', h('span',
          `${s.label}: run ${nd.run_id} is ${nd.status} and is excluded from the comparison.`,
          nd.error ? ` ${nd.error}` : '')));
      }
    }
    for (const p of c.problems || []) {
      mainEl.appendChild(h('div.caveat.crit', h('span', p)));
    }

    /* 2. the same compound on two specimens — flagged */
    mainEl.appendChild(h('div.hr-lbl',
      `2 · Same compound, two specimens · ${METRIC_LABEL[c.metric] || c.metric} — read ONLY where comparable`));
    if (!(c.kill || []).length) {
      mainEl.appendChild(h('div.caveat.info', h('span',
        'No compound is matched across two specimens at an identical RunSpec. ' +
        'Two runs are matched only when everything but the section and the ' +
        'compound is the same — same mesh, same modes, same PK, same window.')));
    } else {
      const ref = c.kill[0].reference;
      // Each ROW holds the compound fixed and varies the specimen. Rows are not
      // an ordering and must never be read down the column as one: the backend
      // emits them alphabetically by compound, not by value.
      mainEl.appendChild(h('div.caveat.info', h('span',
        'One row per compound, sorted by NAME. Each row compares that compound ' +
        'against ITSELF on the other specimen; the Δ and the ratio are ' +
        'cross-specimen. Rows are not ranked against each other and reading ' +
        'down a column as an ordering is exactly the output that was withdrawn.')));
      mainEl.appendChild(h('div.tblwrap', h('table.tbl.zebra',
        h('thead', h('tr',
          h('th', 'compound'),
          ...sids.map((sid) => h('th.right', { title: label(sid) }, sid)),
          h('th.right', 'Δ abs'), h('th.right', 'ratio'),
          h('th', 'comparable'))),
        h('tbody', c.kill.map((row) => {
          const m = row.metrics[c.metric] || { values: {}, delta: {}, ratio: {} };
          const other = sids.find((s) => s !== row.reference);
          return h('tr', { class: row.comparable ? null : 'blocked' },
            h('td.mono.strong', row.drug),
            ...sids.map((sid) => h('td.right.num',
              fmtBy(metricFmt(c.metric), m.values[sid]),
              h('span.unit', METRIC_UNIT[c.metric] || ''))),
            h('td.right.num', other === undefined ? NULL_GLYPH
              : fmtBy(metricFmt(c.metric), m.delta[other])),
            h('td.right.num', other === undefined ? NULL_GLYPH
              : ratioText(m.ratio[other])),
            h('td', row.comparable
              ? h('span.chip.st-done', { title: 'Both panels support a measured cell state and both runs used one. Comparable is a statement about the PAIR, not about the provenance of either number.' }, 'yes')
              : h('span.chip.st-failed', { title: (row.blockers || []).join('\n\n') }, 'NO')));
        })))));
      const blockers = [...new Set(c.kill.flatMap((r) => r.blockers || []))];
      for (const b of blockers) {
        mainEl.appendChild(h('div.caveat.crit', h('span', h('b', 'NOT COMPARABLE: '), b)));
      }
      const caveats = [...new Set(c.kill.flatMap((r) => r.caveats || []))];
      for (const cv of caveats) {
        mainEl.appendChild(h('div.caveat', h('span', cv)));
      }
      mainEl.appendChild(h('div.hint.mono', { style: { padding: '2px 0' } },
        `runs: ${c.kill.map((r) => `${r.drug} [${sids.map((s) => r.runs[s] || '?').join(', ')}]`).join('  ')}`));
      if (ref) {
        mainEl.appendChild(h('div.hint.mono', { style: { padding: '2px 0' } },
          `deltas and ratios are against ${ref} · ${label(ref)}`));
      }
    }
    // The transvascular delivery law produced every delivery number above. It
    // is compound-specific — it changes the top-1 GBM pick — so it does not
    // cancel between two specimens the way a common-mode perturbation does, and
    // it does not cancel out of an absolute level either. This caveat used to
    // hang under the ranking table; the ranking is gone and the disclosure is
    // not. Backend text, read from PDEConfig at request time
    // (app/server/delivery_law.py).
    if (c.delivery_law && c.delivery_law.caveat) {
      mainEl.appendChild(h('div.caveat.crit', { title: c.delivery_law.caveat },
        h('span', h('b', 'LAW-DEPENDENT: '), c.delivery_law.caveat)));
    }

    /* 3. composition */
    const comp = c.composition || { rows: [] };
    mainEl.appendChild(h('div.hr-lbl', '3 · Composition — the tissue difference'));
    if (!comp.rows.length) {
      mainEl.appendChild(h('div.caveat.info', h('span',
        'No completed run on any specimen yet, so there is no reconstruction ' +
        'report to read composition from.')));
    } else {
      const cs = comp.specimens || sids;
      mainEl.appendChild(h('div.tblwrap', h('table.tbl.zebra',
        h('thead', h('tr',
          h('th', 'quantity'),
          ...cs.flatMap((sid) => [
            h('th.right', { title: label(sid) }, sid),
            h('th', 'prov')]),
          h('th.right', 'ratio'), h('th', 'comparable'))),
        h('tbody', comp.rows.map((row) => {
          const other = cs.find((s) => s !== row.reference);
          return h('tr', { class: row.comparable ? null : 'blocked' },
            h('td', row.label, h('span.unit', row.unit || '')),
            ...cs.flatMap((sid) => {
              const p = (row.provenance || {})[sid] || {};
              return [
                h('td.right.num', fmtBy(row.fmt, (row.values || {})[sid])),
                h('td', chip(p.kind, p.why)),
              ];
            }),
            h('td.right.num', other === undefined ? NULL_GLYPH
              : ratioText((row.ratio || {})[other])),
            h('td', row.comparable
              ? h('span.chip.measured', { title: row.reason }, 'yes')
              : h('span.chip.unknown', { title: row.reason }, 'NO')));
        })))));
      const reasons = [...new Set(comp.rows.filter((r) => !r.comparable)
        .map((r) => `${r.label}: ${r.reason}`))];
      for (const r of reasons) {
        mainEl.appendChild(h('div.caveat', h('span', r)));
      }
      mainEl.appendChild(h('div.caveat.info', h('span', comp.note)));
      mainEl.appendChild(h('div.hint.mono', { style: { padding: '2px 0' } },
        'read from ' + cs.map((s) => `${s}=${(comp.from_run || {})[s] || '?'}`).join('  ')));
    }

    /* matched sets, last — provenance for the whole page */
    if ((c.matched_sets || []).length) {
      mainEl.appendChild(h('div.hr-lbl', 'Matched sets — the RunSpec signatures compared'));
      mainEl.appendChild(h('div.tblwrap', h('table.tbl',
        h('thead', h('tr', h('th', 'signature'), h('th', 'specimens'),
          h('th', 'compounds'), h('th.right', 'matched'))),
        h('tbody', c.matched_sets.map((ms) => h('tr',
          h('td.mono', { style: { fontSize: 'var(--fs-micro)' } },
            ms.is_primary ? h('b', 'PRIMARY · ') : null, ms.signature),
          h('td.mono', ms.specimens.join(', ')),
          h('td', ms.drugs.join(', ')),
          h('td.right.num', int(ms.n_matched_drugs))))))));
    }
  }

  render();

  /**
   * Tell the dialog which project the STUDIO is scoped to.
   *
   * Called by js/main.js once `?project=` has actually resolved to a project —
   * never for an id that 404'd, because a context line naming a project that
   * does not exist is worse than no context line.
   */
  function setContext(id) {
    contextId = id || null;
    // Only repaint what is already on screen. The dialog fetches on open, so a
    // context set before it was ever opened is picked up by the first reload().
    if (dlg.open) render();
  }

  return {
    el: dlg,
    button: btn,
    open,
    reload,
    setContext,
    get projectId() { return project ? project.id : null; },
    get contextId() { return contextId; },
  };
}
