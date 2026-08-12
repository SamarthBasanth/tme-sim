// projects.js — the project / specimen browser and the cross-specimen
// comparison view. Self-contained: it injects its own toolbar button and its
// own <dialog>, so it adds nothing to the three-rail workstation layout.
//
// WHY THE COMPARISON IS LAID OUT THIS WAY
// ---------------------------------------
// The order on screen is the order the evidence supports, and it is not
// negotiable by the reader:
//
//   1. SPECIMENS      — what each panel can and cannot support. A number is
//                       never shown above its own capability verdict.
//   2. RANKING        — THE PRIMARY RESULT, and the only one two real patients
//                       have been shown to agree on (heldout_patient.md §3.3;
//                       cell_state_impact.md §3.3 found it invariant across all
//                       three cell-state routes as well). Both axes are
//                       COMMON-MODE. The KILL ranking is NOT invariant across
//                       the vascular delivery law, which changes the GBM top-1
//                       (transvascular_exchange_revalidation.md §4.1), and
//                       delivery-rank invariance does not imply kill-rank
//                       invariance (breast: AUC rho 1.000, kill rho 0.830).
//   3. ABSOLUTE KILL  — shown, but struck through and blocked wherever the
//                       backend says `comparable: false`. Absolute kill moves
//                       0.3252 absolute (1.77x) with the state route ALONE, and
//                       a panel with no hypoxia-response transcript cannot pin
//                       it down under any route (cell_state_impact.md §5).
//   4. COMPOSITION    — where a user actually sees that patient B is more
//                       vascularised than patient A, with the two rows that
//                       moved in OPPOSITE directions on the real pair
//                       (endothelial fraction 0.43x, areal density 2.67x —
//                       heldout_patient.md §1.1) sitting adjacent on purpose.
//
// Every verdict, reason and caveat string comes from the backend verbatim. This
// module decides layout and emphasis; it never decides comparability.

import { h, clear, fill, $ } from '../util/dom.js';
import { chip } from './provenance.js';
import { auto, num, int, pct, sci, ts, NULL_GLYPH } from '../util/fmt.js';
import * as api from '../api.js';
import { toast, notifyError } from './toast.js';

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

export function createProjectsPanel({ store }) {
  let projects = [];
  let project = null;         // the full document for the selected project
  let comparison = null;
  let metric = 'kill_fraction';
  let busy = false;

  /* ---------------------------------------------------------------- shell */

  const railEl = h('div.proj-rail');
  const mainEl = h('div.proj-main');
  const titleEl = h('span.mono.strong.proj-title', { title: '' }, NULL_GLYPH);
  const metricSel = h('select.inp.sm', {
    title: 'Which metric the RANKING is taken over',
    onchange: (e) => { metric = e.target.value; refreshComparison(); },
  }, Object.keys(METRIC_LABEL).map((k) =>
    h('option', { value: k, selected: k === metric }, METRIC_LABEL[k])));

  const dlg = h('dialog#proj-dlg',
    h('div.panel-head',
      h('h2', 'Projects'),
      h('span.hint', 'multi-patient studies'),
      h('span.vrule'),
      h('span.lbl', 'project'), titleEl,
      h('span.grow'),
      h('span.lbl', 'rank by'), metricSel,
      h('button.btn.ghost.sm', { onclick: () => reload(), title: 'GET /api/projects' }, 'refresh'),
      h('button.btn.ghost.sm', { id: 'proj-close', onclick: () => dlg.close() }, 'close')),
    h('div.proj-body', railEl, h('div.proj-div'), mainEl));

  document.body.appendChild(dlg);
  dlg.addEventListener('close', () => { /* keep state for the next open */ });

  /* ------------------------------------------------------------- toolbar */

  const btn = h('button.btn.ghost.sm', {
    id: 'btn-projects',
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
    if (!project && projects.length) await select(projects[0].id);
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
    // the rank-by control and the refresh/close buttons do not all fit on
    // one --h-head line at the dialog's min width) — the full name is still
    // one hover away.
    titleEl.title = project ? label : '';
  }

  /* --- rail: projects, specimens, run association ---------------------- */

  function renderRail() {
    clear(railEl);

    /* projects */
    railEl.appendChild(h('div.hr-lbl', `Projects · ${projects.length}`));
    const nameIn = h('input.inp.sm', { placeholder: 'new project name',
      'aria-label': 'new project name' });
    railEl.appendChild(h('div.proj-new',
      nameIn,
      h('button.btn.sm.primary', {
        onclick: async () => {
          const name = nameIn.value.trim();
          if (!name) { toast('warn', 'Name the project', 'A study needs a name before it can hold specimens.'); return; }
          const p = await guard('POST /api/projects',
            () => api.createProject(name));
          if (!p) return;
          nameIn.value = '';
          toast('ok', 'Project created', `${p.id} · ${p.name}`, 3000);
          await reload();
          await select(p.id);
        },
      }, 'create')));

    if (!projects.length) {
      railEl.appendChild(h('div.caveat.info', h('span',
        'No projects yet. A project groups runs by PATIENT so a therapy panel ' +
        'on one specimen can be compared against the same panel on another.')));
    }
    for (const p of projects) {
      const on = project && p.id === project.id;
      railEl.appendChild(h('div.proj-item', {
        class: on ? 'on' : null,
        tabindex: '0',
        onclick: () => select(p.id),
        onkeydown: (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); select(p.id); } },
      },
        h('div.proj-item-top',
          h('span.mono.strong', p.id),
          h('span.grow', p.name || NULL_GLYPH),
          h('button.btn.ghost.sm.danger', {
            title: 'Delete this project. Runs are NOT deleted.',
            'aria-label': `Delete project ${p.id}`,
            onclick: async (e) => {
              e.stopPropagation();
              const out = await guard(`DELETE /api/projects/${p.id}`,
                () => api.deleteProject(p.id));
              if (!out) return;
              toast('ok', `${p.id} deleted`,
                `${(out.runs_retained || []).length} run(s) retained in app/runs/ — a project owns ids, not results.`, 5000);
              if (project && project.id === p.id) project = null;
              await reload();
            },
          }, '×')),
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
        h('button.btn.ghost.sm.danger', {
          title: 'Remove this specimen. Its runs are NOT deleted.',
          'aria-label': `Remove specimen ${sp.id}`,
          onclick: async () => {
            const out = await guard('DELETE specimen',
              () => api.deleteSpecimen(project.id, sp.id));
            if (!out) return;
            toast('ok', `${sp.id} removed`,
              `${(out.runs_retained || []).length} run(s) retained in app/runs/.`, 4000);
            await select(project.id);
          },
        }, '×')),
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
        h('p', 'Create a project, add two sections as specimens, and attach the ' +
          'same drug panel to each. The comparison then answers the question ' +
          'this whole product is for: do two patients agree?')));
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

    /* 0. headline */
    const agrees = c.ranking && c.ranking.agrees;
    mainEl.appendChild(h('div.verdict',
      { class: agrees === true ? 'ok' : agrees === false ? 'crit' : 'info' },
      h('div.verdict-tag', agrees === true ? 'RANKING AGREES'
        : agrees === false ? 'RANKING DISAGREES' : 'NOT EVALUABLE'),
      h('div.verdict-body', c.headline)));

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

    /* 2. ranking — the primary result */
    const rk = c.ranking || {};
    mainEl.appendChild(h('div.hr-lbl', `2 · Ranking by ${METRIC_LABEL[c.metric] || c.metric} — THE PRIMARY RESULT`));
    if (!rk.drugs || rk.drugs.length < 2) {
      mainEl.appendChild(h('div.caveat', h('span', rk.verdict || 'not evaluable')));
    } else {
      mainEl.appendChild(h('div.tblwrap', h('table.tbl.zebra',
        h('thead', h('tr',
          h('th', 'specimen'),
          ...rk.drugs.map((d) => h('th.right', d)),
          h('th', 'order'))),
        h('tbody', rk.specimens.map((sid) => {
          const blk = rk.by_specimen[sid] || { order: [], values: {} };
          return h('tr',
            h('td', h('span.mono.strong', sid), ' ', label(sid)),
            ...rk.drugs.map((d) => h('td.right.num',
              fmtBy(metricFmt(c.metric), blk.values[d]),
              h('span.unit', METRIC_UNIT[c.metric] || ''))),
            h('td.mono.rank', blk.order.join(' > ')));
        })))));
      mainEl.appendChild(h('div.tblwrap', h('table.tbl',
        h('thead', h('tr', h('th', 'pair'), h('th', 'agrees'),
          h('th.right', 'Kendall tau-b'), h('th.right', 'discordant'),
          h('th', 'swapped drugs'))),
        h('tbody', (rk.pairwise || []).map((p) => h('tr',
          h('td.mono', `${p.a} vs ${p.b}`),
          h('td', p.agrees
            ? h('span.chip.st-done', 'yes') : h('span.chip.st-failed', 'no')),
          h('td.right.num', p.kendall_tau_b === null ? NULL_GLYPH : num(p.kendall_tau_b, 3)),
          h('td.right.num', int(p.discordant_pairs)),
          h('td.mono', (p.swapped || []).length
            ? p.swapped.map((s) => s.join(' ↔ ')).join('; ') : NULL_GLYPH)))))));
      mainEl.appendChild(h('div.caveat.info', h('span', rk.verdict)));
    }

    /* 3. absolute kill — flagged */
    mainEl.appendChild(h('div.hr-lbl',
      `3 · Absolute ${METRIC_LABEL[c.metric] || c.metric} — read ONLY where comparable`));
    if (!(c.kill || []).length) {
      mainEl.appendChild(h('div.caveat.info', h('span',
        'No drug is matched across two specimens at an identical RunSpec. Two ' +
        'runs are matched only when everything but the section and the drug is ' +
        'the same — same mesh, same modes, same PK, same window.')));
    } else {
      const ref = c.kill[0].reference;
      mainEl.appendChild(h('div.tblwrap', h('table.tbl.zebra',
        h('thead', h('tr',
          h('th', 'drug'),
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

    /* 4. composition */
    const comp = c.composition || { rows: [] };
    mainEl.appendChild(h('div.hr-lbl', '4 · Composition — the tissue difference'));
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
          h('th', 'drugs'), h('th.right', 'matched'))),
        h('tbody', c.matched_sets.map((ms) => h('tr',
          h('td.mono', { style: { fontSize: 'var(--fs-micro)' } },
            ms.is_primary ? h('b', 'PRIMARY · ') : null, ms.signature),
          h('td.mono', ms.specimens.join(', ')),
          h('td', ms.drugs.join(', ')),
          h('td.right.num', int(ms.n_matched_drugs))))))));
    }
  }

  render();

  return {
    el: dlg,
    button: btn,
    open,
    reload,
    get projectId() { return project ? project.id : null; },
  };
}
