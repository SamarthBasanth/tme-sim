// compound.js — BRING YOUR OWN COMPOUND: a SMILES in, a compiled compound the
// operator can add to a run, and — the part that matters — a blunt statement of
// what each compiled field actually is.
//
// The shipped library has twelve molecules. A researcher's own series is the
// reason they are here at all: they want to rank THEIR compounds on THEIR
// section. That is the one claim this project defends, so the feature is worth
// building. What is NOT worth building is a compiled compound that looks
// authoritative and is not.
//
// Three rules this panel is built to, in order:
//
//  1. NOTHING IS SHOWN WITHOUT ITS PROVENANCE. Every compiled field renders with
//     the same measured/computed/assumed chip the rest of the Studio uses (from
//     ./provenance.js — one vocabulary, one place). Nothing here is `measured`:
//     no instrument touched this compound. Two fields are `computed` (molecular
//     weight and TPSA are exact functions of the graph the operator typed);
//     everything else is `assumed`, because everything else is downstream of an
//     empirical fit. The headline above the table is the count — "8 of 11 fields
//     are predictions" — so the summary is visible before the detail.
//
//  2. THE REFUSAL IS PART OF THE PRODUCT. The backend refuses antibodies and
//     biologics, salts and mixtures, elements the Crippen logP model has no atom
//     type for, out-of-range molecular weights and unparseable input, each with
//     its own sentence and the name of the field at fault. Those come back as a
//     422 with `{code, message, field}` and are rendered on the field, not as a
//     stack trace and not as a toast that disappears.
//
//  3. WHAT MAY BE CONCLUDED IS ON SCREEN, NOT IN THE DOCS. The can/cannot lists
//     are the backend's own strings (tmesim/molecular/byo.py), rendered verbatim,
//     so the sentence in the module and the sentence on the screen cannot drift.
//
// Stub mode: compiling is REFUSED outright when js/api.js has fallen back to the
// mock. A fabricated compile would be the single most dangerous thing on this
// screen — a DrugParams that looks like chemistry and is a seeded PRNG.

import { h, clear } from '../util/dom.js';
import * as api from '../api.js';
import { notifyError, toast } from './toast.js';
import { confirmable } from './confirm.js';
import { chip } from './provenance.js';
import { caveatNode } from './inspect.js';
import { num, sci } from '../util/fmt.js';

/** Format a compiled value for display, choosing a scale that stays readable. */
function fmtVal(v, unit) {
  if (v === null || v === undefined || !Number.isFinite(v)) return '—';
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e6)) return sci(v, 3);
  if (unit === 'g/mol') return num(v, 3);
  if (a >= 100) return num(v, 1);
  return num(v, 4);
}

export function createCompoundPanel({ onChanged }) {
  /* ------------------------------------------------------------- form -- */

  const fieldEls = new Map();

  const nameIn = h('input.inp', {
    id: 'c-name', type: 'text', placeholder: 'ACME-101',
    maxlength: '64', autocomplete: 'off', spellcheck: 'false',
    oninput: () => clearFieldErrors(),
  });
  const smilesIn = h('textarea.inp.mono', {
    id: 'c-smiles', rows: '3', spellcheck: 'false', autocomplete: 'off',
    autocapitalize: 'off', autocorrect: 'off', wrap: 'soft',
    placeholder: 'CC(=O)Oc1ccccc1C(=O)O',
    style: { resize: 'vertical', minHeight: '52px', lineHeight: '1.4' },
    oninput: () => clearFieldErrors(),
    onkeydown: (e) => { if (e.key === 'Escape') { e.stopPropagation(); smilesIn.blur(); } },
  });
  const ic50In = h('input.inp.num', {
    id: 'c-ic50', type: 'number', step: 'any', min: '0', placeholder: 'optional',
    oninput: () => clearFieldErrors(),
  });
  const noteIn = h('input.inp', {
    id: 'c-note', type: 'text', maxlength: '200', autocomplete: 'off',
    placeholder: 'optional — series, assay, batch',
  });

  const lbl = (text, forId, title) => h('label.lbl', { for: forId, title: title || '' }, text);

  function row(label, forId, control, unitText, title) {
    const holder = unitText ? h('div.with-unit', control, h('span.unit', unitText)) : control;
    const errId = forId + '-err';
    const errEl = h('div.field-err.hidden', { id: errId, role: 'alert' });
    control.setAttribute('aria-describedby', errId);
    fieldEls.set(forId, { control, errEl });
    return [lbl(label, forId, title), holder, errEl];
  }

  const formEl = h('div.fgrid',
    row('Name', 'c-name', nameIn, null,
      'Labels the compound in the picker, in the run record and in any export. It may not shadow a curated library name.'),
    row('SMILES', 'c-smiles', smilesIn, null,
      'One connected small molecule. Salts, mixtures, peptides and antibodies are refused — see the panel below for why.'),
    row('IC50', 'c-ic50', ic50In, 'uM',
      'YOURS if you have one. Structure cannot give potency. Leave blank and a flagged 25 uM placeholder is used instead — and every absolute potency number becomes meaningless.'),
    row('Note', 'c-note', noteIn, null, 'Free text kept with the compound. Not used by the model.'),
  );

  /* --------------------------------------------------------- controls -- */

  const compileBtn = h('button.btn', {
    onclick: () => run(false), title: 'POST /api/compounds/compile — compiles and shows provenance, stores nothing',
  }, 'Compile');
  const addBtn = h('button.btn.primary', {
    disabled: true, onclick: () => run(true), title: 'POST /api/compounds — stores it and adds it to the drug picker',
  }, 'Add to picker');
  const statusChip = h('span.chip.st-queued', 'not compiled');

  const controlsEl = h('div.frow', { style: { padding: 'var(--sp-4) var(--sp-5)' } },
    compileBtn, addBtn, h('span.grow'), statusChip);

  /* ---------------------------------------------------------- output -- */

  const outEl = h('div.stack.tight');
  const storedEl = h('div.stack.tight');

  /* ------------------------------------------------------------ root -- */

  const root = h('details.sect', { open: false, id: 'sect-compound' },
    h('summary', 'Your compound', h('span.hint', 'SMILES → DrugParams')),
    h('div.sect-body',
      h('p.dim', {
        style: {
          margin: '0 0 var(--sp-4)', fontSize: 'var(--fs-xs)',
          lineHeight: 'var(--lh-body)',
        },
      },
      'Compile one of your own molecules and compare it against the others on '
      + 'this section. Everything the compiler derives from a structure is a '
      + 'prediction, and the table below says which — read it before you use the '
      + 'result. In particular: diffusivity is a function of molecular weight '
      + 'alone, and that has been measured against 71 compounds and refuted, so '
      + 'a series of similar-mass analogues is the one case this model is known '
      + 'not to resolve.'),
      formEl),
    controlsEl,
    h('div.sect-body', { style: { padding: '0' } }, outEl),
    h('div.sect-body', { style: { paddingTop: '0' } }, storedEl),
  );

  /* -------------------------------------------------------- behaviour -- */

  function clearFieldErrors() {
    for (const [, { control, errEl }] of fieldEls) {
      control.removeAttribute('aria-invalid');
      control.classList.remove('invalid');
      errEl.classList.add('hidden');
      errEl.textContent = '';
    }
    markStale();
  }

  /**
   * The rendered table describes the inputs it was compiled from. The moment any
   * of them changes it stops describing what is on screen — and a provenance
   * table that belongs to a different molecule is worse than no table, because
   * it is the one thing in this panel the operator is meant to trust. So the
   * output is visibly retired and "Add to picker" is withdrawn until a fresh
   * compile has run.
   */
  function markStale() {
    if (!outEl.firstChild) return;
    addBtn.disabled = true;
    if (statusChip.textContent !== 'stale — recompile') {
      statusChip.className = 'chip st-queued';
      statusChip.textContent = 'stale — recompile';
      outEl.style.opacity = '0.45';
    }
  }

  const FIELD_ID = { name: 'c-name', smiles: 'c-smiles', ic50_uM: 'c-ic50' };

  function showRefusal(detail) {
    clearFieldErrors();
    const code = (detail && detail.code) || 'refused';
    const message = (detail && detail.message)
      || 'The backend refused this compound but sent no reason. That is a bug — do not work around it.';
    const id = FIELD_ID[(detail && detail.field) || 'smiles'];
    const t = fieldEls.get(id);
    if (t) {
      t.control.setAttribute('aria-invalid', 'true');
      t.control.classList.add('invalid');
      t.errEl.textContent = message;
      t.errEl.classList.remove('hidden');
    }
    statusChip.className = 'chip st-failed';
    statusChip.textContent = 'refused';
    addBtn.disabled = true;

    clear(outEl);
    outEl.style.opacity = '';
    outEl.appendChild(h('div.caveats',
      caveatNode(`REFUSED (${code}): ${message}`, 1, 'crit')));
    outEl.appendChild(h('div.dim', {
      style: {
        padding: 'var(--sp-3) var(--sp-5) var(--sp-5)',
        fontSize: 'var(--fs-micro)', lineHeight: 'var(--lh-body)',
      },
    },
    'This is a refusal, not a failure. The compiler stops on input it cannot '
    + 'honestly compile — an antibody or biologic, a salt or mixture, an element '
    + 'the Crippen logP model has no atom type for, a molecular weight outside '
    + 'the 50–2000 g/mol band this bridge asserts, or an unparseable SMILES — '
    + 'because the alternative is a confident-looking number derived from '
    + 'nothing.'));
  }

  /* --------------------------------------------------- the provenance -- */

  function provRow(f) {
    const isPred = !!f.is_prediction;
    return h('div', {
      style: {
        display: 'grid',
        gridTemplateColumns: 'minmax(0, 1fr) auto auto',
        columnGap: 'var(--sp-4)', rowGap: 'var(--sp-1)',
        alignItems: 'baseline',
        padding: 'var(--sp-3) var(--sp-5)',
        borderBottom: '1px solid var(--line-1)',
        borderLeft: isPred ? '2px solid var(--warn)' : '2px solid transparent',
        fontSize: 'var(--fs-xs)',
      },
    },
    h('span', { style: { color: 'var(--fg-0)' } }, f.label),
    h('span.mono', { style: { color: 'var(--fg-1)', whiteSpace: 'nowrap' } },
      fmtVal(f.value, f.unit) + (f.unit ? ' ' + f.unit : '')),
    chip(f.kind, f.method),
    h('div', {
      style: {
        gridColumn: '1 / -1', color: 'var(--fg-2)',
        fontSize: 'var(--fs-micro)', lineHeight: 'var(--lh-body)',
        overflowWrap: 'anywhere',
      },
    }, f.method),
    );
  }

  function bulletList(title, items, tone) {
    return h('div', { style: { padding: 'var(--sp-4) var(--sp-5)' } },
      h('div.hr-lbl', { style: { margin: '0 0 var(--sp-2)' } }, title),
      h('ul', {
        style: {
          margin: '0', paddingLeft: '1.1em', display: 'flex',
          flexDirection: 'column', gap: 'var(--sp-2)',
          fontSize: 'var(--fs-xs)', lineHeight: 'var(--lh-body)',
          color: tone === 'crit' ? 'var(--fg-0)' : 'var(--fg-1)',
        },
      }, items.map((s) => h('li', s))),
    );
  }

  function render(doc) {
    clear(outEl);
    outEl.style.opacity = '';
    const fields = doc.fields || [];
    const nPred = doc.predicted_field_count ?? fields.filter((f) => f.is_prediction).length;

    // The headline. A compiled compound's honest one-line summary is the ratio
    // of predictions to fields, and it goes above the numbers, not below them.
    outEl.appendChild(h('div', {
      style: {
        padding: 'var(--sp-4) var(--sp-5)',
        background: 'var(--warn-bg)',
        borderTop: '1px solid var(--line-2)',
        borderBottom: '1px solid var(--line-2)',
        fontSize: 'var(--fs-xs)', lineHeight: 'var(--lh-body)',
        color: 'var(--fg-0)',
      },
    },
    h('b', `${nPred} of ${fields.length} compiled fields are PREDICTIONS`),
    ' from the structure you typed, not measurements. Only molecular weight and '
    + 'TPSA are exact functions of the graph'
    + (doc.ic50_supplied ? ', and the IC50 is the one you supplied.' : '.'),
    ));

    outEl.appendChild(h('div.hr-lbl',
      { style: { padding: 'var(--sp-4) var(--sp-5) var(--sp-1)' } },
      'Compiled fields'));
    const table = h('div');
    for (const f of fields) table.appendChild(provRow(f));
    outEl.appendChild(table);

    // Caveats, rendered with the same component (and the same hard-stop
    // promotion rule) the result panel uses. A hard stop must look the same
    // wherever it appears.
    const cav = [...(doc.caveats || []), ...(doc.warnings || [])];
    if (cav.length) {
      outEl.appendChild(h('div.hr-lbl',
        { style: { padding: 'var(--sp-5) var(--sp-5) var(--sp-1)' } },
        `Caveats · ${cav.length}`));
      const box = h('div.caveats');
      cav.forEach((c, i) => box.appendChild(caveatNode(String(c), i + 1)));
      outEl.appendChild(box);
    }

    if ((doc.can_conclude || []).length) {
      outEl.appendChild(bulletList('What you CAN conclude from this compound', doc.can_conclude));
    }
    if ((doc.cannot_conclude || []).length) {
      outEl.appendChild(bulletList('What you CANNOT conclude', doc.cannot_conclude, 'crit'));
    }
  }

  /* ------------------------------------------------------------ calls -- */

  async function run(persist) {
    if (api.api && api.api.stub) {
      showRefusal({
        code: 'no_backend',
        field: 'smiles',
        message: 'The Studio is running on stub data because the backend is not '
          + 'answering, and a compound will not be compiled against it. A '
          + 'fabricated DrugParams would look exactly like a real one. Start the '
          + 'backend and try again.',
      });
      return;
    }
    const body = {
      name: nameIn.value.trim(),
      smiles: smilesIn.value.trim(),
      ic50_uM: ic50In.value === '' ? null : Number(ic50In.value),
    };
    if (persist && noteIn.value.trim()) body.note = noteIn.value.trim();

    compileBtn.disabled = true;
    addBtn.disabled = true;
    statusChip.className = 'chip st-running';
    statusChip.textContent = persist ? 'adding' : 'compiling';
    try {
      const doc = persist
        ? await api.postCompound(body)
        : await api.postCompoundCompile(body);
      clearFieldErrors();
      render(doc);
      statusChip.className = 'chip st-done';
      statusChip.textContent = persist ? 'added' : 'compiled';
      addBtn.disabled = persist;
      if (persist) {
        toast('info', `Compound ${doc.name} added`,
          'It is in the drug picker now. Its absolute kill, AUC and penetration '
          + 'are still not defensible outputs — rank it, do not quote it.', 9000);
        await refreshStored();
        if (onChanged) onChanged(doc);
      }
    } catch (err) {
      const detail = err && err.payload && err.payload.detail;
      if (detail && typeof detail === 'object' && detail.message) {
        showRefusal(detail);
      } else {
        statusChip.className = 'chip st-failed';
        statusChip.textContent = 'failed';
        notifyError(persist ? 'POST /api/compounds failed'
          : 'POST /api/compounds/compile failed', err);
      }
    } finally {
      compileBtn.disabled = false;
    }
  }

  async function refreshStored() {
    let list = [];
    try {
      list = await api.listCompounds();
    } catch {
      clear(storedEl);
      return;
    }
    clear(storedEl);
    if (!list.length) return;
    storedEl.appendChild(h('div.hr-lbl', { style: { margin: '0 0 var(--sp-2)' } },
      `Your compounds · ${list.length}`));
    for (const rec of list) {
      storedEl.appendChild(h('div.frow', {
        style: {
          gap: 'var(--sp-3)', alignItems: 'baseline',
          padding: 'var(--sp-2) 0', borderBottom: '1px solid var(--line-1)',
          fontSize: 'var(--fs-xs)',
        },
      },
      h('span', { style: { color: 'var(--fg-0)' } }, rec.name),
      h('span.chip.assumed', {
        title: rec.ic50_uM === null || rec.ic50_uM === undefined
          ? 'No IC50 was supplied for this compound; a flagged 25 uM placeholder stands in.'
          : `IC50 ${rec.ic50_uM} uM, supplied by you. Not measured here.`,
      }, rec.ic50_uM === null || rec.ic50_uM === undefined ? 'no IC50' : `IC50 ${rec.ic50_uM}`),
      h('span.grow'),
      h('span.mono.dim', {
        title: rec.smiles,
        style: { maxWidth: '14ch', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' },
      }, rec.smiles || ''),
      // Two-step, and then undoable. This row holds a SMILES the operator
      // typed and that may exist nowhere else — an unpublished series member
      // is not a cache entry. `confirmable` makes the first click a question;
      // the undo below makes even the answered question reversible, because
      // everything needed to rebuild the record (name, structure, potency,
      // note) is in `rec` and POST /api/compounds is idempotent in effect.
      confirmable(h('button.btn.ghost.sm', {
        title: 'Delete this compound. DELETE /api/compounds/' + rec.id,
        'aria-label': `Remove compound ${rec.name}`,
      }, 'remove'), {
        confirmLabel: `delete ${rec.name}?`,
        armedTitle: `This removes ${rec.name} and its SMILES from the store and `
          + 'from the drug picker. Press again to confirm; Escape cancels.',
        onConfirm: async () => {
          try {
            await api.deleteCompound(rec.id);
            await refreshStored();
            if (onChanged) onChanged(null);
            const t = toast('warn', `${rec.name} removed`,
              'Its SMILES is no longer in the store or the drug picker. '
              + 'Runs already made with it are untouched.', 12000);
            t.querySelector('.toast-body').appendChild(
              h('button.btn.sm', {
                style: { marginTop: 'var(--sp-3)' },
                title: 'Re-compile and re-store this compound from the record just deleted',
                onclick: async () => {
                  t.remove();
                  try {
                    await api.postCompound({
                      name: rec.name,
                      smiles: rec.smiles,
                      ic50_uM: rec.ic50_uM ?? null,
                      ...(rec.note ? { note: rec.note } : {}),
                    });
                    await refreshStored();
                    if (onChanged) onChanged(null);
                    toast('ok', `${rec.name} restored`, '', 4000);
                  } catch (err) { notifyError('Undo failed — POST /api/compounds', err); }
                },
              }, 'undo'));
          } catch (err) { notifyError('DELETE /api/compounds failed', err); }
        },
      }),
      ));
    }
  }

  return {
    el: root,
    refresh: refreshStored,
    open() { root.open = true; nameIn.focus(); },
  };
}
