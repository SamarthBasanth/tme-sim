// account.js — the Account panel: one place for the questions that are about
// THIS DEPLOYMENT rather than about tissue.
//
// WHY IT EXISTS
// -------------
// Access to this Studio was, until now, a single modal that appeared once, took
// a key, and vanished. Everything after that was invisible: which key this
// browser is sending, whether it is remembered here or only for this tab, how
// to call the same API from a script, what this server will refuse before you
// spend a run finding out, and whether anyone is being charged. A user who had
// two keys, or who wanted to sign out on a borrowed laptop, or who wanted to
// automate a sweep, had no surface at all. See app/server/account.py.
//
// FIVE THINGS IT WILL NOT DO
// --------------------------
// 1. It never renders the key. A fingerprint identifies it; the secret itself
//    is write-only from the UI's point of view.
// 2. The copyable examples carry `$VIVUM_API_KEY`, never the live key. A
//    snippet with a real secret baked in is how secrets end up pasted into
//    chat, and "copy" is the button people press without reading.
// 3. It shows no usage number the server did not send, and repeats the
//    server's own word for what those numbers are: an inventory of files that
//    exist right now, not a meter. Deleting a run decrements them.
// 4. It never invents a plan, a price, a quota or an invoice. Billing is not
//    enabled; the section says so, says what is blocking it, and shows nothing
//    that could be mistaken for a charge.
// 5. In stub mode it refuses to render at all rather than describe a backend
//    that is not answering. `getAccount` / `getLimits` are `noFallback` in
//    api.js for the same reason.

import { h, clear, $ } from '../util/dom.js';
import { int, dur, NULL_GLYPH } from '../util/fmt.js';
import * as api from '../api.js';
import { setApiKey, apiKeyStorage } from '../api.js';
import { toast } from './toast.js';

const TABS = [
  ['access', 'Access'],
  ['api', 'API'],
  ['limits', 'Limits'],
  ['usage', 'Usage'],
  ['billing', 'Plan & billing'],
];

/** Where this browser is keeping the key, in the user's terms. */
const STORAGE_WORDS = {
  local: ['this browser', 'Stored on this browser. It survives closing the tab '
    + 'and reopening Studio, and stays until you sign out.'],
  session: ['this tab only', 'Stored for this tab only. Closing the tab '
    + 'forgets it; a reload keeps it.'],
  none: ['not stored', 'No key is set in this browser.'],
  unavailable: ['not stored — storage unavailable', 'The key works for this '
    + 'page load but could not be written anywhere: a private window, or a '
    + 'browser set to block site data. Reloading will ask for it again.'],
};

function row(label, value, title) {
  return h('tr', { title: title || '' },
    h('td.lbl', label),
    h('td.right.mono', value === null || value === undefined
      ? NULL_GLYPH : String(value)));
}

/** A copy button whose payload is fixed text. Mirrors ui/reach.js::copyBtn. */
function copyBtn(getText, what) {
  const b = h('button.btn.ghost.tiny', {
    type: 'button', title: `Copy ${what}`,
  }, 'copy');
  let restore = 0;
  b.addEventListener('click', () => {
    const say = (m) => {
      b.textContent = m;
      clearTimeout(restore);
      restore = setTimeout(() => { b.textContent = 'copy'; }, 2000);
    };
    const cb = navigator.clipboard;
    if (!cb || typeof cb.writeText !== 'function') { say('no clipboard'); return; }
    cb.writeText(getText()).then(() => say('copied'), () => say('failed'));
  });
  return b;
}

function codeBlock(text, what) {
  return h('div.acct-code',
    h('pre', h('code', text)),
    h('div.acct-code-act', copyBtn(() => text, what)));
}

/** The copy-paste examples, as pure text.
 *
 * EXPORTED, and a pure function of (base, header) ONLY, so the invariant that
 * matters here can actually be tested rather than reviewed: the live key is
 * not an input to it, so no example it produces can contain the key. The
 * placeholder is a shell variable the reader exports themselves.
 *
 * The invariant is worth a function of its own because "copy" is the button
 * people press without reading, and a snippet with a real secret baked into it
 * is how secrets end up pasted into a chat window.
 */
export function apiExamples(base, header) {
  const b = String(base || '').replace(/\/+$/, '');
  const hd = String(header || 'X-API-Key');
  const auth = `-H "${hd}: $VIVUM_API_KEY"`;
  return {
    placeholder: '$VIVUM_API_KEY',
    health: `curl -s ${b}/api/health ${auth}`,
    catalogue:
      `curl -s ${b}/api/sections ${auth}\n`
      + `curl -s ${b}/api/drugs   ${auth}\n`
      + `curl -s ${b}/api/limits  ${auth}`,
    validate:
      `curl -s ${b}/api/validate \\\n`
      + `  ${auth} -H 'content-type: application/json' \\\n`
      + `  -d "$(curl -s ${b}/api/defaults ${auth})"`,
    submit:
      `RUN=$(curl -s ${b}/api/runs \\\n`
      + `  ${auth} -H 'content-type: application/json' \\\n`
      + `  -d @spec.json \\\n`
      + `  | python3 -c 'import sys,json;print(json.load(sys.stdin)["run_id"])')\n`
      + `curl -s ${b}/api/runs/$RUN/status ${auth}`,
  };
}

export function createAccountPanel() {
  let acct = null;      // GET /api/account body
  let limits = null;    // GET /api/limits body
  let err = null;
  let tab = 'access';
  let busy = false;

  const navEl = h('div.acct-nav');
  const mainEl = h('div.acct-main');

  const dlg = h('dialog#acct-dlg',
    h('div.panel-head',
      h('h2', 'Account'),
      h('span.hint', 'access, API and billing for this deployment'),
      h('span.grow'),
      h('button.btn.ghost.sm', {
        onclick: () => reload(),
        title: 'GET /api/account and GET /api/limits',
      }, 'refresh'),
      h('button.btn.ghost.sm', { onclick: () => dlg.close() }, 'close')),
    h('div.acct-body', navEl, h('div.acct-div'), mainEl));

  document.body.appendChild(dlg);

  const btn = h('button.btn.ghost.sm', {
    id: 'btn-account',
    title: 'Account: which key this browser sends, how to call the API, '
      + 'what this server refuses, and billing',
    onclick: open,
  }, 'account');
  const slot = $('.topbar-group.last');
  if (slot) slot.insertBefore(btn, slot.firstChild);

  /* ------------------------------------------------------------ data --- */

  function offline() { return !!(api.api && api.api.stub); }

  async function open() {
    if (!dlg.open) dlg.showModal();
    await reload();
  }

  async function reload() {
    if (busy) return;
    busy = true;
    err = null;
    render();
    try {
      // Sequential, not Promise.all: if the first one 401s there is no point
      // asking the second, and a single clear failure beats two.
      acct = await api.getAccount();
      limits = await api.getLimits();
    } catch (e) {
      acct = null; limits = null;
      err = e;
    } finally {
      busy = false;
      render();
    }
  }

  /* ---------------------------------------------------------- render --- */

  function render() {
    clear(navEl);
    for (const [id, label] of TABS) {
      navEl.appendChild(h('button.acct-tab', {
        type: 'button',
        'aria-pressed': String(id === tab),
        onclick: () => { tab = id; render(); },
      }, label));
    }

    clear(mainEl);
    if (offline()) {
      mainEl.appendChild(h('div.caveat.crit', h('span',
        h('b', 'THE BACKEND IS NOT ANSWERING. '),
        'Studio is running on synthetic data, and nothing on this panel can be '
        + 'answered from it: which key you are sending, what this server '
        + 'refuses and whether anyone is charged are facts about a real '
        + 'server. Rather than show you plausible values for them, this panel '
        + 'shows nothing. Reconnect and reopen it.')));
      return;
    }
    if (busy && !acct) { mainEl.appendChild(h('div.hint', 'loading…')); return; }
    if (err) {
      mainEl.appendChild(h('div.caveat.crit', h('span',
        h('b', 'COULD NOT READ THIS DEPLOYMENT. '),
        String((err && err.message) || err),
        err && err.status === 401
          ? ' — the key this browser is sending was rejected. Replace it below.'
          : '')));
      mainEl.appendChild(accessTab());
      return;
    }
    if (!acct) { mainEl.appendChild(h('div.hint', 'no data')); return; }

    if (tab === 'access') mainEl.appendChild(accessTab());
    else if (tab === 'api') mainEl.appendChild(apiTab());
    else if (tab === 'limits') mainEl.appendChild(limitsTab());
    else if (tab === 'usage') mainEl.appendChild(usageTab());
    else mainEl.appendChild(billingTab());
  }

  /* -- access ------------------------------------------------------------ */

  function accessTab() {
    const a = (acct && acct.auth) || {};
    const store = apiKeyStorage();
    const [where, whyWhere] = STORAGE_WORDS[store] || STORAGE_WORDS.unavailable;

    const wrap = h('div.acct-sec');
    wrap.appendChild(h('div.hr-lbl', 'This browser'));
    wrap.appendChild(h('div.tblwrap', h('table.tbl',
      row('key fingerprint', a.caller_fingerprint || NULL_GLYPH,
        'A short one-way digest of the key this browser sends. It is an '
        + 'identifier, not the key: nothing authenticates against it, and it '
        + 'cannot be turned back into the secret.'),
      row('stored', where, whyWhere),
      row('header sent', a.header || 'X-API-Key'),
      row('accepted by server', a.matches === true ? 'yes' : 'no'))));
    wrap.appendChild(h('div.hint', whyWhere));

    wrap.appendChild(h('div.hr-lbl', 'This deployment'));
    wrap.appendChild(h('div.tblwrap', h('table.tbl',
      row('key required', a.required === true ? 'yes' : 'no'),
      row('scheme', a.scheme || NULL_GLYPH),
      row('key set by', a.source || NULL_GLYPH,
        'Where the server got its key: the environment, an --api-key flag, or '
        + 'generated at startup and printed once.'),
      row('server fingerprint', a.key_fingerprint || NULL_GLYPH))));

    // The honest paragraph. Everything a user assumes about an "account" —
    // sign-up, a second key, rotation from the UI — is absent here, and the
    // panel says so rather than leaving them to discover it.
    wrap.appendChild(h('div.caveat.info', h('span',
      h('b', 'THERE IS NO SIGN-UP AND NO SECOND KEY. '),
      'This deployment authenticates with one shared secret. There is no user '
      + 'account behind it, no self-service registration, and no endpoint that '
      + 'issues or rotates a key — issuing keys needs a key store this '
      + 'deployment does not have. If someone shared this Studio with you, the '
      + 'key comes from them. ',
      a.how_to_rotate || '')));

    wrap.appendChild(h('div.hr-lbl', 'Replace or sign out'));
    const input = h('input.inp', {
      type: 'password', autocomplete: 'off', spellcheck: 'false',
      placeholder: 'new X-API-Key value', style: { width: '100%' },
    });
    const remember = h('input', {
      type: 'checkbox', checked: store !== 'session',
    });
    const applyBtn = h('button.btn.primary.sm', {
      onclick: async () => {
        const k = input.value.trim();
        if (!k) { input.focus(); return; }
        const prev = api.api.apiKey;
        setApiKey(k, remember.checked);
        const p = await api.probe();
        if (p.authRequired) {
          setApiKey(prev, store !== 'session');
          toast('err', 'That key was rejected',
            'The server refused it, so the previous key has been kept.');
          input.select();
          return;
        }
        input.value = '';
        toast('ok', 'Key replaced', `Now stored on ${remember.checked
          ? 'this browser' : 'this tab only'}.`);
        await reload();
      },
    }, 'Use this key');
    const outBtn = h('button.btn.ghost.sm', {
      title: 'Forget the key in this browser. Studio will ask for one again.',
      onclick: () => {
        setApiKey('');
        toast('ok', 'Signed out',
          'The key was removed from this browser. Reload to enter another.');
        render();
      },
    }, 'Sign out of this browser');

    wrap.appendChild(h('div.acct-form',
      input,
      h('label.acct-check', remember, 'Stay signed in on this browser'),
      h('div.btn-row', applyBtn, outBtn)));
    return wrap;
  }

  /* -- api --------------------------------------------------------------- */

  function apiTab() {
    const base = (api.api && api.api.base) || window.location.origin;
    const a = (acct && acct.auth) || {};
    const hdr = a.header || 'X-API-Key';
    const ex = apiExamples(base, hdr);
    const wrap = h('div.acct-sec');

    wrap.appendChild(h('div.hr-lbl', 'Calling this deployment from a script'));
    wrap.appendChild(h('div.tblwrap', h('table.tbl',
      row('base URL', base),
      row('auth header', hdr),
      row('contract', (acct && acct.contract) || NULL_GLYPH))));

    // The key is NOT interpolated into any of these. A copy button people
    // press without reading must not put a live secret on the clipboard.
    wrap.appendChild(h('div.caveat.info', h('span',
      h('b', 'THE EXAMPLES BELOW CARRY $VIVUM_API_KEY, NOT YOUR KEY. '),
      'Export it in your shell first (', h('code', 'export VIVUM_API_KEY=…'),
      '). Nothing on this page ever puts the secret itself on your clipboard.')));

    wrap.appendChild(h('div.hr-lbl', 'Check it works'));
    wrap.appendChild(codeBlock(ex.health, 'the health check'));

    wrap.appendChild(h('div.hr-lbl', 'List what this server can run on'));
    wrap.appendChild(codeBlock(ex.catalogue, 'the catalogue calls'));

    wrap.appendChild(h('div.hr-lbl', 'Validate a spec before spending a run'));
    wrap.appendChild(h('div.hint',
      'POST /api/validate answers with the same errors and warnings the run '
      + 'endpoint would refuse with, and costs nothing. It is the call to make '
      + 'first.'));
    wrap.appendChild(codeBlock(ex.validate, 'the validate call'));

    wrap.appendChild(h('div.hr-lbl', 'Submit and poll'));
    wrap.appendChild(codeBlock(ex.submit, 'the submit and poll calls'));

    wrap.appendChild(h('div.caveat', h('span',
      'Runs are asynchronous: POST /api/runs answers 202 with a run id, and '
      + 'the work happens on a single worker shared with the Studio session '
      + 'you are looking at. Two clients submitting at once queue behind each '
      + 'other. The full route list is in API_CONTRACT.md in the repository.')));
    return wrap;
  }

  /* -- limits ------------------------------------------------------------ */

  function limitsTab() {
    const L = limits || {};
    const wrap = h('div.acct-sec');
    wrap.appendChild(h('div.hint',
      'What this server refuses, and where the line between "slow" and '
      + '"refused" sits. The run spec editor reads the same numbers and shows '
      + 'them against your spec before you submit.'));
    wrap.appendChild(h('div.tblwrap', h('table.tbl',
      row('grid ceiling (refused above)', L.max_grid_cells === undefined
        ? NULL_GLYPH : int(L.max_grid_cells)),
      row('ceiling set by', L.max_grid_cells_source || NULL_GLYPH,
        '"default" is the 40 M "no request should ever ask for this" number. '
        + '"TMESIM_MAX_GRID_CELLS" means the operator lowered it to what this '
        + 'machine can actually afford, which a small instance should do.'),
      row('default ceiling', L.max_grid_cells_default === undefined
        ? NULL_GLYPH : int(L.max_grid_cells_default)),
      row('large grid (allowed, slow)', L.large_grid_cells === undefined
        ? NULL_GLYPH : int(L.large_grid_cells)),
      row('minimum dx', L.min_dx_um === undefined
        ? NULL_GLYPH : `${L.min_dx_um} um`),
      row('solver timestep', L.pde_dt_s === undefined
        ? NULL_GLYPH : `${L.pde_dt_s} s`))));
    // The exposure block, ABOVE the formula, because an operator whose server
    // is one request from dying does not need to know how the count is
    // computed first. Rendered only when the server actually sent it, so an
    // older backend degrades to the rest of the tab rather than to a blank.
    const ox = L.oom_exposure;
    if (ox) {
      const m = ox.measured || {};
      const memMB = Number.isFinite(ox.host_memory_bytes)
        ? `${int(Math.round(ox.host_memory_bytes / 1048576))} MB`
        : 'not declared to this process';
      wrap.appendChild(h('div.hr-lbl', 'Can one request kill this server?'));
      wrap.appendChild(h('div.tblwrap', h('table.tbl',
        row('host memory', memMB, ox.host_memory_source || ''),
        row('ceiling permits the measured kill',
          ox.ceiling_permits_measured_oom ? 'YES' : 'no'),
        row('host at or below the measured instance',
          ox.host_at_or_below_measured_instance === null
            ? 'unknown' : (ox.host_at_or_below_measured_instance ? 'yes' : 'no'),
          'Unknown means no limit is declared to this process — which is not '
          + 'the same as having room. A platform can cap memory without '
          + 'exposing it.'))));

      if (ox.ceiling_permits_measured_oom
          && ox.host_at_or_below_measured_instance) {
        wrap.appendChild(h('div.caveat.crit', h('span',
          h('b', 'ONE ORDINARY REQUEST CAN TAKE THIS SERVICE DOWN. '),
          `This host declares ${memMB}, at or below the ${int(m.instance_mb)} `
          + `MB instance where ${int(m.oom_voxels)} voxels OOM-KILLED the `
          + 'container — and the grid ceiling still permits that spec. It is '
          + 'not the run that fails: the PROCESS dies and /healthz stops '
          + 'answering, for everyone. Set TMESIM_MAX_GRID_CELLS on the host.')));
      } else if (ox.ceiling_permits_measured_oom && ox.ceiling_is_default) {
        wrap.appendChild(h('div.caveat', h('span',
          h('b', 'THE CEILING IS THE UNCONFIGURED DEFAULT. '),
          `It permits the ${int(m.oom_voxels)}-voxel spec that OOM-killed a `
          + `${int(m.instance_mb)} MB container, and no memory limit is `
          + 'declared here. On a small host, set TMESIM_MAX_GRID_CELLS.')));
      }

      // Never let the two numbers above be read as a model.
      wrap.appendChild(h('div.hint', m.note || ''));
    }

    wrap.appendChild(h('div.hr-lbl', 'How the grid count is computed'));
    wrap.appendChild(codeBlock(L.grid_cells_formula || '(not published)',
      'the grid formula'));
    wrap.appendChild(h('div.hint',
      `Needs ${L.grid_cells_requires || 'the section extent'}; a section that `
      + 'does not record one cannot be checked before submission at all, and '
      + 'the editor says so rather than assuming the spec is fine.'));
    return wrap;
  }

  /* -- usage ------------------------------------------------------------- */

  function usageTab() {
    const u = (acct && acct.usage) || {};
    const d = (acct && acct.deployment) || {};
    const by = u.runs_by_status || {};
    const wrap = h('div.acct-sec');

    // The server's own framing, first and in full. A page headed "Usage" is
    // read as a meter unless it says otherwise before the numbers.
    // The server's own sentence, READ, not re-authored. This block used to
    // restate `usage.kind` in the client's own words while the module header
    // above claimed it "repeats the server's own word for what those numbers
    // are". It did not, and a caveat written in two places is corrected in
    // one. `kind` leads; the elaboration stays here because the server sends a
    // label, not a paragraph.
    wrap.appendChild(h('div.caveat.info', h('span',
      h('b', `THIS IS ${String(u.kind || 'AN INVENTORY, NOT A METER')
        .toUpperCase()}. `),
      'These are the things that exist on this server right now. Deleting a '
      + 'run makes the number go down. Nothing here accumulates, nothing is '
      + 'billed from it, and it is not a record of what has been done.')));

    wrap.appendChild(h('div.tblwrap', h('table.tbl',
      row('runs stored', int(u.runs_total || 0)),
      ...Object.keys(by).sort().map((k) => row(`  · ${k}`, int(by[k]))),
      row('projects', u.projects === undefined ? NULL_GLYPH : int(u.projects)),
      row('reports', u.reports === null || u.reports === undefined
        ? NULL_GLYPH : int(u.reports)),
      row('process uptime', Number.isFinite(d.uptime_s)
        ? dur(d.uptime_s) : NULL_GLYPH),
      // Served by account.py and, until now, silently dropped here. A field
      // the server computes and the client discards is a field nobody will
      // notice going wrong.
      row('reconstruction cache',
        (u.reconstruction_cache && Number.isFinite(u.reconstruction_cache.size))
          ? `${int(u.reconstruction_cache.size)} / `
            + `${int(u.reconstruction_cache.maxsize)} entries · `
            + `${int(u.reconstruction_cache.hits)} hits, `
            + `${int(u.reconstruction_cache.misses)} misses`
          : NULL_GLYPH,
        'Reconstructions held in memory. Warm entries are why a second run on '
        + 'the same section and dx is fast; it is not a usage record.'))));

    if (u.attribution) {
      wrap.appendChild(h('div.caveat.crit', h('span',
        h('b', 'NO RUN CAN BE ATTRIBUTED TO A PERSON. '), u.attribution,
        '. Anyone holding the key is indistinguishable from anyone else '
        + 'holding it, in the logs and here.')));
    }
    return wrap;
  }

  /* -- billing ----------------------------------------------------------- */

  function billingTab() {
    const b = (acct && acct.billing) || {};
    const wrap = h('div.acct-sec');

    // Loud, first, and unhedged. Every other element on this tab is an empty
    // slot for an integration that does not exist, and an empty slot is only
    // honest while the sentence above it is unmissable.
    wrap.appendChild(h('div.caveat.crit', h('span',
      h('b', 'BILLING IS NOT ENABLED ON THIS DEPLOYMENT. '),
      b.note || 'No charge is made, no payment method is stored, and no usage '
        + 'is metered.')));

    wrap.appendChild(h('div.tblwrap', h('table.tbl',
      row('status', b.status || NULL_GLYPH),
      row('provider', b.provider || NULL_GLYPH),
      row('plan', b.plan || NULL_GLYPH),
      row('seats', b.seats === null || b.seats === undefined
        ? NULL_GLYPH : int(b.seats)),
      row('billing period', b.period || NULL_GLYPH),
      row('metered quantities', (b.metered_quantities || []).length
        ? b.metered_quantities.join(', ') : NULL_GLYPH),
      row('invoices', (b.invoices || []).length
        ? int(b.invoices.length) : NULL_GLYPH))));

    if (b.status === 'configured-but-not-implemented') {
      wrap.appendChild(h('div.caveat.crit', h('span',
        h('b', 'A PROVIDER IS NAMED BUT NOTHING IS WIRED TO IT. '),
        'TMESIM_BILLING_PROVIDER is set on this host, and no implementation is '
        + 'mounted. Nothing is charged. This state is reported rather than '
        + 'hidden so it cannot be mistaken for a working integration.')));
    }

    if ((b.blocked_on || []).length) {
      wrap.appendChild(h('div.hr-lbl', 'What has to exist before usage can be billed'));
      wrap.appendChild(h('ul.acct-list',
        ...b.blocked_on.map((r) => h('li', r))));
      wrap.appendChild(h('div.hint',
        'These are preconditions, not a schedule. Neither is a small change: '
        + 'the first needs per-caller identity in place of one shared secret, '
        + 'and the second needs a usage record that does not go down when a '
        + 'file is deleted.'));
    }
    return wrap;
  }

  return { node: dlg, open, reload, get element() { return dlg; } };
}
