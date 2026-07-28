// state.js — a 40-line observable store. No framework.

export function createStore(initial) {
  let state = initial;
  const subs = new Set();
  const keySubs = new Map();

  function get() { return state; }

  function set(patch) {
    const next = typeof patch === 'function' ? patch(state) : { ...state, ...patch };
    const changed = [];
    for (const k in next) if (next[k] !== state[k]) changed.push(k);
    for (const k in state) if (!(k in next)) changed.push(k);
    if (!changed.length) return state;
    const prev = state;
    state = next;
    for (const fn of subs) fn(state, prev, changed);
    for (const k of changed) {
      const s = keySubs.get(k);
      if (s) for (const fn of s) fn(state[k], prev[k], state);
    }
    return state;
  }

  function subscribe(fn) { subs.add(fn); return () => subs.delete(fn); }

  function on(key, fn) {
    if (!keySubs.has(key)) keySubs.set(key, new Set());
    keySubs.get(key).add(fn);
    return () => keySubs.get(key).delete(fn);
  }

  return { get, set, subscribe, on };
}

export const INITIAL = {
  sections: [],
  drugs: [],
  spec: null,             // the live RunSpec
  specText: '',           // the editor buffer
  specTextValid: true,    // JSON parses
  validation: null,       // {ok, errors[], warnings[]}
  validating: false,
  runs: [],
  activeRunId: null,      // the run whose geometry is loaded
  activeRun: null,
  submitting: false,
  status: null,           // live poll payload for the tracked run
  trackedRunId: null,
  compareA: null,
  compareB: null,
  fieldName: null,        // null = no overlay
  field: null,
  geometryMeta: null,
  hovered: null,
  pinned: null,
  theme: null,
};
