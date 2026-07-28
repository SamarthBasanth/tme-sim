// projects_static.js — the STATIC EXPORT's stand-in for js/ui/projects.js.
//
// Projects and specimens are mutable server state: creating a project, attaching
// a run to a specimen and running a cross-specimen comparison are all writes to
// app/server, and the comparison is computed there. A file cannot hold any of
// it, so the export ships none of it — including the toolbar button, which
// would otherwise be a control that looks live and answers 501.
//
// The static-export banner names this among the surfaces that were removed, so
// the absence is stated rather than silent.
//
// Same constructor signature and returned interface as the real module, so
// js/main.js is bundled unmodified.

export function createProjectsPanel() {
  return {
    el: null,
    button: null,
    async open() {},
    async reload() {},
    get projectId() { return null; },
  };
}
