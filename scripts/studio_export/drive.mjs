// drive.mjs — interaction test for results/studio_export.html.
//
// Launches the installed Chrome headless, opens the export from file://, and
// drives it over the DevTools Protocol: switch runs, select each field
// overlay, toggle a channel, zoom, and pick a cell. Anything that only works
// because a server was answering fails here.
//
// Run: node scripts/studio_export/drive.mjs [outdir] [file]

import { spawn } from 'node:child_process';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const OUT = process.argv[2] || path.join(HERE, '..', '..', '.export-drive');
fs.mkdirSync(OUT, { recursive: true });

let FILE = path.resolve(process.argv[3]
  || path.join(HERE, '..', '..', 'results', 'studio_export.html'));

// A body-level fragment (the Artifact build) has no <!doctype>/<html>/<body>.
// Browsers parse it fine on its own, but the honest test is the shape the
// publisher actually serves, so wrap it in the same skeleton it uses and drive
// THAT. The wrapper is a throwaway written next to the screenshots.
const raw = fs.readFileSync(FILE, 'utf8');
const isFragment = !/^\s*<!doctype/i.test(raw);
if (isFragment) {
  const wrapped = path.join(OUT, 'wrapped.html');
  fs.writeFileSync(wrapped,
    '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
    + '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
    + '<style>*,*::before,*::after{box-sizing:border-box}'
    + 'body{margin:0}</style>\n</head>\n<body>\n' + raw + '\n</body>\n</html>\n');
  console.log(`fragment detected -> wrapped in a publisher-style skeleton at ${wrapped}`);
  FILE = wrapped;
}

const CHROME = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const PORT = 9333 + (process.pid % 200);

const chrome = spawn(CHROME, [
  '--headless=new', '--no-sandbox', `--remote-debugging-port=${PORT}`,
  `--user-data-dir=${path.join(OUT, 'profile')}`,
  '--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader',
  '--window-size=1700,1050', '--hide-scrollbars',
  'about:blank',
], { stdio: ['ignore', 'ignore', 'pipe'] });

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function endpoint() {
  for (let i = 0; i < 80; i++) {
    try {
      const r = await fetch(`http://127.0.0.1:${PORT}/json/version`);
      return (await r.json()).webSocketDebuggerUrl;
    } catch { await sleep(250); }
  }
  throw new Error('chrome did not expose a debugger port');
}

let failures = 0;
const ok = (c, m) => { if (c) console.log('  ok   ' + m); else { failures++; console.log('  FAIL ' + m); } };

const url = await endpoint();
const ws = new WebSocket(url);
await new Promise((r, j) => { ws.onopen = r; ws.onerror = j; });

let id = 0;
const pending = new Map();
const netRequests = [];
const consoleErrors = [];
ws.onmessage = (ev) => {
  const m = JSON.parse(ev.data);
  if (m.id && pending.has(m.id)) { pending.get(m.id)(m); pending.delete(m.id); return; }
  if (m.method === 'Network.requestWillBeSent') netRequests.push(m.params.request.url);
  if (m.method === 'Runtime.exceptionThrown') {
    consoleErrors.push(m.params.exceptionDetails.text + ' '
      + (m.params.exceptionDetails.exception || {}).description);
  }
};
const send = (method, params = {}, sessionId) => new Promise((res, rej) => {
  const n = ++id;
  pending.set(n, (m) => (m.error ? rej(new Error(method + ': ' + m.error.message)) : res(m.result)));
  ws.send(JSON.stringify({ id: n, method, params, sessionId }));
});

const { targetId } = await send('Target.createTarget', { url: 'about:blank' });
const { sessionId } = await send('Target.attachToTarget', { targetId, flatten: true });
const S = (m, p) => send(m, p, sessionId);
await S('Page.enable'); await S('Runtime.enable'); await S('Network.enable');
await S('Emulation.setDeviceMetricsOverride',
  { width: 1700, height: 1050, deviceScaleFactor: 1, mobile: false });

await S('Page.navigate', { url: 'file://' + FILE });
await sleep(6000);

const evalIn = async (expr) => {
  const r = await S('Runtime.evaluate', { expression: expr, awaitPromise: true, returnByValue: true });
  if (r.exceptionDetails) throw new Error(r.exceptionDetails.text + ' :: ' + expr.slice(0, 80));
  return r.result.value;
};
const shot = async (name) => {
  const { data } = await S('Page.captureScreenshot', { format: 'png' });
  fs.writeFileSync(path.join(OUT, name), Buffer.from(data, 'base64'));
};

console.log('boot');
ok(await evalIn('!!window.__tmesim && !!window.__tmesim.viewer'), 'viewer constructed (WebGL2 context acquired)');
ok(await evalIn('window.__tmesim.viewer.geom.n') === 40887, 'geometry loaded: 40,887 cells in the VBO');
ok(await evalIn('document.querySelector(".xbanner") !== null'), 'static-export banner present');
ok(await evalIn('/static export/.test(document.getElementById("hdr-backend").textContent)'),
  'top bar says "static export · no backend"');
ok(await evalIn('/1:1/.test(document.querySelector(".topbar").textContent)'), 'top bar states the 1:1 sampling rate');
ok(await evalIn('document.querySelectorAll(".caveat").length') >= 7, 'caveats rendered in full');
ok(await evalIn('document.body.textContent.includes("FALSIFIED AS VALUE-ADD")'),
  'the falsification caveat is on screen');

console.log('dead controls');
ok(await evalIn('!Array.from(document.querySelectorAll("button")).some(b=>/^Submit run$/.test(b.textContent.trim()))'),
  'no Submit button anywhere in the DOM');
ok(await evalIn('!Array.from(document.querySelectorAll("button")).some(b=>/^Cancel/.test(b.textContent.trim()))'),
  'no Cancel button anywhere in the DOM');
ok(await evalIn('!document.querySelector("#help-dlg tbody").textContent.includes("submit the current RunSpec")'),
  'the Ctrl+Enter submit shortcut is gone from the keyboard help');
ok(await evalIn('document.body.textContent.includes("Run submission is disabled")'),
  'the submit row is replaced by an explanation');
ok(await evalIn('document.querySelector(".xlog").textContent.includes("[stage] simulate")'),
  'the recorded engine log is shown for the open run');

console.log('field overlays');
for (const f of ['vascular_density', 'ecm_density', 'oxygen', 'drug_auc']) {
  await evalIn(`(()=>{const s=document.getElementById("f-overlay");s.value=${JSON.stringify(f)};`
    + 's.dispatchEvent(new Event("change"));return 1;})()');
  await sleep(900);
  const st = await evalIn('(()=>{const f=window.__tmesim.store.get().field;'
    + 'return f?{name:window.__tmesim.store.get().fieldName,rows:f.rows,cols:f.cols,'
    + 'min:f.min,max:f.max,n:f.values.length}:null;})()');
  ok(st && st.name === f && st.rows === 424 && st.cols === 675 && st.n === 424 * 675,
    `${f}: 424x675 grid decoded and uploaded (min ${st && st.min.toExponential(2)}, max ${st && st.max.toExponential(2)})`);
  ok(await evalIn('!document.querySelector(".legend").classList.contains("hidden")'),
    `${f}: legend visible with units`);
}
await shot('overlay-drug-auc.png');

console.log('interaction');
await evalIn('window.__tmesim.viewer.setColorMode(1); window.__tmesim.viewer.dirty=true;');
await sleep(500);
ok(await evalIn('window.__tmesim.viewer.colorMode') === 1, 'cells recoloured by field value');
await shot('colour-by-field.png');

// Channel visibility is culled in the vertex stage (a zero-area point), so
// stats.drawn does not move — the only honest check is that the pixels do.
const px = async () => {
  const { data } = await S('Page.captureScreenshot', { format: 'png' });
  return createHash('sha256').update(data).digest('hex');
};
const before = await px();
await evalIn('window.__tmesim.viewer.setVisible(0,false);');
await sleep(700);
const after = await px();
ok(await evalIn('window.__tmesim.viewer.visible[0]') === 0 && before !== after,
  'TUMOR channel toggled off and the rendered image changed');
await evalIn('window.__tmesim.viewer.setVisible(0,true);');
await sleep(500);

await evalIn('window.__tmesim.viewer.zoomAt(600,400,4); window.__tmesim.viewer.dirty=true;');
await sleep(600);
ok(await evalIn('window.__tmesim.viewer.cam.zoom') > 0.3, 'zoom applied');
await shot('zoomed.png');

const info = await evalIn(
  '(()=>{const v=window.__tmesim.viewer;let i=-1;'
  + 'for(let px=200;px<900&&i<0;px+=17)for(let py=150;py<700&&i<0;py+=17)i=v.pick(px,py);'
  + 'if(i<0)return null;const c=v.cellInfo(i);'
  + 'return{index:c.index,type:c.typeName,x:c.x,y:c.y,field:c.fieldValue,total:c.composition.total};})()');
ok(!!info, `cell picking works: #${info && info.index} ${info && info.type} at `
  + `(${info && info.x.toFixed(0)}, ${info && info.y.toFixed(0)}) um, drug_auc `
  + `${info && info.field !== null ? info.field.toExponential(2) : 'null'}, `
  + `${info && info.total} neighbours in radius`);

console.log('run switching');
const ids = await evalIn('window.__tmesim.store.get().runs.map(r=>r.run_id)');
for (const rid of ids) {
  await evalIn(`(()=>{const t=Array.from(document.querySelectorAll(".dock tbody tr"))`
    + `.find(tr=>tr.textContent.includes(${JSON.stringify(rid)}));t.click();return 1;})()`);
  await sleep(1200);
  const cur = await evalIn('(()=>{const r=window.__tmesim.store.get().activeRun;'
    + 'return r?{id:r.run_id,drug:r.spec.drug,kill:r.metrics.kill_fraction,'
    + 'cav:r.provenance.caveats.length,log:(r.log_tail||[]).length}:null;})()');
  ok(cur && cur.id === rid && cur.cav >= 7,
    `${rid} (${cur && cur.drug}): kill ${(cur.kill * 100).toFixed(2)}%, `
    + `${cur.cav} caveats, ${cur.log} log lines`);
}
await shot('final.png');

console.log('theme');
await evalIn('document.getElementById("btn-theme").click()');
await sleep(600);
ok(await evalIn('document.documentElement.getAttribute("data-theme")') === 'light', 'light theme applies');
await shot('light.png');
await evalIn('document.getElementById("btn-theme").click()');

console.log('network');
const external = netRequests.filter((u) => !u.startsWith('file://') && !u.startsWith('data:'));
ok(external.length === 0,
  `${netRequests.length} network events, 0 non-file/data requests`
  + (external.length ? ' -> ' + external.slice(0, 5).join(', ') : ''));
ok(consoleErrors.length === 0, 'no uncaught exceptions'
  + (consoleErrors.length ? ': ' + consoleErrors.slice(0, 2).join(' | ') : ''));

ws.close();
chrome.kill();
console.log(failures ? `\n${failures} FAILURE(S)` : `\nall checks passed — screenshots in ${OUT}`);
process.exit(failures ? 1 : 0);
