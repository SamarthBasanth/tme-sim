# Auth contract — the login portal

**Frozen 2026-10-09.** Written before any code, so several agents can build
against it at once without reading each other's diffs. If you need to deviate,
say so in your report rather than deviating quietly.

This contract AMENDS `docs/WORKSPACE_CONTRACT.md` §9, which says sign-up and
user accounts are not authorised. That clause recorded a design decision, not a
law; the product owner has since asked for a login portal, so the clause is
superseded for identity and NOT for the rest of §9.

---

## 0. Why there is a key at all, and why it does not go away

The instinct behind this work is right: a shared secret pasted into a prompt is
not a login, and a visitor should not meet one. But the key is NOT only a login
substitute — it is the credential for every non-browser caller:

* `.github/workflows/docker.yml` drives the built image with `X-API-Key`.
* `.github/workflows/product.yml` asserts the gate with it.
* `render.yaml` generates it (`generateValue: true`) and the server's own
  startup banner prints it.
* `docs/DEPLOYMENT.md` and `docs/RUNBOOK_DEPLOY.md` document `curl` flows.

Deleting it would break all of that and leave scripts with no way in. So:

> **Two credentials, one gate. A browser session for people, `X-API-Key` for
> programs. Neither is removed, and the gate accepts either.**

A reader who only ever uses the key notices no change. That is this contract's
compatibility bar, and it is the same bar §8 of the workspace contract set.

---

## 1. Routes

| Route | Method | Credential | Purpose |
|---|---|---|---|
| `/login` | GET | none | the page. Serves sign-in and the first-run claim. |
| `/api/auth/me` | GET | session **or** key **or** none | who am I; `{"authenticated": false}` when nobody. Never 401 — the page asks this to decide what to render. |
| `/api/auth/login` | POST | none | email + password → sets the cookie. |
| `/api/auth/logout` | POST | session | clears it. |
| `/api/auth/signup` | POST | **bootstrap token** | creates an account. See §4. |
| `/api/auth/password` | POST | session | change own password (old + new). |

`/login` and `/api/auth/*` are the only new ungated paths. Everything else in
`/api/` keeps the existing gate, now satisfied by **either** credential.

## 2. The session cookie

Name `vivum_session`. `HttpOnly`, `SameSite=Lax`, `Path=/`, and `Secure`
whenever the request did not arrive on loopback. Value is
`<user_id>.<expiry_unix>.<hmac>` where the HMAC is SHA-256 over
`"<user_id>.<expiry>"` keyed on a server secret.

* **The secret** is `TMESIM_SESSION_SECRET` if set, else derived from the user
  store's own random `instance_secret`, generated once on first write. It must
  NOT be derived from `TMESIM_STUDIO_API_KEY` — rotating the key would then
  silently invalidate every session, and the two secrets have different
  lifetimes on purpose.
* **Lifetime** 14 days. There is no refresh-on-use in v1; re-login is cheap and
  a sliding window is a thing to get wrong quietly.
* Compared with `secrets.compare_digest`, always. Never `==`.

## 3. Passwords

`hashlib.scrypt` from the standard library — **no new dependency**, which
matters because `constraints.txt` exists to pin numerically load-bearing
packages and adding `bcrypt`/`argon2` to a scientific pin set for a login form
is a poor trade.

Stored as `scrypt$<n>$<r>$<p>$<salt_b64>$<hash_b64>` with `n=2**15, r=8, p=1`
and a 16-byte `secrets.token_bytes` salt. The parameters live IN the string so
they can be raised later without invalidating existing hashes.

Minimum length 12 characters. No composition rules — they push people toward
`Password1!` and the length floor does more. Reject anything over 1024 bytes
before hashing (scrypt on a megabyte of input is a denial of service).

## 4. The first account, and why sign-up is not open

An internet-reachable instance with open sign-up is claimed by whoever finds it
first. So `POST /api/auth/signup` requires a **bootstrap token**, and the
bootstrap token is the operator's existing `TMESIM_STUDIO_API_KEY`.

This is deliberate reuse rather than a new mechanism: the operator already has
that value (Render prints it in the dashboard; the banner prints it locally), it
is already the thing that proves "I run this server", and §9 of the workspace
contract forbids *a new API-key mechanism* — this adds none.

After the first account exists, sign-up still requires the token. There is no
invite-by-email in v1 (see §7).

## 5. What a session is NOT

**Projects, runs and specimens stay SHARED.** A session identifies the person
at the keyboard; it does not partition the stores. `GET /api/workspace` returns
the same workspace to every authenticated user.

Stated this plainly because the half-built version is the dangerous one: a UI
that says "your projects" over a shared store teaches a false belief about who
can see what. So, until a separate decision says otherwise:

* no "your projects", no owner column, no avatars — the §9 prohibition on
  implying per-user data **still stands**
* `GET /api/account` keeps telling the truth about what is shared
* the login page says, in its own words, that an account is a way in rather
  than a private space

## 6. Storage

`<state_dir>/auth/users.json`, created `0700`/`0600` via
`cache.private_makedirs` and `PRIVATE_FILE_MODE`. Under `TMESIM_STATE_DIR` so
accounts survive a redeploy on a host that has a disk — and so that a host
WITHOUT one loses accounts exactly as it loses projects, which the home page
already discloses.

Shape:

```json
{
  "version": 1,
  "instance_secret": "<hex>",
  "users": [
    {"id": "u_00001", "email": "...", "password": "scrypt$...",
     "created_at": 0.0, "last_login_at": null, "disabled": false}
  ]
}
```

Email is stored casefolded and matched casefolded. One file, one lock, same
discipline as `ProjectStore`.

## 7. What this contract does NOT authorise

* **No email.** No verification mail, no reset link, no SMTP config. A
  forgotten password is reset by the operator editing the store or by a CLI
  subcommand. Half-built email is worse than none.
* **No OAuth / SSO.**
* **No roles or permissions.** Every account can do everything the key can do.
  A role model with one role is just a column that lies later.
* **No per-user data partitioning** (see §5).
* **No rate-limit-free login.** Failed attempts per email are throttled; see
  the implementation's own note.
* **Deleting or weakening the `X-API-Key` path.**

## 8. Honesty requirements, carried over

Unchanged and absolute: nothing here may imply clinical use, diagnostic use,
treatment planning, a "digital twin", personalised medicine, or AI/ML
capability. Tissue registers **by filesystem path**, never "upload". Nothing
under `results/` or `archive/` is rewritten.

---

## 9. Amendment — 2026-10-10: the operator reset path, and the password form

**APPENDED, NOT AN EDIT.** Everything above this line is the frozen text of
2026-10-09 and is byte-identical to it. This section records two things the
frozen contract described without providing, and it adds no new mechanism: it
names where §7's operator reset now lives and pins the body shape §1 left
open. Nothing in §0–§8 is superseded.

### 9.1 §7's operator reset is `python -m app.server account`

§7 says a forgotten password "is reset by the operator editing the store or by
a CLI subcommand", and there was no subcommand — so the only route was
hand-editing `users.json`, which means hand-computing a password hash. These
are that subcommand:

```
python -m app.server account list
python -m app.server account reset-password --email <address> [--stdin]
python -m app.server account disable --email <address>
python -m app.server account enable  --email <address>
python -m app.server account rotate-session-key --yes
```

No subcommand still starts the server, so every existing caller — `vivum
studio`, the `Dockerfile`, `render.yaml`, both documented deploy flows — is
unchanged.

**The authority is filesystem ownership, and that is the whole of it.** There
is no route, because the only authority an HTTP caller could present is
`TMESIM_STUDIO_API_KEY`, and §4 already spends that as the sign-up bootstrap
token: a value printed in a hosting dashboard and pasted into CI must not also
become a master credential over identities. §6 already makes "can edit the
account store" and "is the operator on this host" the same sentence (`0600`
under `0700`), so these commands inherit that and add nothing to it.

Three consequences of that, stated so they are not rediscovered:

* **no `--password` flag, ever.** A password on an argv is readable by every
  user on the host through `ps`, is appended to the shell's history file, and
  is printed into any orchestration tool's job log. The password is read from
  a no-echo prompt (asked twice) or from a pipe with `--stdin`.
* **these commands do not create accounts.** §4 governs account creation and
  requires the bootstrap token; a reset that created an account on a typo'd
  address would be a second way in that §4 never authorised.
* **the store is named on every invocation.** §6's fallback when
  `TMESIM_STATE_DIR` is unset is `app/auth/`, so an operator whose server has
  that variable set and whose shell does not would otherwise reset a password
  in a file no server reads — and find out at the sign-in that refuses it.

### 9.2 What a reset does and does not end

§2 says there is no refresh-on-use and no per-session record. The consequence
is sharper than it looks and the CLI says it in full:

| action | ends live sessions? |
|---|---|
| `reset-password` | **No.** None. Every `vivum_session` cookie already issued for that account keeps opening the gate until its own 14-day expiry. |
| `disable` | **Yes, that one account, at once** — `disabled` is re-read on every request, so the next request a held cookie makes is refused. |
| `enable` | Lets those same unexpired cookies back in, so it is not a clean undo. |
| `rotate-session-key` + restart | **Yes, every session on the deployment**, every account. There is no narrower version. |

A reset therefore does **not** evict somebody who already has a session. If
that is the reason for the reset, `disable` is the lever. This is a property of
§2's stateless signed token, not a defect to be fixed in the CLI, and the
refusal to pretend otherwise is the same refusal `POST /api/auth/password`
already makes with `other_sessions_invalidated: false`.

`rotate-session-key` is refused outright when `TMESIM_SESSION_SECRET` is set:
§2 gives that variable precedence, so the store's `instance_secret` would never
be consulted and the command would print a success having ended nothing.

### 9.3 `POST /api/auth/password`'s body shape

§1 pinned the route and not the body, which is why the portal shipped without
this form — a form built on a guess is how two agents disagree in production.
The route is written, so the shape is now read off it rather than guessed, and
pinned here:

```
POST /api/auth/password          credential: session cookie ONLY, never the key
{"old_password": "<string>", "new_password": "<string>"}

200 {"ok": true, "other_sessions_invalidated": false,
     "other_sessions_note": "<the sentence the UI renders verbatim>"}
4xx {"detail": {"message": "<sentence>", "reason"?: "...",
                "retry_after_s"?: <int>}}
```

Both fields are required and must be non-empty strings; any other naming
(`old`, `current_password`, `password`) is a 400. The front end is
`app/web/js/auth/password.js`, behind the existing `account` control beside
logout, and it renders `other_sessions_note` verbatim rather than re-typing it.

### 9.4 One correction to a storage assumption, under §6

§6 puts the store under `TMESIM_STATE_DIR` and the implementation falls back to
`app/auth/` when that is unset. That fallback path was **not** in `.gitignore`,
while every sibling Studio store was — so `users.json`, which holds every
password hash *and* the `instance_secret` §2 derives the session signing key
from, was stageable by `git add -A`. It is ignored as of this amendment and
guarded by `tests/test_account_recovery.py`. No such file was ever committed.
