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
