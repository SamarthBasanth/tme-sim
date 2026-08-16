import { defineConfig, globalIgnores } from "eslint/config";
import nextVitals from "eslint-config-next/core-web-vitals";
import nextTs from "eslint-config-next/typescript";

const eslintConfig = defineConfig([
  ...nextVitals,
  ...nextTs,
  // Override default ignores of eslint-config-next.
  globalIgnores([
    // Default ignores of eslint-config-next:
    ".next/**",
    "out/**",
    "build/**",
    "next-env.d.ts",
    // Vendored third-party decoder: Emscripten output, redistributed
    // unmodified under Apache-2.0 (see public/draco/README.md). Linting it
    // reports on Google's build tooling rather than on this codebase, and the
    // one thing worth checking about it -- that it is byte-for-byte the file
    // we reviewed -- is checked by digest in scripts/audit-origins.mjs.
    // Scoped to this directory, not to public/ as a whole, so that a script
    // dropped anywhere else under public/ is still linted.
    "public/draco/**",
  ]),
]);

export default eslintConfig;
