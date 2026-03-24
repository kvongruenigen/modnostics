# Modnostics Improvement Backlog

This backlog turns the audit into a sequence of small, safe, professional-grade pull requests.

Guiding principles:
- Preserve current functionality unless a task explicitly documents a behavior change.
- Prefer incremental refactors over rewrites.
- Land release-readiness and correctness fixes before structural cleanup.
- Add tests before or alongside risky refactors.

## Current audit summary

Observed issues:
- `diagnose_lmm()` is a large monolithic function that mixes model parsing, data shaping, plotting, UI, and server logic.
- The function relies on brittle assumptions about model structure and data lookup.
- Package hygiene is incomplete: `R CMD check --no-manual .` currently errors, `LICENSE` is missing, and local check artifacts are present.
- Namespace usage is inconsistent: many `library()` calls happen inside the exported function and `NAMESPACE` uses `exportPattern`.
- Test coverage is minimal and misses important edge cases.
- There is dead code and deprecated plotting API usage.

## Recommended PR order

### PR 1: Package hygiene and release baseline

Goal:
- Make the repository buildable, checkable, and easier to work with.

Why it matters:
- This removes avoidable packaging noise before code refactors start.
- A passing package baseline reduces the risk of shipping regressions.

Suggested changes:
- Fix `DESCRIPTION` metadata so `R CMD check` passes cleanly.
- Add the missing `LICENSE` file or adjust the declared license.
- Ignore local artifacts like `..Rcheck/` and `.DS_Store`.
- Tighten `.Rbuildignore` as needed for non-package files.
- Expand `README.md` from a short note into a proper package entry point.

Affected files:
- `DESCRIPTION`
- `README.md`
- `.gitignore`
- `.Rbuildignore`
- `LICENSE`

Risk:
- Low

Estimated scope:
- Small

Acceptance criteria:
- `R CMD check --no-manual .` passes locally.
- No local check artifacts or editor junk are intended for version control.
- README explains purpose, installation, and basic usage.

### PR 2: Define supported input contract for `diagnose_lmm()`

Goal:
- Explicitly document what model shapes and inputs are supported.

Why it matters:
- The current function behaves as if it supports general `lmer` models, but the implementation assumes a narrower shape.
- Clear support boundaries prevent hidden bugs and ambiguous maintenance work.

Suggested changes:
- Validate `lmm` class and fail with clear messages for unsupported inputs.
- Decide whether the function supports:
  - only single-group random-intercept LMMs, or
  - broader `lmer` structures such as random slopes and multiple grouping factors.
- Update docs and tests to match the chosen contract.

Affected files:
- `R/linear-mixed-model.R`
- `man/diagnose_lmm.Rd`
- `README.md`
- `tests/testthat/test-diagnose_lmm.R`

Risk:
- Medium

Estimated scope:
- Small

Acceptance criteria:
- Unsupported inputs produce explicit, user-facing errors.
- Supported inputs are clearly documented and tested.

### PR 3: Fix brittle model and data extraction

Goal:
- Remove implementation assumptions that currently break valid model objects.

Why it matters:
- The current implementation breaks when the model uses namespaced data references like `lme4::sleepstudy`.
- Direct slot access and positional term assumptions are fragile and difficult to extend safely.

Suggested changes:
- Replace direct use of `lmm@call` and `lmm@frame` with more stable access patterns where possible.
- Derive response, fixed predictors, grouping variables, and model frame from dedicated helper functions.
- Remove assumptions that the last term in the formula is the grouping variable.
- Make random-effects and influence calculations robust to the supported model shapes.

Affected files:
- `R/linear-mixed-model.R`
- new helper files under `R/`
- `tests/testthat/test-diagnose_lmm.R`

Risk:
- Medium

Estimated scope:
- Medium

Acceptance criteria:
- The current failing namespaced dataset case passes.
- Supported model metadata is extracted through helper functions rather than inline slot parsing.
- No new regressions in the basic dashboard flow.

### PR 4: Decompose `diagnose_lmm()` into testable helpers

Goal:
- Improve readability, reuse, and maintainability without changing the public API.

Why it matters:
- A single large function is hard to test and discourages safe iteration.
- Splitting responsibilities will make future functionality additions much lower risk.

Suggested changes:
- Keep `diagnose_lmm()` as the public entry point.
- Extract internal helpers for:
  - metadata extraction
  - summary table construction
  - plot construction
  - Shiny UI layout
  - server output wiring
- Move nested helper functions to top-level internal functions.

Affected files:
- `R/linear-mixed-model.R`
- new helper files such as:
  - `R/diagnose-lmm-helpers.R`
  - `R/diagnose-lmm-plots.R`
  - `R/diagnose-lmm-ui.R`

Risk:
- Medium

Estimated scope:
- Medium

Acceptance criteria:
- `diagnose_lmm()` is short and orchestration-focused.
- Internal helpers have single responsibilities and readable names.
- Tests cover the extracted helpers where practical.

### PR 5: Standardize namespace and dependency usage

Goal:
- Make imports explicit and remove hidden side effects.

Why it matters:
- Calling `library()` inside an exported function changes runtime behavior in surprising ways.
- `exportPattern` is broad and becomes risky as more internal helpers are added.

Suggested changes:
- Remove `library()` calls from package code.
- Use explicit namespace qualification or roxygen-managed imports consistently.
- Replace `exportPattern("^[[:alpha:]]+")` with explicit exports.
- Review `DESCRIPTION` imports and remove anything unused.

Affected files:
- `R/linear-mixed-model.R`
- all new or split files under `R/`
- `DESCRIPTION`
- `NAMESPACE`

Risk:
- Low

Estimated scope:
- Small

Acceptance criteria:
- No package code depends on attaching packages at runtime.
- Exported functions are explicit.
- Imports reflect actual usage.

### PR 6: Remove dead code and deprecated API usage

Goal:
- Reduce noise and keep the package compatible with current package versions.

Why it matters:
- Dead code hides intent and makes maintenance harder.
- Deprecation warnings are future failures in waiting.

Suggested changes:
- Remove unused variables and stale comments.
- Replace deprecated `ggplot2` calls such as `geom_errorbarh()` and `aes_string()`.
- Standardize local naming for tables, plots, and helper outputs.

Affected files:
- `R/linear-mixed-model.R`
- helper files under `R/`
- tests as needed

Risk:
- Low

Estimated scope:
- Small

Acceptance criteria:
- No obvious unused local variables remain.
- Tests run without deprecation warnings from package code.
- Naming is consistent across helpers and outputs.

### PR 7: Strengthen test coverage and contributor workflow

Goal:
- Make refactors safer and improve developer experience.

Why it matters:
- The current single happy-path test is not enough to support ongoing maintenance.
- Basic contributor automation pays off quickly even in a small package.

Suggested changes:
- Add tests for:
  - namespaced dataset references
  - invalid inputs
  - supported model variants
  - expected return type
- Add a simple CI workflow that runs package checks and tests.
- Optionally add lint/style tooling if the team wants enforced consistency.

Affected files:
- `tests/testthat/test-diagnose_lmm.R`
- new test files under `tests/testthat/`
- `.github/workflows/`
- optional lint config

Risk:
- Low

Estimated scope:
- Medium

Acceptance criteria:
- CI runs package checks automatically.
- Test coverage reflects the documented support contract.
- Contributors can tell quickly whether a change is safe.

### PR 8: UX and functionality improvements after stabilization

Goal:
- Improve usefulness once the baseline code quality work is complete.

Why it matters:
- There are clear opportunities to improve the dashboard, but they should follow stabilization work.
- Some of these changes require product decisions, not just engineering effort.

Candidate improvements:
- Add interpretive guidance for diagnostics.
- Add report or HTML export.
- Add additional model diagnostics where appropriate.
- Improve layout and presentation for larger models.

Affected files:
- Shiny UI/server files after extraction
- docs and README
- tests

Risk:
- Medium

Estimated scope:
- Medium to Large

Acceptance criteria:
- Any new functionality is documented, tested, and explicitly scoped.
- Product-facing behavior changes are approved before implementation.

## Quick wins

- Fix package metadata and license setup.
- Ignore `..Rcheck/` and `.DS_Store`.
- Remove `library()` calls from package code.
- Replace deprecated plotting calls.
- Remove unused variables and stale comments.
- Add tests for the currently failing namespaced data case.

## Deeper architectural improvements

- Split the package into analysis helpers, plot builders, and Shiny composition.
- Establish a clear support contract for model complexity.
- Add CI and optional linting/style gates.
- Introduce export/reporting only after product decisions are made.

## Product or design decisions needed

- Should support be limited to simple random-intercept LMMs, or broadened to more general `lmer` models?
- Should the dashboard remain a technical diagnostic tool, or include interpretive guidance for users?
- Is saving/exporting the dashboard or a report a requirement?

## Definition of done

A professional-quality baseline for this repository means:
- Package checks pass cleanly.
- Public behavior is documented and tested.
- The exported API is explicit and stable.
- The codebase avoids hidden side effects and dead code.
- Core logic is split into maintainable units.
- Repo hygiene supports normal development and review workflows.
- Future enhancements can be added without rewriting the package again.
