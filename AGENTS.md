# AGENTS.md

This file gives repo-specific instructions to coding agents working in `modnostics`.

## Project overview

- `modnostics` is an R package for diagnostic inspection of linear mixed models.
- The current public surface is very small: the main exported function is `diagnose_lmm()`.
- Most logic currently lives in a single file: `R/linear-mixed-model.R`.
- The package includes Shiny UI code, statistical summaries, and plot generation in the same function.

## Primary goals for work in this repo

When making changes, optimize for:
- maintainability
- readability
- consistency with normal R package conventions
- minimal disruption to current behavior
- small, safe, reviewable pull requests

Prefer incremental refactors over broad rewrites.

## Current repository priorities

The current recommended order of work is:
1. package hygiene and `R CMD check` baseline
2. input validation and supported-model contract
3. robust model/data extraction
4. helper extraction and modularization
5. namespace/import cleanup
6. dead-code and deprecation cleanup
7. stronger tests and CI
8. optional UX and functionality improvements

See [BACKLOG.md](/Users/kayvongrunigen/Work/Projects/modnostics/BACKLOG.md) for the full plan.

## Rules for agents

### 1. Preserve behavior unless a task explicitly allows changes

- Do not change the visible behavior of `diagnose_lmm()` unless the task clearly calls for it.
- If you need to narrow or clarify supported behavior, document it and add tests.
- Avoid speculative redesigns.

### 2. Fix packaging and correctness before polish

- Prioritize broken checks, brittle behavior, and weak tests ahead of visual or stylistic improvements.
- Do not spend time on new features while core package hygiene is still failing.

### 3. Follow standard R package practices

- Do not call `library()` inside exported package functions.
- Prefer explicit namespace usage like `pkg::fun()` or roxygen-managed imports.
- Avoid direct slot access like `@call` or `@frame` when a more stable accessor exists.
- Keep exports explicit; do not rely on broad `exportPattern()` behavior as the package grows.

### 4. Keep refactors small and test-backed

- Break larger work into PR-sized steps.
- For risky refactors, add or update tests before fully reshaping the implementation.
- Do not mix unrelated cleanup into the same change unless it is necessary to complete the task safely.

### 5. Respect local user changes

- The worktree may contain unrelated files such as `Progress.Rmd`, `.DS_Store`, or `..Rcheck/`.
- Do not revert user changes unless explicitly asked.
- Ignore unrelated dirty-state files when possible.

## Known codebase pitfalls

Agents should be aware of these existing issues:

- `diagnose_lmm()` currently assumes model structure in brittle ways.
- Namespaced data references such as `lme4::sleepstudy` can fail in the current implementation.
- The last formula term is currently treated as the grouping variable; that is not a safe general assumption.
- The code currently mixes package attachment, data extraction, plotting, and Shiny composition in one place.
- There is minimal test coverage and at least one real functional failure discovered under `testthat::test_local(".")`.
- Some code uses deprecated `ggplot2` patterns such as `geom_errorbarh()` and `aes_string()`.

## Preferred change patterns

### For package-code changes

- Keep `diagnose_lmm()` as the public entry point unless there is a compelling reason not to.
- Extract pure internal helpers for:
  - model metadata extraction
  - table building
  - plot building
  - UI composition
  - server wiring
- Name helpers clearly and keep each helper focused on one responsibility.

### For tests

When changing behavior around supported inputs, add tests for:
- a valid baseline LMM
- invalid input types
- namespaced data references
- supported model structures
- return type and basic Shiny object creation

### For documentation

- Keep `README.md` useful for someone seeing the package for the first time.
- Update man docs when the function contract changes.
- If behavior is intentionally unsupported, say so clearly rather than leaving it implicit.

## Validation checklist

Before finishing a change, prefer running the smallest relevant validation set:

- targeted tests during iteration
- `Rscript -e 'testthat::test_local(\".\")'` for package-level tests
- `R CMD check --no-manual .` for release-readiness changes

If you cannot run a command, say so clearly in the handoff.

## File-level guidance

- `R/linear-mixed-model.R`
  Treat as the current integration point. Avoid making it even larger.

- `tests/testthat/`
  Expand coverage incrementally as refactors land.

- `DESCRIPTION` and `NAMESPACE`
  Treat changes here as important package-surface changes; keep them deliberate and minimal.

- `README.md`
  Prefer concise, practical package documentation over research notes.

## What not to do

- Do not rewrite the package into a different architecture in one pass.
- Do not introduce new frameworks or tooling without a clear payoff.
- Do not add feature work before stabilizing packaging and tests.
- Do not silently broaden support for model shapes without tests and documentation.

## Good end state

A professional-quality change in this repository should leave the package:
- easier to read
- easier to test
- more explicit about supported inputs
- more standard in its package conventions
- safer to extend in future PRs
