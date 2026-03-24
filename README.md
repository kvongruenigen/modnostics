# modnostics

`modnostics` is an R package for exploring diagnostics from linear mixed-effects
models fit with `lme4`.

The package currently centers on `diagnose_lmm()`, which returns a Shiny app
that brings together:
- model fit statistics
- fixed and random effect summaries
- residual diagnostic plots
- influence diagnostics
- predictor effect visualizations

Current supported model scope:
- fitted `lmer` models
- exactly one random-effects term
- intercept-only random effects of the form `(1 | group)`
- a model `data` argument supplied as a named object

## Installation

The package is not published on CRAN yet. For local development:

```r
renv::restore()
devtools::load_all()
```

## Example

```r
library(lme4)
library(modnostics)

model <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

app <- diagnose_lmm(model)
app
```

## Development

Useful local checks:

```r
testthat::test_local(".")
```

```sh
R CMD check --no-manual .
```

## Current status

The package is under active refactoring. See `BACKLOG.md` for the current
improvement plan and `AGENTS.md` for repo-specific implementation guidance.
