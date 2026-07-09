test_that("diagnose_lmm works", {
  assign("sleepstudy_test", lme4::sleepstudy, envir = .GlobalEnv)
  on.exit(rm("sleepstudy_test", envir = .GlobalEnv), add = TRUE)

  lmm <- lme4::lmer(Reaction ~ Days + (1|Subject), data = sleepstudy_test)
  expect_s3_class(diagnose_lmm(lmm), "shiny.appobj")
})

test_that("diagnose_lmm works with namespaced datasets", {
  lmm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  expect_s3_class(diagnose_lmm(lmm), "shiny.appobj")
})

test_that("diagnose_lmm works with inline data expressions", {
  lmm <- lme4::lmer(
    Reaction ~ Days + (1 | Subject),
    data = subset(lme4::sleepstudy, Days >= 0)
  )

  expect_s3_class(diagnose_lmm(lmm), "shiny.appobj")
})

test_that("diagnose_lmm rejects unsupported input types", {
  expect_error(
    diagnose_lmm(mtcars),
    "requires a fitted `lmer` model"
  )
})

test_that("diagnose_lmm rejects models with multiple random-effects terms", {
  sleepstudy <- lme4::sleepstudy
  sleepstudy$Batch <- factor(rep(seq_len(18), each = 10))

  lmm <- suppressWarnings(
    lme4::lmer(
      Reaction ~ Days + (1 | Subject) + (1 | Batch),
      data = sleepstudy
    )
  )

  expect_error(
    diagnose_lmm(lmm),
    "exactly one random-effects term"
  )
})

test_that("diagnose_lmm rejects random-slope models", {
  sleepstudy <- lme4::sleepstudy
  lmm <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)

  expect_error(
    diagnose_lmm(lmm),
    "intercept-only random-effects terms"
  )
})

test_that("diagnose_lmm accepts lmerTest models", {
  sleepstudy <- lme4::sleepstudy
  lmm <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

  expect_s3_class(diagnose_lmm(lmm), "shiny.appobj")
})

test_that("diagnose_lmm rejects models without grouping-variable metadata", {
  sleepstudy <- lme4::sleepstudy
  lmm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

  broken_lmm <- lmm
  broken_lmm@flist <- list()

  expect_error(
    diagnose_lmm(broken_lmm),
    "single grouping variable named"
  )
})
