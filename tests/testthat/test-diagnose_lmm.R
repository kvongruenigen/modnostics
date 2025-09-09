test_that("diagnose_lmm works", {
  lmm <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
  expect_s3_class(diagnose_lmm(lmm), "shiny.appobj")
})
