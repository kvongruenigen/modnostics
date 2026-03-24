test_that("diagnose_lmm works", {
  sleepstudy <- lme4::sleepstudy
  lmm <- lme4::lmer(Reaction ~ Days + (1|Subject), data = sleepstudy)
  expect_s3_class(diagnose_lmm(lmm), "shiny.appobj")
})
