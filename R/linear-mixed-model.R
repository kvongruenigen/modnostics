#' Diagnose a Linear Mixed Model
#'
#' This function launches a Shiny dashboard to inspect diagnostics of a
#' supported linear mixed-effects model.
#'
#' Supported inputs are currently limited to fitted `lmer` models with:
#' \itemize{
#'   \item exactly one random-effects grouping term,
#'   \item an intercept-only random-effects structure of the form
#'   \code{(1 | group)}.
#' }
#'
#' @param lmm A fitted `lmer` model that matches the currently supported model
#'   structure.
#' @return A Shiny app object.
#' @export

# Extract model metadata using stable accessors rather than direct slot parsing.
.get_diagnose_lmm_metadata <- function(lmm) {
  model_formula <- stats::formula(lmm)
  fixed_formula <- suppressWarnings(lme4::nobars(model_formula))
  model_frame <- stats::model.frame(lmm)
  random_groups <- names(lme4::getME(lmm, "flist"))
  response_var <- all.vars(fixed_formula)[1]
  fixed_vars <- setdiff(all.vars(fixed_formula), response_var)

  list(
    model_command = paste(deparse(getCall(lmm)), collapse = "\n"),
    model_formula = model_formula,
    model_frame = model_frame,
    response_var = response_var,
    fixed_vars = fixed_vars,
    grouping_var = random_groups[[1]]
  )
}

# Build a model copy that uses influence.ME's built-in `data.update` fallback,
# which refits from `model.frame(model)` instead of trying to resolve the
# original data expression.
.build_diagnose_lmm_influence_model <- function(lmm) {
  lmm_influence <- lmm
  lmm_influence@call$data <- quote(data.update)
  lmm_influence
}

# Validate the subset of lmer models currently supported by the dashboard.
.validate_diagnose_lmm_input <- function(lmm) {
  if (!inherits(lmm, c("lmerMod", "lmerModLmerTest"))) {
    stop(
      paste(
        "`diagnose_lmm()` requires a fitted `lmer` model from `lme4` or",
        "`lmerTest`."
      ),
      call. = FALSE
    )
  }

  model_formula <- stats::formula(lmm)
  random_terms <- suppressWarnings(lme4::findbars(model_formula))

  if (length(random_terms) != 1) {
    stop(
      paste(
        "`diagnose_lmm()` currently supports models with exactly one",
        "random-effects term."
      ),
      call. = FALSE
    )
  }

  random_term <- random_terms[[1]]

  if (!identical(random_term[[2]], quote(1))) {
    stop(
      paste(
        "`diagnose_lmm()` currently supports intercept-only random-effects",
        "terms of the form `(1 | group)`."
      ),
      call. = FALSE
    )
  }

  random_groups <- names(lme4::getME(lmm, "flist"))

  if (length(random_groups) != 1 || is.na(random_groups[[1]]) ||
      identical(random_groups[[1]], "")) {
    stop(
      paste(
        "`diagnose_lmm()` currently supports a single grouping variable named",
        "in the fitted model."
      ),
      call. = FALSE
    )
  }

  model_frame <- try(stats::model.frame(lmm), silent = TRUE)

  if (inherits(model_frame, "try-error")) {
    stop(
      paste(
        "`diagnose_lmm()` currently requires a fitted model with an",
        "accessible model frame."
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.get_diagnose_lmm_thresholds <- function() {
  list(
    p_value_highlight = 0.05,
    residual_sd_highlight = 10
  )
}

.get_diagnose_lmm_lmer_test_model <- function(lmm) {
  if (!inherits(lmm, "lmerModLmerTest")) {
    suppressMessages({
      lmm_lmerTest <- try(lmerTest::as_lmerModLmerTest(lmm), silent = TRUE)
    })
    # If conversion fails, just use the original model (no p-values)
    if (inherits(lmm_lmerTest, "try-error")) {
      lmm_lmerTest <- lmm
    }
  } else {
    lmm_lmerTest <- lmm
  }

  lmm_lmerTest
}

.build_diagnose_lmm_tables <- function(lmm, metadata, thresholds) {
  lmm_lmerTest <- .get_diagnose_lmm_lmer_test_model(lmm)

  fixed <- broom.mixed::tidy(lmm_lmerTest, effects = "fixed", conf.int = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  fixed_df <- fixed %>%
    mutate(conf_int = paste0("[", conf.low, ", ", conf.high, "]"))

  if ("p.value" %in% names(fixed_df)) {
    fixed_df <- fixed_df %>% select(term, estimate, conf_int, p.value)
  } else {
    fixed_df <- fixed_df %>% select(term, estimate, conf_int)
  }

  rand <- as.data.frame(VarCorr(lmm))[,c("vcov", "sdcor")] %>%
    mutate(
      group = c(paste0(metadata$grouping_var, " (Intercept)"), "Residual"),
      Variance = round(vcov, 3),
      Std_Dev = round(sdcor, 3)
    ) %>%
    select(group, Variance, Std_Dev)

  fit_stats <- performance::model_performance(lmm) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
  fit_stats_df <- as.data.frame(fit_stats) %>%
    dplyr::mutate(across(where(is.numeric), ~ round(.x, 2)))

  overfitting <- lme4::isSingular(lmm, tol = 1e-4)

  list(
    fixed = fixed,
    fixed_df = fixed_df,
    rand = rand,
    fit_stats_df = fit_stats_df,
    overfitting = overfitting,
    thresholds = thresholds
  )
}

.build_diagnose_lmm_effect_plots <- function(model, response_var) {
  eff <- predictorEffects(model)

  plots <- lapply(names(eff), function(var) {
    df <- as.data.frame(eff[[var]])

    if (is.numeric(df[[var]])) {
      p <- ggplot(df, aes_string(x = var, y = "fit")) +
        geom_line(color = "#B163FF", linewidth = 1) +
        geom_ribbon(
          aes(ymin = lower, ymax = upper),
          alpha = 0.2,
          fill = "#CCCCFF"
        ) +
        labs(
          title = paste("Effect of", var),
          x = var,
          y = paste0("Predicted ", response_var)
        )
    } else {
      p <- ggplot(df, aes_string(x = var, y = "fit")) +
        geom_errorbar(
          aes(ymin = lower, ymax = upper),
          width = 0.1,
          color = "#FF63D3"
        ) +
        geom_point(size = 3, color = "#B163FF") +
        labs(
          title = paste("Effect of", var),
          x = var,
          y = paste0("Predicted ", response_var)
        )
    }

    p + theme_minimal()
  })

  names(plots) <- names(eff)
  plots
}

.build_diagnose_lmm_plots <- function(lmm, metadata, fixed) {
  influence_model <- .build_diagnose_lmm_influence_model(lmm)

  plot_resid_fitted <- ggplot(data.frame(
    Fitted = fitted(lmm),
    Residuals = scale(resid(lmm))
  ), aes(x = Fitted, y = Residuals)) +
    geom_point(color = "#B163FF", alpha = 0.6) +
    geom_smooth(method = "loess", color = "#FF63D3", se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(x = "Fitted Values", y = "Scaled Residuals")

  plot_resid_qq <-  ggplot(
    data.frame(
      z.resid = scale(resid(lmm)),
      fitted = fitted(lmm)
    ),
    aes(sample = z.resid)
  ) +
    stat_qq() +
    geom_abline(
      intercept = 0,
      slope = 1,
      col = "red"
    ) +
    labs(x = "Theoretical quantiles", y = "Sample Quantiles") +
    theme_minimal()

  ranef_df <- broom.mixed::tidy(lmm, effects = "ran_vals", conf.int = TRUE)

  plot_random <- ggplot(ranef_df, aes(x = estimate, y = level)) +
    geom_point(color = "#B163FF") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", height = 0.2) +
    facet_wrap(~ term, scales = "free_x") +
    theme_minimal() +
    labs(x = "Estimate", y = metadata$grouping_var)

  fixed_plot <- ggplot(fixed, aes(x = estimate, y = term)) +
    geom_point(color = "black") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    theme_minimal() +
    labs(x = "Estimate", y = "Term")

  plot_pairs <- ggpairs(metadata$model_frame[, metadata$fixed_vars, drop = FALSE]) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 6)
    )

  infl <- influence(influence_model, group = metadata$grouping_var)

  cooks <- cooks.distance(infl)
  cooks_df <- data.frame(
    Group = rownames(cooks),
    CookD = as.numeric(cooks)
  )
  cutoff <- 4 / length(cooks)

  plot_influence <- ggplot(cooks_df, aes(x = CookD, y = Group)) +
    geom_point(aes(color = CookD > cutoff), size = 3) +
    scale_color_manual(values = c("FALSE" = "#B163FF", "TRUE" = "red")) +
    geom_vline(xintercept = cutoff, linetype = "dashed", color = "red") +
    geom_text_repel(
      data = subset(cooks_df, CookD > cutoff),
      aes(label = Group),
      color = "black",
      nudge_x = 0.01
    ) +
    labs(
      title = "Cook's Distance by Cluster",
      x = "Cook's Distance",
      y = paste0("Cluster (", metadata$grouping_var, ")")
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_flip()

  list(
    plot_resid_fitted = plot_resid_fitted,
    plot_resid_qq = plot_resid_qq,
    plot_random = plot_random,
    fixed_plot = fixed_plot,
    plot_pairs = plot_pairs,
    plot_influence = plot_influence,
    plots_effects = .build_diagnose_lmm_effect_plots(lmm, metadata$response_var)
  )
}

.build_diagnose_lmm_ui <- function(overfitting) {
  fluidPage(

    titlePanel("Linear Mixed Model Evaluation"),

    fluidRow(
      column(12, h3("Command")),
      column(12, verbatimTextOutput("model_command")),
    ),

    fluidRow(
      column(12, h3("Model Fit Statistics"), DTOutput("fit_stats")),
      if (overfitting) {
        column(12, strong("Warning: The model may be overfitting (singular fit)."), style = "color: red;")
      }
    ),

    fluidRow(
      column(6, h3("Fixed Effects"), DTOutput("fixed_effects")),
      column(6, h3("Random Effects"), DTOutput("random_effects"))
    ),

    fluidRow(
      column(6, h3("Fixed Effects Coefficients"), plotOutput("plot_fixed")),
      column(6, h3("Random Effects Caterpillar"), plotOutput("plot_random"))
    ),

    fluidRow(
      column(6, h3("Residuals vs Fitted"), plotOutput("plot_resid_fitted")),
      column(6, h3("Residuals Q-Q Plot"), plotOutput("plot_resid_qq"))
    ),


    fluidRow(
      column(6, h3("Cook's D"), plotOutput("plot_influence")),
      column(6, h3("Variable Correlations"), plotOutput("plot_pairs")),
    ),

    fluidRow(
      column(12, h3("Predictor Effects"), plotOutput("all_effects"))
    )
  )
}

.build_diagnose_lmm_server <- function(model_command, tables, plots) {
  function(input, output, session) {

    output$model_command <- renderText({ model_command })

    if ("p.value" %in% names(tables$fixed_df)) {
      output$fixed_effects <- renderDT({
        datatable(tables$fixed_df, options = list(pageLength = 5)) %>%
          formatStyle(
            "p.value",
            backgroundColor = styleInterval(
              tables$thresholds$p_value_highlight,
              c("lightgreen", "")
            )
          )
      })
    } else {
      output$fixed_effects <- renderDT({
        datatable(tables$fixed_df, options = list(pageLength = 5))
      })
    }

    output$random_effects <- renderDT({
      datatable(tables$rand, options = list(pageLength = 5)) %>%
        formatStyle(
          "Std_Dev",
          backgroundColor = styleInterval(
            tables$thresholds$residual_sd_highlight,
            c("", "lightcoral")
          )
        )
    })

    output$fit_stats <- DT::renderDataTable({ tables$fit_stats_df })

    output$plot_resid_fitted <- renderPlot({ plots$plot_resid_fitted })
    output$plot_resid_qq     <- renderPlot({ plots$plot_resid_qq })
    output$plot_random       <- renderPlot({ plots$plot_random })
    output$plot_fixed        <- renderPlot({ plots$fixed_plot })
    output$plot_pairs        <- renderPlot({ plots$plot_pairs })
    output$plot_influence    <- renderPlot({ plots$plot_influence })
    output$all_effects <- renderPlot({
      wrap_plots(
        plots$plots_effects,
        ncol = length(plots$plots_effects)
      )
    })
  }
}

# Make function for dashboard
diagnose_lmm <- function(lmm) {
  .validate_diagnose_lmm_input(lmm)

  library(shiny)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(performance)
  library(DT)
  library(dplyr)
  library(ggplot2)
  library(GGally)
  library(influence.ME)
  library(ggrepel)
  library(patchwork)
  library(effects)

  metadata <- .get_diagnose_lmm_metadata(lmm)
  thresholds <- .get_diagnose_lmm_thresholds()
  tables <- .build_diagnose_lmm_tables(lmm, metadata, thresholds)
  plots <- .build_diagnose_lmm_plots(lmm, metadata, tables$fixed)
  ui <- .build_diagnose_lmm_ui(tables$overfitting)
  server <- .build_diagnose_lmm_server(metadata$model_command, tables, plots)

  shinyApp(ui, server)
}
