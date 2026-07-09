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

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "CookD", "Fitted", "Group", "Residuals", "Std_Dev", "Variance",
      "Tolerance", "VIF", "above_cutoff", "conf.high", "conf.low",
      "conf_int", "density", "estimate", "fit", "group", "level",
      "lower", "p.value", "predictor_value", "rank", "sdcor", "term",
      "upper", "vcov", "z.resid"
    )
  )
}

# Extract model metadata using stable accessors rather than direct slot parsing.
.get_diagnose_lmm_metadata <- function(lmm) {
  model_formula <- stats::formula(lmm)
  fixed_formula <- suppressWarnings(lme4::nobars(model_formula))
  model_frame <- stats::model.frame(lmm)
  random_groups <- names(lme4::getME(lmm, "flist"))
  response_var <- all.vars(fixed_formula)[1]
  fixed_vars <- setdiff(all.vars(fixed_formula), response_var)

  list(
    model_command = paste(deparse(stats::getCall(lmm)), collapse = "\n"),
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

.with_diagnose_lmm_data_update <- function(model, code) {
  had_data_update <- exists("data.update", envir = .GlobalEnv, inherits = FALSE)

  if (had_data_update) {
    old_data_update <- get("data.update", envir = .GlobalEnv, inherits = FALSE)
  }

  assign("data.update", stats::model.frame(model), envir = .GlobalEnv)

  on.exit(
    {
      if (had_data_update) {
        assign("data.update", old_data_update, envir = .GlobalEnv)
      } else if (exists("data.update", envir = .GlobalEnv, inherits = FALSE)) {
        rm("data.update", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )

  force(code)
}

.run_diagnose_lmm_influence <- function(model, grouping_var) {
  if (!("package:nlme" %in% search())) {
    base::attachNamespace(asNamespace("nlme"))
  }

  .with_diagnose_lmm_data_update(
    model,
    influence.ME::influence(model, group = grouping_var)
  )
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
    residual_sd_highlight = 10,
    vif_moderate = 5,
    vif_high = 10
  )
}

.build_diagnose_lmm_collinearity_table <- function(lmm, thresholds) {
  collinearity <- suppressMessages(
    suppressWarnings(
      try(as.data.frame(performance::check_collinearity(lmm)), silent = TRUE)
    )
  )

  if (inherits(collinearity, "try-error") || nrow(collinearity) == 0) {
    return(data.frame(
      Term = "Not available",
      VIF = NA_real_,
      Tolerance = NA_real_,
      Status = "Need at least two fixed-effect predictors for VIF."
    ))
  }

  collinearity_df <- data.frame(
    Term = collinearity[["Term"]],
    VIF = round(collinearity[["VIF"]], 2),
    Tolerance = round(collinearity[["Tolerance"]], 2)
  )

  dplyr::mutate(
    collinearity_df,
    Status = dplyr::case_when(
      VIF >= thresholds$vif_high ~ "High",
      VIF >= thresholds$vif_moderate ~ "Moderate",
      TRUE ~ "Low"
    )
  )
}

.build_diagnose_lmm_influence_summary <- function(lmm, metadata) {
  legacy_data_model <- .build_diagnose_lmm_influence_model(lmm)
  infl <- .run_diagnose_lmm_influence(legacy_data_model, metadata$grouping_var)

  cooks <- stats::cooks.distance(infl)
  cooks_df <- dplyr::arrange(
    dplyr::mutate(
      data.frame(
        Group = rownames(cooks),
        CookD = as.numeric(cooks)
      ),
      CookD = round(CookD, 3)
    ),
    dplyr::desc(CookD)
  )
  cutoff <- 4 / max(nrow(cooks_df), 1)

  dplyr::mutate(
    cooks_df,
    above_cutoff = CookD > cutoff,
    rank = dplyr::row_number(),
    cutoff = cutoff,
    Flag = ifelse(above_cutoff, "Review", "OK")
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
  influence_df <- .build_diagnose_lmm_influence_summary(lmm, metadata)
  collinearity_df <- .build_diagnose_lmm_collinearity_table(lmm, thresholds)

  fixed <- dplyr::mutate(
    broom.mixed::tidy(lmm_lmerTest, effects = "fixed", conf.int = TRUE),
    dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))
  )

  fixed_df <- dplyr::mutate(
    fixed,
    conf_int = paste0("[", conf.low, ", ", conf.high, "]")
  )

  if ("p.value" %in% names(fixed_df)) {
    fixed_df <- dplyr::select(fixed_df, term, estimate, conf_int, p.value)
  } else {
    fixed_df <- dplyr::select(fixed_df, term, estimate, conf_int)
  }

  rand <- dplyr::select(
    dplyr::mutate(
      as.data.frame(lme4::VarCorr(lmm))[, c("vcov", "sdcor")],
      group = c(paste0(metadata$grouping_var, " (Intercept)"), "Residual"),
      Variance = round(vcov, 3),
      Std_Dev = round(sdcor, 3)
    ),
    group,
    Variance,
    Std_Dev
  )

  fit_stats <- dplyr::mutate(
    performance::model_performance(lmm),
    dplyr::across(dplyr::where(is.numeric), ~ round(.x, 2))
  )
  fit_stats_df <- dplyr::mutate(
    as.data.frame(fit_stats),
    dplyr::across(dplyr::where(is.numeric), ~ round(.x, 2))
  )

  overfitting <- lme4::isSingular(lmm, tol = 1e-4)

  list(
    fixed = fixed,
    fixed_df = fixed_df,
    rand = rand,
    fit_stats_df = fit_stats_df,
    collinearity_df = collinearity_df,
    influence_df = influence_df,
    overfitting = overfitting,
    thresholds = thresholds
  )
}

.build_diagnose_lmm_effect_plots <- function(model, response_var) {
  eff <- .with_diagnose_lmm_data_update(
    model,
    effects::predictorEffects(model)
  )

  plots <- lapply(names(eff), function(var) {
    df <- as.data.frame(eff[[var]])
    df$predictor_value <- df[[var]]

    if (is.numeric(df[[var]])) {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = predictor_value, y = fit)
      ) +
        ggplot2::geom_line(color = "#B163FF", linewidth = 1) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower, ymax = upper),
          alpha = 0.2,
          fill = "#CCCCFF"
        ) +
        ggplot2::labs(
          title = paste("Effect of", var),
          x = var,
          y = paste0("Predicted ", response_var)
        )
    } else {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = predictor_value, y = fit)
      ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = lower, ymax = upper),
          width = 0.1,
          color = "#FF63D3"
        ) +
        ggplot2::geom_point(size = 3, color = "#B163FF") +
        ggplot2::labs(
          title = paste("Effect of", var),
          x = var,
          y = paste0("Predicted ", response_var)
        )
    }

    p + ggplot2::theme_minimal()
  })

  names(plots) <- names(eff)
  plots
}

.build_diagnose_lmm_plots <- function(lmm, metadata, fixed, influence_df) {
  legacy_data_model <- .build_diagnose_lmm_influence_model(lmm)
  residual_df <- data.frame(
    Fitted = stats::fitted(lmm),
    Residuals = as.numeric(scale(stats::resid(lmm)))
  )

  plot_resid_fitted <- ggplot2::ggplot(
    residual_df,
    ggplot2::aes(x = Fitted, y = Residuals)
  ) +
    ggplot2::geom_point(color = "#B163FF", alpha = 0.6) +
    ggplot2::geom_smooth(method = "loess", color = "#FF63D3", se = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Fitted Values", y = "Scaled Residuals")

  plot_resid_qq <-  ggplot2::ggplot(
    residual_df,
    ggplot2::aes(sample = Residuals)
  ) +
    ggplot2::stat_qq() +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      col = "red"
    ) +
    ggplot2::labs(x = "Theoretical quantiles", y = "Sample Quantiles") +
    ggplot2::theme_minimal()

  plot_resid_distribution <- ggplot2::ggplot(
    residual_df,
    ggplot2::aes(x = Residuals)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 20,
      fill = "#CCCCFF",
      color = "white"
    ) +
    ggplot2::geom_density(color = "#B163FF", linewidth = 1) +
    ggplot2::geom_vline(
      xintercept = c(-2, 0, 2),
      linetype = c("dotted", "dashed", "dotted"),
      color = c("#999999", "red", "#999999")
    ) +
    ggplot2::labs(
      x = "Scaled Residuals",
      y = "Density"
    ) +
    ggplot2::theme_minimal()

  ranef_df <- broom.mixed::tidy(lmm, effects = "ran_vals", conf.int = TRUE)

  plot_random <- ggplot2::ggplot(ranef_df, ggplot2::aes(x = estimate, y = level)) +
    ggplot2::geom_point(color = "#B163FF") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      height = 0.2
    ) +
    ggplot2::facet_wrap(~ term, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Estimate", y = metadata$grouping_var)

  fixed_plot <- ggplot2::ggplot(fixed, ggplot2::aes(x = estimate, y = term)) +
    ggplot2::geom_point(color = "black") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      height = 0.2
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Estimate", y = "Term")

  plot_pairs <- GGally::ggpairs(metadata$model_frame[, metadata$fixed_vars, drop = FALSE]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(size = 6)
    )
  cutoff <- unique(influence_df$cutoff)[1]

  plot_influence <- ggplot2::ggplot(
    influence_df,
    ggplot2::aes(x = CookD, y = stats::reorder(Group, CookD))
  ) +
    ggplot2::geom_point(ggplot2::aes(color = above_cutoff), size = 3) +
    ggplot2::scale_color_manual(values = c("FALSE" = "#B163FF", "TRUE" = "red")) +
    ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed", color = "red") +
    ggrepel::geom_text_repel(
      data = base::subset(influence_df, above_cutoff),
      ggplot2::aes(label = Group),
      color = "black",
      nudge_x = 0.01
    ) +
    ggplot2::labs(
      title = "Cook's Distance by Cluster",
      x = "Cook's Distance",
      y = paste0("Cluster (", metadata$grouping_var, ")")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::coord_flip()

  list(
    plot_resid_fitted = plot_resid_fitted,
    plot_resid_qq = plot_resid_qq,
    plot_resid_distribution = plot_resid_distribution,
    plot_random = plot_random,
    fixed_plot = fixed_plot,
    plot_pairs = plot_pairs,
    plot_influence = plot_influence,
    plots_effects = .build_diagnose_lmm_effect_plots(
      legacy_data_model,
      metadata$response_var
    )
  )
}

.build_diagnose_lmm_ui <- function(overfitting) {
  guidance_text <- .get_diagnose_lmm_guidance(overfitting)

  shiny::fluidPage(

    shiny::titlePanel("Linear Mixed Model Evaluation"),

    shiny::fluidRow(
      shiny::column(
        12,
        shiny::tags$div(
          style = "margin-bottom: 18px; padding: 12px; background-color: #f7f7f7; border-left: 4px solid #4a6fa5;",
          shiny::h4("How To Read This Dashboard"),
          shiny::p(guidance_text$overview)
        )
      )
    ),

    shiny::fluidRow(
      shiny::column(12, shiny::h3("Command")),
      shiny::column(12, shiny::verbatimTextOutput("model_command")),
    ),

    shiny::fluidRow(
      shiny::column(12, shiny::h3("Model Fit Statistics"), DT::DTOutput("fit_stats")),
      shiny::column(
        12,
        shiny::p(
          style = "color: #555;",
          guidance_text$model_fit
        )
      ),
      if (overfitting) {
        shiny::column(
          12,
          shiny::strong("Warning: The model may be overfitting (singular fit)."),
          style = "color: red;"
        )
      }
    ),

    shiny::fluidRow(
      shiny::column(
        6,
        shiny::h3("Fixed Effects"),
        shiny::p(style = "color: #555;", guidance_text$fixed_effects),
        DT::DTOutput("fixed_effects")
      ),
      shiny::column(
        6,
        shiny::h3("Random Effects"),
        shiny::p(style = "color: #555;", guidance_text$random_effects),
        DT::DTOutput("random_effects")
      )
    ),

    shiny::fluidRow(
      shiny::column(
        6,
        shiny::h3("Fixed Effects Coefficients"),
        shiny::p(style = "color: #555;", guidance_text$fixed_plot),
        shiny::plotOutput("plot_fixed")
      ),
      shiny::column(
        6,
        shiny::h3("Random Effects Caterpillar"),
        shiny::p(style = "color: #555;", guidance_text$random_plot),
        shiny::plotOutput("plot_random")
      )
    ),

    shiny::fluidRow(
      shiny::column(
        4,
        shiny::h3("Residuals vs Fitted"),
        shiny::p(style = "color: #555;", guidance_text$residual_fitted),
        shiny::plotOutput("plot_resid_fitted")
      ),
      shiny::column(
        4,
        shiny::h3("Residuals Q-Q Plot"),
        shiny::p(style = "color: #555;", guidance_text$residual_qq),
        shiny::plotOutput("plot_resid_qq")
      ),
      shiny::column(
        4,
        shiny::h3("Residual Distribution"),
        shiny::p(style = "color: #555;", guidance_text$residual_distribution),
        shiny::plotOutput("plot_resid_distribution")
      )
    ),


    shiny::fluidRow(
      shiny::column(
        8,
        shiny::h3("Cook's D"),
        shiny::p(style = "color: #555;", guidance_text$influence),
        shiny::plotOutput("plot_influence")
      ),
      shiny::column(
        4,
        shiny::h3("Influence Summary"),
        shiny::p(style = "color: #555;", guidance_text$influence_table),
        DT::DTOutput("influence_summary")
      )
    ),

    shiny::fluidRow(
      shiny::column(
        6,
        shiny::h3("Variable Correlations"),
        shiny::p(style = "color: #555;", guidance_text$correlations),
        shiny::plotOutput("plot_pairs")
      ),
      shiny::column(
        6,
        shiny::h3("Collinearity (VIF)"),
        shiny::p(style = "color: #555;", guidance_text$collinearity),
        DT::DTOutput("collinearity")
      ),
    ),

    shiny::fluidRow(
      shiny::column(
        12,
        shiny::h3("Predictor Effects"),
        shiny::p(style = "color: #555;", guidance_text$predictor_effects),
        shiny::plotOutput("all_effects")
      )
    )
  )
}

.get_diagnose_lmm_guidance <- function(overfitting) {
  singular_note <- if (overfitting) {
    paste(
      "This model was flagged as singular, so treat variance components and",
      "random-effects interpretation cautiously."
    )
  } else {
    "No singular-fit warning was detected, but model assumptions should still be checked visually."
  }

  list(
    overview = paste(
      "Use this dashboard as a structured review rather than a pass/fail checklist.",
      "Look for patterns that suggest follow-up modelling decisions, then confirm those",
      "decisions in the context of your study design and domain knowledge."
    ),
    model_fit = paste(
      "These summaries describe overall model fit and complexity.",
      singular_note
    ),
    fixed_effects = paste(
      "Review estimates, confidence intervals, and p-values together.",
      "Small p-values can be useful, but effect size, direction, and uncertainty often matter more."
    ),
    random_effects = paste(
      "Random-effects estimates help you judge how much variation is attributed to grouping structure.",
      "Very small variance estimates can indicate limited group-level signal."
    ),
    fixed_plot = paste(
      "This coefficient plot is a quick visual check of direction and uncertainty.",
      "Effects with intervals far from zero usually provide stronger evidence of a stable association."
    ),
    random_plot = paste(
      "Use the caterpillar plot to see how group-level deviations are distributed.",
      "Large spread suggests meaningful between-group variation."
    ),
    residual_fitted = paste(
      "A good residuals-vs-fitted plot usually looks patternless and centered around zero.",
      "Curvature or funnel shapes can suggest nonlinearity or non-constant variance."
    ),
    residual_qq = paste(
      "Points close to the reference line are more consistent with normal residuals.",
      "Systematic departures in the tails may suggest skewness, heavy tails, or influential observations."
    ),
    residual_distribution = paste(
      "This histogram and density curve make it easier to see skewness, heavy tails, or unusually wide spread.",
      "Large piles in the tails or strong asymmetry can support what you see in the Q-Q plot."
    ),
    influence = paste(
      "Clusters above the Cook's D cutoff deserve a second look.",
      "Influential cases are not automatically wrong, but they may justify sensitivity analyses."
    ),
    influence_table = paste(
      "Use this table to rank the most influential clusters and see which ones cross the heuristic cutoff.",
      "Start by reviewing flagged groups, then compare conclusions with and without them if needed."
    ),
    correlations = paste(
      "Strong predictor correlations can complicate interpretation and inflate uncertainty.",
      "If you see strong overlap, check whether collinearity is affecting the model."
    ),
    collinearity = paste(
      "Variance inflation factors quantify how strongly predictors overlap in the fitted model.",
      "Values around 5 deserve attention, and values around 10 usually warrant closer review."
    ),
    predictor_effects = paste(
      "These plots show how the fitted model translates predictors into expected outcomes.",
      "Use them to explain practical meaning, not just statistical significance."
    )
  )
}

.build_diagnose_lmm_server <- function(model_command, tables, plots) {
  function(input, output, session) {

    output$model_command <- shiny::renderText({ model_command })

    if ("p.value" %in% names(tables$fixed_df)) {
      output$fixed_effects <- DT::renderDT({
        DT::formatStyle(
          DT::datatable(tables$fixed_df, options = list(pageLength = 5)),
            "p.value",
            backgroundColor = DT::styleInterval(
              tables$thresholds$p_value_highlight,
              c("lightgreen", "")
            )
          )
      })
    } else {
      output$fixed_effects <- DT::renderDT({
        DT::datatable(tables$fixed_df, options = list(pageLength = 5))
      })
    }

    output$random_effects <- DT::renderDT({
      DT::formatStyle(
        DT::datatable(tables$rand, options = list(pageLength = 5)),
          "Std_Dev",
          backgroundColor = DT::styleInterval(
            tables$thresholds$residual_sd_highlight,
            c("", "lightcoral")
          )
        )
    })

    output$fit_stats <- DT::renderDataTable({ tables$fit_stats_df })
    output$collinearity <- DT::renderDT({
      collinearity_table <- DT::datatable(
        tables$collinearity_df,
        options = list(pageLength = 5, dom = "t"),
        rownames = FALSE
      )

      if ("VIF" %in% names(tables$collinearity_df) &&
          all(is.na(tables$collinearity_df$VIF))) {
        return(collinearity_table)
      }

      DT::formatStyle(
        collinearity_table,
        "VIF",
        backgroundColor = DT::styleInterval(
          c(
            tables$thresholds$vif_moderate,
            tables$thresholds$vif_high
          ),
          c("", "#fff3cd", "#f8d7da")
        )
      )
    })
    output$influence_summary <- DT::renderDT({
      DT::formatStyle(
        DT::datatable(
          tables$influence_df[, c("rank", "Group", "CookD", "Flag")],
          options = list(pageLength = 8, dom = "t"),
          rownames = FALSE
        ),
        "Flag",
        backgroundColor = DT::styleEqual(
          c("OK", "Review"),
          c("", "#f8d7da")
        )
      )
    })

    output$plot_resid_fitted <- shiny::renderPlot({ plots$plot_resid_fitted })
    output$plot_resid_qq     <- shiny::renderPlot({ plots$plot_resid_qq })
    output$plot_resid_distribution <- shiny::renderPlot({
      plots$plot_resid_distribution
    })
    output$plot_random       <- shiny::renderPlot({ plots$plot_random })
    output$plot_fixed        <- shiny::renderPlot({ plots$fixed_plot })
    output$plot_pairs        <- shiny::renderPlot({ plots$plot_pairs })
    output$plot_influence    <- shiny::renderPlot({ plots$plot_influence })
    output$all_effects <- shiny::renderPlot({
      patchwork::wrap_plots(
        plots$plots_effects,
        ncol = length(plots$plots_effects)
      )
    })
  }
}

# Make function for dashboard
diagnose_lmm <- function(lmm) {
  .validate_diagnose_lmm_input(lmm)

  metadata <- .get_diagnose_lmm_metadata(lmm)
  thresholds <- .get_diagnose_lmm_thresholds()
  tables <- .build_diagnose_lmm_tables(lmm, metadata, thresholds)
  plots <- .build_diagnose_lmm_plots(lmm, metadata, tables$fixed, tables$influence_df)
  ui <- .build_diagnose_lmm_ui(tables$overfitting)
  server <- .build_diagnose_lmm_server(metadata$model_command, tables, plots)

  shiny::shinyApp(ui, server)
}
