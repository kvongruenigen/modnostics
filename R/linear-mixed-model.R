#' Diagnose a Linear Mixed Model
#'
#' This function launches a Shiny dashboard to inspect diagnostics of an LMM.
#'
#' @param lmm A fitted `lmer` model.
#' @return A Shiny app object.
#' @export

# Make function for dashboard
diagnose_lmm <- function(lmm) {

  # Load required packages
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

  ##############################################################################

  # Create a string for the model command
  model_command <- deparse(lmm@call) # "lmer(formula = weight ~ Time + Diet + (1 | Chick), data = cw)"

  # Create a string for the model formula
  model_formula <- as.character(lmm@call)[2] # "weight ~ Time + Diet + (1 | Chick)"

  # Make a list of the terms used in the model
  model_terms <- all.vars(as.formula(model_formula)) # [1] "weight" "Time"   "Diet"   "Chick"

  # Make a list of the data classed of the data used in the model
  data_classes <- sapply(eval(getCall(lmm)$data), class) # $weight [1] "numeric" $Time [1] "numeric" $Chick [1] "ordered" "factor" $Diet [1] "factor"

  ##############################################################################

  # Thresholds
  p_threshold_norm <- 0.05
  p_threshold_very <- 0.01
  p_threshold_highly <- 0.001
  resid_sd_threshold <- 10

  ##############################################################################

  # Extract model stats

  # Create a copy that can safely provide p-values
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

  ## Fixed effects (from lmerTest model)
  fixed <- broom.mixed::tidy(lmm_lmerTest, effects = "fixed", conf.int = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  ## Curate dataframe
  fixed_df <- fixed %>%
    mutate(conf_int = paste0("[", conf.low, ", ", conf.high, "]"))

  # Include p.value only if it exists
  if ("p.value" %in% names(fixed_df)) {
    fixed_df <- fixed_df %>% select(term, estimate, conf_int, p.value)
  } else {
    fixed_df <- fixed_df %>% select(term, estimate, conf_int)
  }

  ## Random effects
  rand <- as.data.frame(VarCorr(lmm))[,c("vcov", "sdcor")] %>%
    mutate(
      group = c(paste0(model_terms[length(model_terms)], " (Intercept)"), "Residual"),
      Variance = round(vcov, 3),
      Std_Dev = round(sdcor, 3)
    ) %>%
    select(group, Variance, Std_Dev)

  ## Model fit statistics
  fit_stats <- performance::model_performance(lmm) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
  # Convert to dataframe for DT
  fit_stats_df <- as.data.frame(fit_stats) %>%
    dplyr::mutate(across(where(is.numeric), ~ round(.x, 2)))

  # Residual SD flag
  resid_sd <- sigma(lmm)
  resid_flag <- ifelse(resid_sd > resid_sd_threshold, "⚠", "")

  overfitting <- lme4::isSingular(lmm, tol = 1e-4)

  ##############################################################################

  # Plots

  ## Residuals vs Fitted Plot}
  plot_resid_fitted <- ggplot(data.frame(
    Fitted = fitted(lmm),
    Residuals = scale(resid(lmm))
  ), aes(x = Fitted, y = Residuals)) +
    geom_point(color = "#B163FF", alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(x = "Fitted Values", y = "Scaled Residuals")


  ## Q-Q Plot of Residuals}
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

  ## Random Effects Caterpillar Plot}
  ranef_df <- broom.mixed::tidy(lmm, effects = "ran_vals", conf.int = TRUE)

  plot_random <- ggplot(ranef_df, aes(x = estimate, y = level)) +
    geom_point(color = "#CCCCFF") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", height = 0.2) +
    facet_wrap(~ term, scales = "free_x") +
    theme_minimal() +
    labs(x = "Estimate", y = model_terms[length(model_terms)])


  ## Fixed Effects Plot}
  fixed_plot <- ggplot(fixed, aes(x = estimate, y = term)) +
    geom_point(color = "black") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    theme_minimal() +
    labs(x = "Estimate", y = "Term")

  ## Variable Correlations
  # Create pairwise plots for all variables in the model except the random effect
  plot_pairs <- ggpairs(lmm@frame[, model_terms[-length(model_terms)]]) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 6)
    )

  ## Influence Index Plot
  # plot_influence <- recordPlot({
  #   influenceIndexPlot(influence(lmm), groups = model_terms[length(model_terms)], grid = FALSE, main = "")
  # })
  # Influence analysis at cluster level (Chick)
  infl <- influence(lmm, group = model_terms[length(model_terms)])

  # Extract Cook's distance
  cooks <- cooks.distance(infl)

  # Make into a dataframe for ggplot
  cooks_df <- data.frame(
    Group = rownames(cooks),
    CookD = as.numeric(cooks)
  )

  # Define cutoff (common rule of thumb: 4/nclusters)
  cutoff <- 4 / length(cooks)

  plot_influence <- ggplot(cooks_df, aes(x = CookD, y = Group)) +
    # Base points: flagged (red) vs safe (steelblue)
    geom_point(aes(color = CookD > cutoff), size = 3) +
    scale_color_manual(values = c("FALSE" = "#B163FF", "TRUE" = "red")) +
    # Cutoff line
    geom_vline(xintercept = cutoff, linetype = "dashed", color = "red") +
    # Labels for flagged clusters
    geom_text_repel(
      data = subset(cooks_df, CookD > cutoff),
      aes(label = Group),
      color = "black",
      nudge_x = 0.01
    ) +
    labs(
      title = "Cook's Distance by Cluster",
      x = "Cook's Distance",
      y = paste0("Cluster (", model_terms[length(model_terms)], ")")
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_flip()

  gg_predictor_effects <- function(model) {
    eff <- predictorEffects(model)

    plots <- lapply(names(eff), function(var) {
      df <- as.data.frame(eff[[var]])

      if (is.numeric(df[[var]])) {

        # Continuous predictor → line + ribbon
        p <- ggplot(df, aes_string(x = var, y = "fit")) +
          geom_line(color = "#B163FF", linewidth = 1) +
          geom_ribbon(aes(ymin = lower, ymax = upper),
                      alpha = 0.2, fill = "#CCCCFF") +
          labs(
            title = paste("Effect of", var),
            x = var,
            y = paste0("Predicted ", model_terms[1])
          )
      } else {

        # Categorical predictor → points + error bars
        p <- ggplot(df, aes_string(x = var, y = "fit")) +
          geom_point(size = 3, color = "#B163FF") +
          geom_errorbar(aes(ymin = lower, ymax = upper),
                        width = 0.1, color = "") +
          labs(
            title = paste("Effect of", var),
            x = var,
            y = paste0("Predicted ", model_terms[1])
          )
      }

      p + theme_minimal()
    })
    names(plots) <- names(eff)  # 👈 ensures list is keyed by predictor names
    plots
  }
  plots_effects <- gg_predictor_effects(lmm)
  ##############################################################################

  # UI
  ui <- fluidPage(

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

  # Server
  server <- function(input, output, session) {

    output$model_command <- renderText({ model_command })

    if ("p.value" %in% names(fixed_df)) {
      output$fixed_effects <- renderDT({
        datatable(fixed_df, options = list(pageLength = 5)) %>%
          formatStyle("p.value",
                      backgroundColor = styleInterval(p_threshold_norm, c("lightgreen", "")))
      })
    } else {
      output$fixed_effects <- renderDT({
        datatable(fixed_df, options = list(pageLength = 5))
      })
    }


    output$random_effects <- renderDT({
      datatable(rand, options = list(pageLength = 5)) %>%
        formatStyle("Std_Dev",
                    backgroundColor = styleInterval(resid_sd_threshold, c("", "lightcoral")))
    })

    output$fit_stats <- DT::renderDataTable({fit_stats_df})

    # Plots
    output$plot_resid_fitted <- renderPlot({ plot_resid_fitted })
    output$plot_resid_qq     <- renderPlot({ plot_resid_qq })
    output$plot_random       <- renderPlot({ plot_random })
    output$plot_fixed        <- renderPlot({ fixed_plot })
    output$plot_pairs        <- renderPlot({ plot_pairs })
    output$plot_influence    <- renderPlot({ plot_influence })
    output$all_effects <- renderPlot({
      wrap_plots(plots_effects, ncol = length(plots_effects))  # from patchwork
    })
  }

  shinyApp(ui, server)
}
