# plotting_functions.R
# Load necessary libraries for plotting
library(dplyr)
library(ggplot2)
library(gridExtra)

residualplot_step1 <- function(model, nameofoutcome, nameofage, data, 
                               ylimit = c(-50, 50), save_path = NULL, model_name = NULL) {
  require(dplyr)
  require(ggplot2)
  require(gridExtra)
  
  if (is.null(model_name)) {
    stop("Model name must be provided.")
  }
  
  # Check if the model has ng (number of classes), handle one-class model separately
  if (is.null(model$ng)) {
    # For one-class model, we can set ng to 1 (or handle it differently if needed)
    k <- 1
    preds <- model$pred
    names(preds)[6] <- nameofoutcome
    nameofid <- names(model$pred)[1]
    test <- dplyr::left_join(preds, model$pprob, by = nameofid)
    test <- dplyr::left_join(test, data, by = c(nameofid, nameofoutcome))
    
    # No need for class-based plotting in one-class model
    newplotvalues <- test %>% 
      mutate(Residuals = get(nameofoutcome) - pred_ss1)
    
    plotvaluessub <- newplotvalues
    p <- ggplot2::ggplot(data = plotvaluessub, aes(x = get(nameofage), y = Residuals)) +
      theme(axis.text = element_text(size = 8), text = element_text(size = 8)) + 
      geom_point(size = 0.1) + 
      stat_summary(fun = mean, geom = "line", size = 1, col = "CadetBlue", group = 1) + 
      ggtitle("Residuals in one class") + 
      ylim(ylimit) + 
      xlim(range(data[[nameofage]], na.rm = TRUE)) + 
      labs(x = "Year") + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
      )
    
    # Return the plot for one-class model
    return(p)
  }
  
  # Otherwise, proceed with the multi-class model code
  k <- model$ng
  preds <- model$pred
  names(preds)[6] <- nameofoutcome
  nameofid <- names(model$pred)[1]
  test <- dplyr::left_join(preds, model$pprob, by = nameofid)
  test <- dplyr::left_join(test, data, by = c(nameofid, nameofoutcome))
  
  plot_list <- list()
  xlim_range <- range(data[[nameofage]], na.rm = TRUE)
  
  for (i in 1:k) {
    newplotvalues <- test %>% 
      filter(class == i) %>% 
      mutate(Residuals = get(nameofoutcome) - eval(parse(text = paste0("pred_ss", i))))
    
    plotvaluessub <- newplotvalues
    p <- ggplot2::ggplot(data = plotvaluessub, aes(x = get(nameofage), y = Residuals, group = class)) +
      theme(axis.text = element_text(size = 8), text = element_text(size = 8)) + 
      geom_point(size = 0.1) + 
      stat_summary(fun = mean, geom = "line", size = 1, col = "CadetBlue", group = 1) + 
      ggtitle(paste("Residuals in class", i)) + 
      ylim(ylimit) + 
      xlim(xlim_range) + 
      labs(x = "Year") + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
      )
    
    plot_list[[i]] <- p
  }
  
  # Handle empty space if necessary
  total_plots <- length(plot_list)
  max_cols <- 7
  if (total_plots %% max_cols != 0) {
    empty_plots <- max_cols - (total_plots %% max_cols)
    for (j in seq_len(empty_plots)) {
      plot_list[[total_plots + j]] <- ggplot() + 
        theme_void() + 
        theme(panel.border = element_rect(colour = "white"))
    }
  }
  
  # Combine plots into a grid
  combined_plot <- grid.arrange(grobs = plot_list, ncol = max_cols, nrow = ceiling(total_plots / max_cols))
  
  if (!is.null(save_path)) {
    g <- arrangeGrob(grobs = plot_list, ncol = max_cols, nrow = ceiling(total_plots / max_cols))
    ggsave(filename = file.path(save_path, paste0(model_name, ".jpeg")), plot = g,
           width = 21, height = 3, dpi = 1200)
  }
  
  return(combined_plot)
}

