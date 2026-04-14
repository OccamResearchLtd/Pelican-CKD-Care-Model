# The function to generate a tornado plot
tornado_plot <- function(lower.mat, upper.mat, base, index) {
  # variables for plotting
  outcome <- switch(index, "Incremental cost", "Incremental QALY", "ICER")
  title_plot <- switch(index, "Incremental cost ($)", "Incremental QALY", "ICER ($/QALY)")
  format_hover <- switch(index, "$,.0f", ".2f", "$,.0f")
  format_tick <- switch(index, "$,.0f", ".1f", "$,.0f")
  base_case <- base[outcome]
  
  # Create a data frame
  df <- data.frame(
    Variable = c("Cost of CKD health states",
                 "Utility of CKD health states",
                 "RR HF - CV Meds",
                 "RR MI - CV Meds",
                 "RR Stroke - CV Meds",
                 "Annual price - CV Meds",
                 "RR HF - SGLT2i",
                 "RR MI - SGLT2i",
                 "RR Stroke - SGLT2i",
                 "eGFR effect - SGLT2i",
                 "Annual price - SGLT2i"),
    Lower = lower.mat[outcome,] - base_case,
    Upper = upper.mat[outcome,] - base_case,
    Difference = abs(lower.mat[outcome,] - upper.mat[outcome,])
  )
  
  # Order the variables by the absolute difference between the lower and upper results
  df <- df[order(df$Difference), ]
  
  # Set up the range of x axis
  minrange <- min(lower.mat[outcome,],upper.mat[outcome,])
  maxrange <- max(lower.mat[outcome,],upper.mat[outcome,])
  xrange <- c(minrange - 0.36*(maxrange - minrange), maxrange + 0.36*(maxrange - minrange))
  
  # Create the Tornado diagram
  fig <- plot_ly(df, x = ~Lower, y = ~Variable, type = 'bar', name = 'Lower value', base = ~base_case, marker = list(color = "#AFABAB"), hoverinfo = 'x+text', xhoverformat = format_hover,
                 hovertext = ~paste(outcome, "at lower value")) %>%
    add_trace(x = ~Upper, name = 'Upper value', base = ~base_case, marker = list(color = "#00B050"),
              hovertext = ~paste(outcome, "at upper value")) %>%
    layout(barmode = 'relative', 
           xaxis = list(title = title_plot, range = xrange, showticklabels = T, tickformat = format_tick, zeroline = FALSE), 
           yaxis = list(title = "", categoryorder = "trace"),
           shapes = list(list(type = "line", x0 = base_case, x1 = base_case, y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1))),
           bargap = 0.4) %>% 
    config(displayModeBar = F)
}