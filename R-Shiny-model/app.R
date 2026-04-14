############################
# CKD model with Shiny UI  #
# Programmed by Ziyi Lin   #
############################

# Use a loop to check whether multiple packages are installed in R before loading them ----
packages <- c("shiny", "plotly", "readr", "ggplot2", "rstudioapi", "htmlTable", "dplyr")
for (package in packages) {
  if (!requireNamespace(package, quietly = T)) {
    install.packages(package)
  }
}

# Load the necessary libraries
library(shiny)
library(plotly)
library(readr)
library(tidyverse)
library(rstudioapi)
library(htmlTable)

# Load the CKD model as a function ----
source("function/CKD.Excel.R")

# Load the global function for tornado plots
source("function/Tornado.R")

# Define the events of interest
states <- c("No event", "Non-fatal ESRD", "Non-fatal stroke", "Non-fatal MI", "Non-fatal HF", "CR death", "Other death")
states_factor <- factor(states, levels = rev(states))
components <- c("No event", "ESRD", "Stroke", "MI", "HF", "Medication")
components_q <- c("No event", "ESRD", "Stroke", "MI", "HF")

# Define colors for the charts
colors_pie <- c("#105671", "#9E480E")
colors_bar <- c("#70AD47", "#5B9BD5", "#C00000", "#ED7D31", "#A5A5A5", "#FFC000")

color_palette_com <- c("No event" = "#105671", "ESRD" = "#FFC000", "Stroke" = "#A5A5A5", "MI" = "#ED7D31", "HF" = "#C00000", "Medication" = "#7030A0")
color_palette_int <- c("No event" = "#8497B0", "ESRD" = "#FFD966", "Stroke" = "#C9C9C9", "MI" = "#F4B084", "HF" = "#E07F7F", "Medication" = "#B889DB")

color_palette2_com <- c("ESRD" = "#FFC000", "Stroke" = "#A5A5A5", "MI" = "#ED7D31", "HF" = "#C00000","CR death" = "#5B9BD5", "Other death" = "#70AD47")
color_palette2_int <- c("ESRD" = "#FFD966", "Stroke" = "#C9C9C9", "MI" = "#F4B084", "HF" = "#E07F7F","CR death" = "#9BC2E6", "Other death" = "#A9D08E")

# Define the size of the charts
h_num <- 330
w_num <- 583
h_px <- paste0(h_num,"px")
w_px <- paste0(w_num,"px")

# Define UI for application
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .container-fluid {
        min-width: 1460px;
      }
    "))
  ),
  titlePanel("CKD Care Model"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 # Input: Selector for choosing intervention ----
                 h3("Comparison"),
                 selectInput(inputId = "com",
                             label = "Choose a comparator:",
                             choices = c(
                               "Usual care",
                               "CV risk factor management", 
                               "CV risk factor management for all eligible", 
                               "CV risk factor management for all eligible and 50% CKD management",
                               "CV risk factor management for all eligible and full CKD management"),
                             selected = "Usual care"),
                 
                 selectInput(inputId = "int",
                             label = "Choose an intervention:",
                             choices = c(
                               "CV risk factor management", 
                               "CV risk factor management for all eligible", 
                               "CV risk factor management for all eligible and 50% CKD management",
                               "CV risk factor management for all eligible and full CKD management",
                               "CV risk factor management for all eligible and full CKD management, no delay"),
                             selected = "CV risk factor management for all eligible and full CKD management, no delay"),
                 
                 br(),
                 h3("Time horizon"),
                 selectInput(inputId = "timehorizon",
                             label = "Choose a time horizon:",
                             choices = c(
                               "1 year", 
                               "2 years",
                               "3 years", 
                               "4 years",
                               "5 years", 
                               "10 years",
                               "15 years",
                               "20 years",
                               "40 years",
                               "70 years (lifetime)"),
                             selected = "10 years"),
                 br(),
                 h3("Model population"),
                 selectInput(inputId = "egfrstage",
                             label = "eGFR stage (mL/min/1.73m^2)",
                             choices = c(
                               "G1 (≥90)",
                               "G2 (60-89)", 
                               "G3a (45-59)", 
                               "G3b (30-44)",
                               "G4 (15-29)",
                               "G5 (<15)"),
                             selected = "G3a (45-59)"),
                 
                 selectInput(inputId = "albstage",
                             label = "Albuminuria uACR stage (mg/g)",
                             choices = c(
                               "A1 (<30)",
                               "A2 (30-300)", 
                               "A3 (>300)"),
                             selected = "A1 (<30)"),
                 
                 sliderInput(inputId = "age",
                             label = "Age (year)",
                             min = 18,
                             max = 99,
                             value = 57,
                             step = 1),
                 
                 sliderInput(inputId = "dia",
                             label = "% Diabetes",
                             min = 0,
                             max = 100,
                             value = 44.5,
                             step = 0.1),
                 
                 sliderInput(inputId = "sex",
                             label = "% Female",
                             min = 0,
                             max = 100,
                             value = 40.9,
                             step = 0.1),
                 
                 sliderInput(inputId = "hyp",
                             label = "% Hypertensive (BP ≥140/90 mmHg or baseline BP treatment)",
                             min = 0,
                             max = 100,
                             value = 86.1,
                             step = 0.1),
                 br(),
                 
                 
                 
    ),
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Dashboard", 
                 br(),
                 br(),
                 fluidRow(column(width = 7,plotlyOutput("chart_com", height = h_px, width = w_px)),
                          column(width = 5,plotlyOutput("chart_int", height = h_px, width = w_px))),
                 br(),
                 tags$hr(style = "border-top: 1px solid #C9C9C9"),
                 fluidRow(column(width = 7,plotlyOutput("chart_cost_s", height = h_px, width = w_px)),
                          column(width = 5,plotlyOutput("chart_qaly_s", height = h_px, width = w_px))),
                 br(),
                 tags$hr(style = "border-top: 1px solid #C9C9C9"),
                 fluidRow(column(width = 7,plotlyOutput("chart_event", height = h_px, width = w_px)),
                          column(width = 5,plotlyOutput("chart_tte", height = h_px, width = w_px))),
                 br(),
                 tags$hr(style = "border-top: 1px solid #C9C9C9"),
                 fluidRow(column(width = 7,plotlyOutput("chart_costlt", height = h_px, width = w_px)),
                          column(width = 5,plotlyOutput("chart_qalylt", height = h_px, width = w_px))),
                 br(),
                 br()
        ),
        tabPanel("About",
                 div(style = "padding-left: 30px;",
                     br(),
                     h1("About the Model"),
                     div(style = "padding-left: 30px;",
                         h3("Model Schematic"),
                         img(src="Model_Schematic.png", height = 350, width = 750),
                         p("Utility: function of health state, eGFR, and diabetes"),
                         p("Cost: function of health state, eGFR, uACR, and diabetes"),
                         p("Post-event survival: function of age at event and eGFR"),
                         br(),
                         
                         h3("Intervention Schematic"),
                         img(src="Intervention_Schematic.png", height = 300, width = 700),
                         HTML("<p><span style='color:#1B8FBD;'>eGFR</span>: Function of baseline eGFR and SGLT2i medication use</p>"),
                         HTML("<p><span style='color:#00B050;'>RR(x)</span>: Statin relative risk for CVD event (x) per new statin user (adjusted for current statin utilization)</p>"),
                         HTML("<p><span style='color:#7030A0;'>RR(x)</span>: Blood pressure medication relative risk for CVD event (x) per additional medication (adjusted for current utilization)</p>"),
                         HTML("<p><span style='color:#FF0000;'>RR(x)</span>: Aspirin relative risk  for CVD event (x) per new aspirin user (adjusted for current aspirin utilization)</p>"),
                         HTML("<p><span style='color:#1B8FBD;'>RR(x)</span>: SGLT2i relative risk for CVD event (x) per new SGLT2i user</p>"),
                         
                         br()
                     )
                 )
                 

                 
        ),
        tabPanel("CEA result",
                 div(style = "padding-left: 30px;",
                     br(),
                     h1("Output of the Core Model"),
                     div(style = "padding-left: 30px;",
                         plotOutput(outputId = "distPlot1", height = 360, width = 700),
                         plotOutput(outputId = "distPlot2", height = 360, width = 700),
                         br()
                     ),
                     uiOutput("dynamic_title"),
                     div(style = "padding-left: 30px;",
                         uiOutput("CEA_s")
                         ),
                     br(),
                     br()
                 )
                 
        ),
        tabPanel("OWSA",
                 div(style = "padding-left: 40px;",
                     br(),
                     h1("One-Way Sensitivity Analysis"),
                     actionButton(inputId = "run_DSA", label = "Run One-way SA"),
                     h3("Tornado Diagram for Incremental Cost"),
                     fluidRow(plotlyOutput("DSAPlot1", height = h_px, width = "800px")),
                     br(),
                     h3("Tornado Diagram for Incremental QALY"),
                     fluidRow(plotlyOutput("DSAPlot2", height = h_px, width = "800px")),
                     br(),
                     h3("Tornado Diagram for ICER"),
                     fluidRow(plotlyOutput("DSAPlot3", height = h_px, width = "800px")),
                     br(), # for TSA
                     br()
                 )
        ),
        tabPanel("PSA",
                 div(style = "padding-left: 40px;",
                     br(),
                     h1("Probabilistic Sensitivity Analysis"),
                     actionButton(inputId = "run_PSA", label = "Run PSA"),
                     h3("Cost-Effectiveness Scatter Plot"),
                     fluidRow(plotlyOutput("PSAPlot1", height = h_px, width = "800px")),
                     br(),
                     h3("Cost-Effectiveness Acceptability Curve (CEAC)"),
                     fluidRow(plotlyOutput("PSAPlot2", height = h_px, width = "800px")),
                     br(),
                     h3("Expected Value of Perfect Information (EVPI)"),
                     fluidRow(plotlyOutput("PSAPlot3", height = h_px, width = "800px")),
                     br(), # for PSA
                     br()
                 )
        )
      )
    )
  )
)




# Define server logic
server <- function(input, output) {
  
  # For dashboard and base CEA results
    results <- reactive({
      timehorizon <- input$timehorizon
      egfrstage <- input$egfrstage
      albstage <- input$albstage
      age <- input$age
      dia <- input$dia
      sex <- input$sex
      hyp <- input$hyp
      int <- input$int
      com <- input$com
      update_sub <- T
      
      CKD_model(analysis = 1,
                input_timehorizon = timehorizon,
                input_update_sub = update_sub,
                input_egfrstage = egfrstage,
                input_albstage = albstage,
                input_age = age,
                input_dia = dia / 100,
                input_sex = sex / 100,
                input_hyp = hyp / 100,
                intervention = int, comparator = com)
    }) 
    
    output$dynamic_title <- renderUI({
      h1(paste0("Cost and Effectiveness Results at ", input$timehorizon))
    })
    
    output$distPlot1 <- renderPlot({
      
      results()$plot1
      
    })
    
    output$distPlot2 <- renderPlot({
      
      results()$plot2
      
    })
    
    output$CEA_s <- renderUI({
      results()$CEAtable_s
    })
    
    output$CEA <- renderUI({
      results()$CEAtable
    })
    
    output$chart_com <- renderPlotly({
      new_event_short_com <- results()$event_short_com
      new_event_short_int <- results()$event_short_int
      
      # Group "ESRD", "Stroke", "MI", "HF", "CR death", "Other death" into "Events"
      states_grouped <- c("No event", "Events")
      new_event_short_com_grouped <- c(new_event_short_com[1], sum(new_event_short_com[2:7]))
      new_event_short_int_grouped <- c(new_event_short_int[1], sum(new_event_short_int[2:7]))
      
      # Calculate the rotation angles
      rotation_com <- 90 + 0.5 * new_event_short_com_grouped[2] / sum(new_event_short_com) * 360
      rotation_int <- 90 + 0.5 * new_event_short_int_grouped[2] / sum(new_event_short_int) * 360
      
      # Define the maximum y-axis range for the column charts
      max_range <- min(max(sum(new_event_short_com[2:7]), sum(new_event_short_int[2:7])) * 1.2,1)
      
      fig <- plot_ly(domain = list(x = c(0, 0.7), y = c(0, 1)))
      fig <- fig %>% add_pie(data = data.frame(states_grouped, new_event_short_com_grouped), labels = ~states_grouped, values = ~new_event_short_com_grouped, name = "SOC", hole = 0.3, marker = list(colors = colors_pie), sort = FALSE, 
                             rotation = rotation_com, hoverinfo = 'label+percent+name', pull = c(0, 0.1))
      fig <- fig %>% layout(title = paste0("Events at ", input$timehorizon, ": SOC"), 
                            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), 
                            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      
      fig3 <- plot_ly(x = c('Event Breakdown'), y = new_event_short_com[2:7], type = 'bar', name = states_factor[2:7], 
                      marker = list(color = colors_bar), hoverinfo = 'y+name', yhoverformat = '.1%') %>%
        layout(yaxis = list(title = 'Percentage', range = c(0, max_range), showgrid = FALSE, zeroline = FALSE, showticklabels = T, tickformat = '.0%', side = 'right'), barmode = 'stack',
               legend = list(x = 1.1, y = 0.9))
      
      subplot(fig, fig3, nrows = 1, widths = c(0.8, 0.2)) %>% config(displayModeBar = F)
    })
    
    output$chart_int <- renderPlotly({
      new_event_short_com <- results()$event_short_com
      new_event_short_int <- results()$event_short_int
      
      # Group "ESRD", "Stroke", "MI", "HF", "CR death", "Other death" into "Events"
      states_grouped <- c("No event", "Events")
      new_event_short_com_grouped <- c(new_event_short_com[1], sum(new_event_short_com[2:7]))
      new_event_short_int_grouped <- c(new_event_short_int[1], sum(new_event_short_int[2:7]))
      
      # Calculate the rotation angles
      rotation_com <- 90 + 0.5 * new_event_short_com_grouped[2] / sum(new_event_short_com) * 360
      rotation_int <- 90 + 0.5 * new_event_short_int_grouped[2] / sum(new_event_short_int) * 360
      
      # Define the maximum y-axis range for the column charts
      max_range <- min(max(sum(new_event_short_com[2:7]), sum(new_event_short_int[2:7])) * 1.2,1)
      
      fig <- plot_ly(domain = list(x = c(0, 0.7), y = c(0, 1)))
      fig <- fig %>% add_pie(data = data.frame(states_grouped, new_event_short_int_grouped), labels = ~states_grouped, values = ~new_event_short_int_grouped, name = "intervention", hole = 0.3, marker = list(colors = colors_pie), sort = FALSE, 
                             rotation = rotation_int, hoverinfo = 'label+percent+name', pull = c(0, 0.1))
      fig <- fig %>% layout(title = paste0("Events at ", input$timehorizon, ": CVD optimized"),
                            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), 
                            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      
      fig3 <- plot_ly(x = c('Event Breakdown'), y = new_event_short_int[2:7], type = 'bar', name = states_factor[2:7], 
                      marker = list(color = colors_bar), hoverinfo = 'y+name', yhoverformat = '.1%') %>%
        layout(yaxis = list(title = 'Percentage', range = c(0, max_range), showgrid = FALSE, zeroline = FALSE, showticklabels = T, tickformat = '.0%', side = 'right'), barmode = 'stack',
               legend = list(x = 1.1, y = 0.9))
      
      subplot(fig, fig3, nrows = 1, widths = c(0.8, 0.2)) %>% config(displayModeBar = F)
    })
  
    output$chart_cost_s <- renderPlotly({
      df <- results()$cost_short
      
      # Add a new column to the data frame for colors
      df$color <- ifelse(df$category == "SOC", color_palette_com[df$component], color_palette_int[df$component])
      
      # Convert the 'component' column to a factor
      df$component <- factor(df$component, levels = components)
      
      # Create the plot
      fig <- plot_ly(df, x = ~cost, y = ~category, type = 'bar', orientation = 'h', marker = list(color = ~color), hoverinfo = 'y+text',
                     hovertext = ~paste("Cost component: ", component, "<br>Value: $", format(cost, digits = 2, big.mark = ","))) %>%
        layout(title = paste0("Disagregated Health Care Costs at ", input$timehorizon),
               barmode = 'stack', showlegend = F, yaxis = list(title = ""), 
               xaxis = list(title = 'Cost', showticklabels = T, tickformat = '$,.0f'),
               bargap = 0.4,margin = list(t = 80)) %>% config(displayModeBar = F)
      
    })
    
    output$chart_qaly_s <- renderPlotly({
      df <- results()$qaly_short
      
      # Add a new column to the data frame for colors
      df$color <- ifelse(df$category == "SOC", color_palette_com[df$component], color_palette_int[df$component])
      
      # Convert the 'component' column to a factor
      df$component <- factor(df$component, levels = components_q)
      
      # Create the plot
      fig <- plot_ly(df, x = ~qaly, y = ~category, type = 'bar', orientation = 'h', marker = list(color = ~color), hoverinfo = 'y+text',
                     hovertext = ~paste("QALY component: ", component, "<br>Value: ", format(qaly, digits = 3, big.mark = ","))) %>%
        layout(title = paste0("Disagregated QALYs at ", input$timehorizon),
               barmode = 'stack', showlegend = F, yaxis = list(title = ""), 
               xaxis = list(title = 'QALY', showticklabels = T),
               bargap = 0.4,margin = list(t = 80)) %>% config(displayModeBar = F)
      
    })
    
    output$chart_event <- renderPlotly({
      event_com <- results()$event_lt_com
      event_int <- results()$event_lt_int
      
      health_states <- c("Non-fatal ESRD", "Non-fatal stroke", "Non-fatal MI", "Non-fatal HF", "CR death", "Other death")
      
      # Create data frames for SOC and Intervention
      df_soc <- data.frame(State = health_states, Proportion = event_com, color = color_palette2_com)
      df_int <- data.frame(State = health_states, Proportion = event_int, color = color_palette2_int)
      
      # Create the plot
      p <- plot_ly() %>%
        add_trace(data = df_soc, x = ~Proportion, y = ~State, type = 'bar', orientation = 'h', marker = list(color = ~color), xhoverformat = '.1%', name = "SOC") %>%
        add_trace(data = df_int, x = ~Proportion, y = ~State, type = 'bar', orientation = 'h', marker = list(color = ~color), xhoverformat = '.1%', name = "Intervention") %>%
        layout(title = "Proportion of Events (lifetime)", barmode = 'group', showlegend = F, 
               yaxis = list(title = "", categoryorder = "trace"), 
               xaxis = list(title = "", showticklabels = T, tickformat = '.0%'),
               margin = list(t = 50)) %>% config(displayModeBar = F)

      
    })
    
    output$chart_tte <- renderPlotly({
      df_soc <- results()$tte_com
      df_int <- results()$tte_int
      
      df_soc$color <- color_palette2_com
      df_int$color <- color_palette2_int
      
      # Create the plot
      p <- plot_ly() %>%
        add_trace(data = df_soc, x = ~TimeToEvent, y = ~Event, type = 'bar', orientation = 'h', marker = list(color = ~color), xhoverformat = '.2f', name = "SOC") %>%
        add_trace(data = df_int, x = ~TimeToEvent, y = ~Event, type = 'bar', orientation = 'h', marker = list(color = ~color), xhoverformat = '.2f', name = "Intervention") %>%
        layout(title = "Time to Event (lifetime)", barmode = 'group', showlegend = F, 
               yaxis = list(title = "", categoryorder = "trace"), 
               xaxis = list(title = "", showticklabels = T, tickformat = '.1f'),
               margin = list(t = 50)) %>% config(displayModeBar = F)
      
    })
    
    output$chart_costlt <- renderPlotly({
      df <- results()$cost_lt
      
      # Add a new column to the data frame for colors
      df$color <- ifelse(df$category == "SOC", color_palette_com[df$component], color_palette_int[df$component])
      
      # Convert the 'component' column to a factor
      df$component <- factor(df$component, levels = components)
      
      # Create the plot
      fig <- plot_ly(df, x = ~cost, y = ~category, type = 'bar', orientation = 'h', marker = list(color = ~color), hoverinfo = 'y+text',
                     hovertext = ~paste("Cost component: ", component, "<br>Value: $", format(cost, digits = 2, big.mark = ","))) %>%
        layout(title = "Disagregated Health Care Costs (lifetime)",
               barmode = 'stack', showlegend = F, yaxis = list(title = ""), 
               xaxis = list(title = 'Cost', showticklabels = T, tickformat = '$,.0f'),
               bargap = 0.4,margin = list(t = 80)) %>% config(displayModeBar = F)
      
    })
    
    output$chart_qalylt <- renderPlotly({
      df <- results()$qaly_lt
      
      # Add a new column to the data frame for colors
      df$color <- ifelse(df$category == "SOC", color_palette_com[df$component], color_palette_int[df$component])
      
      # Convert the 'component' column to a factor
      df$component <- factor(df$component, levels = components_q)
      
      # Create the plot
      fig <- plot_ly(df, x = ~qaly, y = ~category, type = 'bar', orientation = 'h', marker = list(color = ~color), hoverinfo = 'y+text',
                     hovertext = ~paste("QALY component: ", component, "<br>Value: ", format(qaly, digits = 3, big.mark = ","))) %>%
        layout(title = "Disagregated QALYs (lifetime)",
               barmode = 'stack', showlegend = F, yaxis = list(title = ""), 
               xaxis = list(title = 'QALY', showticklabels = T),
               bargap = 0.4,margin = list(t = 80)) %>% config(displayModeBar = F)
      
    })
    
    # For DSA
    observeEvent(input$run_DSA, {
      # Show the "Please wait..." message
      showModal(modalDialog(
        title = "Please wait...",
        "The One-way Sensitivity Analysis is running. This can take a couple of minutes.",
        easyClose = FALSE,
        footer = NULL
      ))

      timehorizon <- input$timehorizon
      egfrstage <- input$egfrstage
      albstage <- input$albstage
      age <- input$age
      dia <- input$dia
      sex <- input$sex
      hyp <- input$hyp
      int <- input$int
      com <- input$com
      update_sub <- F

      DSA_lower <- sapply(2*(1:11)+1,
                          CKD_model,
                          input_timehorizon = timehorizon,
                          input_update_sub = update_sub,
                          input_egfrstage = egfrstage,
                          input_albstage = albstage,
                          input_age = age,
                          input_dia = dia / 100,
                          input_sex = sex / 100,
                          input_hyp = hyp / 100,
                          intervention = int, comparator = com)

      DSA_upper <- sapply(2*(1:11)+2,
                          CKD_model,
                          input_timehorizon = timehorizon,
                          input_update_sub = update_sub,
                          input_egfrstage = egfrstage,
                          input_albstage = albstage,
                          input_age = age,
                          input_dia = dia / 100,
                          input_sex = sex / 100,
                          input_hyp = hyp / 100,
                          intervention = int, comparator = com)

      reset_submodel <- CKD_model(1,
                                  input_timehorizon = timehorizon,
                                  input_update_sub = update_sub,
                                  input_egfrstage = egfrstage,
                                  input_albstage = albstage,
                                  input_age = age,
                                  input_dia = dia / 100,
                                  input_sex = sex / 100,
                                  input_hyp = hyp / 100,
                                  intervention = int, comparator = com)$CEAresults_s
      

      output$DSAPlot1 <- renderPlotly({
        tornado_plot(DSA_lower, DSA_upper, reset_submodel, 1) # For incremental cost
      })
      
      output$DSAPlot2 <- renderPlotly({
        tornado_plot(DSA_lower, DSA_upper, reset_submodel, 2) # For incremental QALY
      })
      
      output$DSAPlot3 <- renderPlotly({
        tornado_plot(DSA_lower, DSA_upper, reset_submodel, 3) # For ICER
      })
      
      # Remove the "Please wait..." message
      removeModal()

    })

    # For PSA
    observeEvent(input$run_PSA, {
      # Show the "Please wait..." message
      showModal(modalDialog(
        title = "Please wait...",
        "The Probabilistic Sensitivity Analysis is running. This can take several minutes.",
        easyClose = FALSE,
        footer = NULL
      ))
      
      timehorizon <- input$timehorizon
      egfrstage <- input$egfrstage
      albstage <- input$albstage
      age <- input$age
      dia <- input$dia
      sex <- input$sex
      hyp <- input$hyp
      int <- input$int
      com <- input$com
      update_sub <- F
      
      # # for debug
      # timehorizon <- "10 years"
      # egfrstage <- "G3a (45-59)" # unit: mL/min/1.73m^2
      # albstage <- "A1 (<30)" # unit: mg/g
      # age <- 57
      # dia <- 44.5 # % diabetes
      # sex <- 40.9 # % female
      # hyp <- 86.1 # % hypertensive
      # int <- "CV risk factor management for all eligible and full CKD management, no delay"
      # com <- "Usual care"
      # update_sub <- F
      
      PSA_CCQQ <- sapply(rep(2,1000),
                          CKD_model,
                          input_timehorizon = timehorizon,
                          input_update_sub = update_sub,
                          input_egfrstage = egfrstage,
                          input_albstage = albstage,
                          input_age = age,
                          input_dia = dia / 100,
                          input_sex = sex / 100,
                          input_hyp = hyp / 100,
                          intervention = int, comparator = com)
      
      reset_submodel <- CKD_model(1,
                                  input_timehorizon = timehorizon,
                                  input_update_sub = update_sub,
                                  input_egfrstage = egfrstage,
                                  input_albstage = albstage,
                                  input_age = age,
                                  input_dia = dia / 100,
                                  input_sex = sex / 100,
                                  input_hyp = hyp / 100,
                                  intervention = int, comparator = com)$CCQQresults_s
      
      # Create a data frame
      PSA_CCQQ_df <- data.frame(t(PSA_CCQQ))
      base_case_df <- data.frame(delta_cost = reset_submodel["Cost_int"] - reset_submodel["Cost_comp"],
                                 delta_qaly = reset_submodel["QALY_int"] - reset_submodel["QALY_comp"])
      
      # Calculate the differences in cost and QALY
      PSA_CCQQ_df$delta_cost <- PSA_CCQQ_df$Cost_int - PSA_CCQQ_df$Cost_comp
      PSA_CCQQ_df$delta_qaly <- PSA_CCQQ_df$QALY_int - PSA_CCQQ_df$QALY_comp
      
      # Define the range of willingness-to-pay thresholds
      max_wtp <- 400000
      min_wtp <- 150000
      wtp_thresholds <- seq(0, max_wtp, by = 1000)
      
      # Calculate the probability of each arm being cost-effective at each WTP threshold
      NMB <- sapply(wtp_thresholds, function(wtp) PSA_CCQQ_df$delta_qaly * wtp - PSA_CCQQ_df$delta_cost)
      prob_ce_int <- colMeans(NMB > 0)
      prob_ce_comp <- 1 - prob_ce_int
      
      # Calculate the expected value of perfect information (EVPI) at each WTP threshold
      EVPI <- sapply(wtp_thresholds, function(wtp) {
        NMB_each <- PSA_CCQQ_df$delta_qaly * wtp - PSA_CCQQ_df$delta_cost
        mean(pmax(NMB_each, 0)) - max(mean(NMB_each), 0)
      })
      
      # Find the range of WTP presented in the plot
      ## Apply LOWESS to smooth the CEAC
      lowess_fit <- lowess(wtp_thresholds, prob_ce_int)
      ## Find the minimum WTP where 
      ### 1) the probability of int is larger than 50%, and 
      ### 2) the probability is steady (diff < 0.0002)
      diff_lowess <- c(0, diff(lowess_fit$y))
      wtp_12 <- wtp_thresholds[lowess_fit$y > 0.5 & diff_lowess < 0.0002]
      wtp_found <- ifelse(length(wtp_12) == 0, max_wtp-50000, min(wtp_12))
      
      output$PSAPlot1 <- renderPlotly({
        # Set up the ranges
        ## Determine right-or-left-skewed
        space <- c(0.46, 0.32) # default: left-skewed
        if (mean(PSA_CCQQ_df$delta_cost) > median(PSA_CCQQ_df$delta_cost)) { # if right-skewed, turn around
          space <- space[c(2,1)]
        }
        ylower <- quantile(PSA_CCQQ_df$delta_cost, 0.005) 
        yupper <- quantile(PSA_CCQQ_df$delta_cost, 0.995) 
        yrange <- c(ylower - space[1] * (yupper - ylower), yupper + space[2] * (yupper - ylower))
        
        space <- c(0.46, 0.32) # default: left-skewed
        if (mean(PSA_CCQQ_df$delta_qaly) > median(PSA_CCQQ_df$delta_qaly)) { # if right-skewed, turn around
          space <- space[c(2,1)]
        }
        xlower <- quantile(PSA_CCQQ_df$delta_qaly, 0.005) 
        xupper <- quantile(PSA_CCQQ_df$delta_qaly, 0.995) 
        xrange <- c(xlower - space[1] * (xupper - xlower), xupper + space[1] * (xupper - xlower))
        
        # Create the scatter plot
        p <- plot_ly(PSA_CCQQ_df, x = ~delta_qaly, y = ~delta_cost, name = "PSA iterations", type = "scatter", mode = "markers", hoverinfo = "none", 
                     marker = list(color = "rgba(0, 176, 80, 0.5)", line = list(color = '#00B050', width = 1))) %>%
          add_trace(data = base_case_df, name = "Base case", marker = list(color = 'black', size = 10, line = NULL)) %>%
          layout(xaxis = list(title = "Incremental QALY", range = xrange, tickformat = ".1f"),
                 yaxis = list(title = "Incremental Cost (USD 2023)", range = yrange, tickformat = "$,.0f")) %>% 
          config(displayModeBar = F) 
      }) # For Cost-Effectiveness Scatter Plot
      
      output$PSAPlot2 <- renderPlotly({
        # Set up the ranges
        xmax <- max(wtp_found, min_wtp)
        xrange <- c(0, (ceiling(xmax/50000)+1)*50000+1000)
        
        # Create the CEAC
        p <- plot_ly(x = ~wtp_thresholds, y = ~prob_ce_int, name = "Intervention", type = 'scatter', mode = 'lines', yhoverformat = ".3f",
                line = list(color = '#00B050')) %>%
          add_trace(x = ~wtp_thresholds, y = ~prob_ce_comp, name = "Comparator", line = list(color = 'black')) %>%
          layout(xaxis = list(title = "Willingness-to-Pay per QALY (USD 2023)", range = xrange, tickformat = "$,.0f"),
                 yaxis = list(title = "Probability of Being Cost-Effective", range = c(0, 1), tickformat = ".1f")) %>% 
          config(displayModeBar = F) 
        
      }) # For CEAC
      
      output$PSAPlot3 <- renderPlotly({
        # Set up the ranges
        xmax <- max(wtp_found, min_wtp)
        xrange <- c(0, (ceiling(xmax/50000)+1)*50000+1000)
        
        ymax <- max(EVPI)
        yrange <- c(0, (ceiling(ymax/5000)+1)*5000+100)
        
        # Create the EVPI graph
        p <- plot_ly(x = ~wtp_thresholds, y = ~EVPI, type = 'scatter', 
                     mode = 'lines', name = "EVPI at WTP", hoverinfo = "name+y", 
                     line = list(color = '#00B050')) %>%
          layout(xaxis = list(title = "Willingness-to-Pay per QALY (USD 2023)", 
                              range = xrange, tickformat = "$,.0f"),
                 yaxis = list(title = "EVPI", 
                              range = yrange, tickformat = "$,.0f"))%>% 
          config(displayModeBar = F) 
      }) # For EVPI
      
      # Remove the "Please wait..." message
      removeModal()
      
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
