#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

shinyUI(
  fluidPage(
  
  # Application title
  titlePanel("TOFU"),
  sidebarLayout(
    sidebarPanel(
       selectInput("needReplicate","whether need replicate in result",choices = c("TRUE","FALSE")),
       selectInput("dataGE","GE result is one summary or multiple results",choices = c("summary","multipile results")),
       checkboxGroupInput("analysis_type","GE,ASE or both",choices = c("Gene Expression","ASE")),
       fileInput('sample_list','upload sample list ',accept = c('text/csv', 
                                                                'text/comma-separated-values,text/plain', 
                                                                '.csv')),
       #tags$hr(),      
       submitButton("Run")
       
      )
    ,
  mainPanel(
    verbatimTextOutput("replicateStaus"),
    verbatimTextOutput("dataGEstaus"),
    verbatimTextOutput("analysisType"),
    tableOutput('contents')
    )
  
)))
