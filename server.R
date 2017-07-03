#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
get_group=function(sample_name){
  groupname=substr(as.character(sample_name),1,stop = as.numeric(gregexpr("\\#",sample_name))-1)
  return(groupname)
}
get_condition=function(sample_name){
  condition=substr(as.character(sample_name), as.numeric(gregexpr("\\#",sample_name))+1,nchar(as.character(sample_name)))
  return(condition)  
}

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  needReplicate <- reactive({
    switch(input$needReplicate,
           "TRUE" =TRUE,
           "FALSE"=FALSE)
  })
  dataGE <- reactive({
    switch(input$dataGE,
           "multipile results" ="multipile results",
           "summary"="summary")
  })

  
  output$replicateStaus<- renderPrint({
    replicateStaus<-needReplicate()
    print(replicateStaus)
  })
  output$dataGEstaus<- renderPrint({
    dataGEstaus<-dataGE()
    print(dataGEstaus)
  })
  output$analysisType<- renderPrint({
    print(input$analysis_type)
  })
  output$contents <- renderText({
    inFile <- input$sample_list
    if (is.null(inFile))
      return(NULL)
    sample_list_file=read.csv(inFile$datapath,header = TRUE,stringsAsFactors = FALSE)
    group_list_foreverysample=sapply(sample_list_file$sample_name,get_group)
    group_list=group_list_foreverysample[!duplicated(group_list_foreverysample)]
    condition_list_foreverysample=sapply(sample_list_file$sample_name,get_condition)
    condition_list=condition_list_foreverysample[!duplicated(condition_list_foreverysample)]
    print(group_list)
    #print(condition_list)
    
  })

  
})
