library(shiny)
options(shiny.host = '0.0.0.0')
options(shiny.port = 6204)

infectedDataPath = "time_series_covid19_confirmed_global.csv"
recoveredDataPath = "time_series_covid19_recovered_global.csv"

ui <- fluidPage(
    
    # Application title
    titlePanel("COVID-19 Spread Prediction"),
    textOutput("appInfo1"),
    textOutput("appInfo2"),
    textOutput("appInfo3"),
    tags$hr(),
    
    fluidRow(
        
        column(4,
               wellPanel(
                   selectInput(inputId = "Country", label = strong("Country"),
                               choices = c("All"),
                               selected = "All"),
                   selectInput(inputId = "Province", label = strong("State/Province"),
                               choices = c("All"),
                               selected = "All"
                   ),
                   
                   # Eventually deal with Exponential regression
                   # radioButtons("modelType", "Model type:",
                   #              c("Linear" = "linear",
                   #                "Exponential" = "exp")),
                   
                   verbatimTextOutput("noRecovered"),
                   
                   h5(textOutput("maxDate")),
                   verbatimTextOutput("info2"),
                   
                   h5(textOutput("predictionData")),
                   tableOutput('predictionDataTable'),
                   
                   tags$hr(),
                   h5(textOutput("dataOrigin")),
                   tags$a(href="https://github.com/CSSEGISandData/COVID-19", "COVID-19 Data")
               )      
        ),
        
        column(8,
               plotOutput("distPlot2", brush = "plot_brush"),
               plotOutput("distPlot"),
               plotOutput("distPlot3")
        )
    )
)

server <- function(input, output, session) {
    observe({
        updateSelectInput(session, "Country",
                          choices = c("All", unique(coronaRaw()$Country.Region))
        )
    })
    
    observe({
        updateSelectInput(session, "Province",
                          choices = c("All", unique(coronaRaw()$Province.State[which(coronaRaw()$Country.Region == input$Country)]))
        )
    })
    
    coronaRaw <- reactive({
        invalidateLater(1000 * 60 * 60 * 3)
        
        read.csv(infectedDataPath, stringsAsFactors=FALSE)
    })
    
    coronaRec <- reactive({
        invalidateLater(1000 * 60 * 60 * 3)
        
        read.csv(recoveredDataPath, stringsAsFactors=FALSE)
    })
    
    infectedData <- function(){
        infected = 0
        if(input$Country == "All"){
            infected = apply(coronaRaw()[, -c(1:4)], 2, sum)
            recovered = apply(coronaRec()[, -c(1:4)], 2, sum)
        }
        else if(input$Province == "All"){
            infected = apply(coronaRaw()[coronaRaw()$Country.Region == input$Country, -c(1:4)],2, sum)
            recovered = apply(coronaRec()[coronaRec()$Country.Region == input$Country, -c(1:4)],2, sum)
        }
        else {
            infected = apply(coronaRaw()[coronaRaw()$Province.State == input$Province, -c(1:4)],2, function(a){return(a)})
            recovered = apply(coronaRec()[coronaRec()$Province.State == input$Province, -c(1:4)],2, function(a){return(a)})
        }
        
        isRecoveredDataAvailable = ""
        if(length(recovered) == 0){
            isRecoveredDataAvailable = "no recovered data available"
            data = as.vector(infected)
        }
        else{
            data = as.vector(infected - recovered)
        }
        output$noRecovered = function(){return(isRecoveredDataAvailable)}
        return(data)
    }
    
    infectionRate <- function(infected){
        nextInfect = c(infected, 1)
        prevInfect = c(1, infected)
        rate = (nextInfect / prevInfect)[-c(1,length(nextInfect))]
        
        return(rate)
    }
    
    infectedPlotData <- function(infected){
        dates = as.Date(paste0(substring(names(coronaRaw()), 2)[-c(1:4)], "20"), "%m.%d.%Y")
        data = data.frame(Date = dates, Infected = infected)
        
        return(data)
    }
    
    ratePlotData <- function(infected){
        rates = infectionRate(infected)
        dates = as.Date(paste0(substring(names(coronaRaw()), 2)[-c(1:4)], "20"), "%m.%d.%Y")
        data = data.frame(Date = dates[-1], Rate = rates)
        
        return(data)
    }
    
    getSelectedData <- function(){
        infected = infectedData()
        data = ratePlotData(infected)
        
        selectedData = brushedPoints(data, input$plot_brush, xvar = "Date", yvar = "Rate")
        
        return(selectedData)
    }
    
    getDownwardsModel <- function(){
        selectedData = getSelectedData()
        
        if(dim(selectedData)[1] == 0){
            return(NULL)
        }
        
        # Maybe adding exponential back into the mix
        #if(input$modelType == "linear"){
        return(getDownwardsLinearModel(selectedData))
        #}
        # else if(input$modelType == "exp"){
        #     return(getDownwardsodelExponentialModel(selectedData))
        # }
    }
    
    getDownwardsLinearModel <- function(selectedData){
        return(list(model = lm(selectedData$Rate ~ selectedData$Date), selectedData = selectedData))
    }
    
    # Currently not used
    getDownwardsodelExponentialModel <- function(selectedData){
        selectedData = cbind(index = 1:length(selectedData$Date), selectedData)
        
        theta.0 <- min(selectedData$Rate) * 0.5  
        
        # Estimate the rest parameters using a linear model
        model.0 <- lm(log(Infected - theta.0) ~ index, data=selectedData)  
        alpha.0 <- exp(coef(model.0)[1])
        beta.0 <- coef(model.0)[2]
        
        start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
        
        model <- nls(Infected ~ alpha * exp(beta * index) + theta , data = selectedData, start = start)
        
        
        return(list(model = model, selectedData = selectedData))
    }
    
    getDownwardsModelPlotData <- function(modelData){
        modelPlot = predict(modelData[[1]],list(Date=modelData[[2]]$Date))
        
        return(data.frame(Date = modelData[[2]]$Date, Rate = modelPlot))
    }
    
    # Currently not used
    predictFutureExp <- function(infected, rates, modelData){
        currentRate = rates[length(rates)]
        currentInfected = tail(infected, n=1)
        
        prediction = infected;
        maxIterations = 100
        iterationCounter = 0;
        while (tail(prediction, n=1)$Infected > 1 && iterationCounter < maxIterations) {
            iterationCounter = iterationCounter + 1
            currentInfectedData = tail(prediction, n=1)
            nextDate = currentInfectedData$Date + 1
            nextRate = predict(modelData[[1]],data.frame(Date=nextDate))
            nextInfected = currentInfectedData$Infected * nextRate
            prediction = rbind(prediction, data.frame(Date = nextDate, Infected = nextInfected))
        }
        
        return(prediction)
    }
    
    predictFuture <- function(infected, rates, rateOfChange){
        if(is.null(rateOfChange) || rateOfChange > 0)
            return(NULL)
        currentRate = rates[length(rates)]
        currentInfected = tail(infected, n=1)
        
        prediction = infected;
        
        while (tail(prediction, n=1)$Infected > 1) {
            currentRate = currentRate + rateOfChange
            lastRow = tail(prediction, n=1)
            lastRow$Date = lastRow$Date + 1
            lastRow$Infected = lastRow$Infected * currentRate
            prediction = rbind(prediction, lastRow)
        }
        
        maxDate = max(prediction$Date[which.max(prediction$Infected)])
        
        output$info2 = function(){return(maxDate)}
        
        return(list(prediction, maxDate))
    }
    
    output$distPlot <- renderPlot({
        infected = infectedData()
        data = infectedPlotData(infected)
        plot(data, main="Infected")
    })
    
    output$distPlot2 <- renderPlot({
        infected = infectedData()
        data = ratePlotData(infected)
        plot(data$Date, data$Rate, main="Rate of Infection - Select data to use in the prediction", xlab = "Date", ylab = "Rate")
        model = getDownwardsModel()
        if(!is.null(model)){
            downwardsModelPlotData = getDownwardsModelPlotData(model)
            lines(downwardsModelPlotData$Date, downwardsModelPlotData$Rate, col = "red", lwd=2)
        }
    })
    
    output$distPlot3 <- renderPlot({
        infected = infectedData()
        infectedData = infectedPlotData(infected)
        rates = infectionRate(infected)
        model = getDownwardsModel()
        if(!is.null(model)){
          predictionData = predictFuture(infectedData, rates, model[[1]][1]$coefficients[2])
          if(!is.null(predictionData)){
            plot(predictionData[[1]], main="Infection Prediction")
          }
        }
    })
    
    formatDates <- function(data){
      data$Date = as.character(data$Date)
      return(data)
    }
    
    output$predictionDataTable <- renderTable(formatDates(getSelectedData()))
    
    output$maxDate <- function(){return("Max Date")}
    
    output$predictionData <- function(){return("Prediction Data")}
    
    output$dataOrigin <- function(){
        return("Data Source")
    }
    
    output$info <- function(){
        return("The data used on this page comes from: ")
    }
    
    output$appInfo1 <- function(){
        return("This app's purpose is to try to determine COVID-19's progression. Ultimatly, it tries to estimate when will it hit it's perak.")
    }
    output$appInfo2 <- function(){
        return("In the second plot, select the data that you want to use in order to estimate the infection rates progression.")
    }
    output$appInfo3 <- function(){
        return("A linear regression is applied to the selected rates in order to estimate the progression. NB: Other regressions might be better.")
    }
}

# Run the application 
shinyApp(ui = ui, server = server)
