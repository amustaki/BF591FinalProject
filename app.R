library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(dplyr)
library(tidyverse)

options(shiny.maxRequestSize=30*1024^2)


ui <- fluidPage(
  titlePanel("Final Assignment BF591"), 
    tabsetPanel(
      tabPanel("Samples", 
               sidebarPanel(
                 fileInput("csv", paste0("Load Sample Information Matrix")),
                 radioButtons("button1", "Choose a column to sort the dataframe by", c("Sample_title", "Geo_accession", "Sample_type", "Sample_channel_count", "Source_name","Organism", "Tissue", "Diagnosis", "PMI","Age_of_death","RIN","Seq_reads")),
                 #radioButtons("button2", "Choose the column for the x-axis", c("Sample_title", "Geo_accession", "Sample_type", "Sample_channel_count", "Source_name","Organism", "Tissue", "Diagnosis", "PMI","Age_of_death","RIN","Seq_reads","PMI")), 
                 #radioButtons("button3", "Choose the column for the y-axis", c("Sample_title", "Geo_accession", "Sample_type", "Sample_channel_count", "Source_name","Organism", "Tissue", "Diagnosis", "PMI","Age_of_death","RIN","Seq_reads","Diagnosis")), 
                 radioButtons("button10", "Chose a continous variable for the histogram", c("Age_of_death","PMI","RIN", "Seq_reads")),
                 sliderInput(inputId = "bins", label = "Number of Bins for Histogram:", min = 1, max = 50, value = 30),
                 radioButtons("button11", "Chose a continous variable for the desnity plot", c("Age_of_death","PMI","RIN", "Seq_reads")),
                 colourInput("color1", "Color", "purple"),
                 colourInput("color2", "Color for ", "pink"),
                 
                 ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", tableOutput(outputId = "summary")),
                   tabPanel("Table", tableOutput(outputId = "table1")),
                   tabPanel("Plots", plotOutput(outputId = "histogram"),plotOutput(outputId = "density"))) 
                 ),
                ), 
               
      tabPanel("Counts",
               sidebarPanel(
                 fileInput("csv2", paste0("Load Normalized Counts")), 
                 sliderInput("threshold1", label = "Choose threshold on percentile variance", min = 0, max =100, value = 10),
                 sliderInput("threshold2", label = "Choose threshold for number of samples that are non-zero", min = 0, max = 100, value = 25), 
                 radioButtons("button4", "Choose the PC for the x-axis", c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                 radioButtons("button5", "Choose the PC for the y-axis", c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                ),
               
               mainPanel(
               tabsetPanel(
                 tabPanel("Summary", textOutput("text")), 
                 tabPanel("Diagonistic Scatter Plots", plotOutput("scatter1"), plotOutput("scatter2")), 
                 tabPanel("Clustered Heatmap", plotOutput("heatmap1")), 
                 tabPanel("Principal Component Analysis", plotOutput("PCscatter")))
               ),
              ),
              
      tabPanel("DE",
               sidebarPanel(
                 fileInput("csv3", paste0("Load in Differential Expression Results")), 
                 helpText("Here are the options for formatting the restults tab"),
                 radioButtons("button8", "Chosse a column to sort the results table by", c("baseMean", "log2FoldChange", "lfcSE","stat","pvalue","padj","lfcSE")),
                 helpText("Here are the options for formatting the plot and table tab"),
                 helpText("A volcano plot can be generated with log2 fold-change on the x-axis and p-adjusted on the y-axis"),
                 radioButtons("button", "Choose the column for the x-axis", c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                 radioButtons("but", "Choose the column for the y-axis", c("baseMean","log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                 colourInput("col","Base Point Color", "purple"),
                 colourInput("color", "Highlight Point Color", "blue"), 
                 sliderInput("range", label = "Select the magnitude of the p adjusted coloring:", min = -39, max = 0,value = -15),
                 #actionButton("b","Plot"),
                 #actionButton("c", "Table")
                 ),
                
               mainPanel(
               tabsetPanel(
                 tabPanel("Results", tableOutput("tabledata")),
                 tabPanel("Table", tableOutput("table2")), 
                 tabPanel("Plot", plotOutput("volcano")))
               ), 
              ),
      tabPanel("GSEA",
               sidebarPanel(
                 fileInput("csv4", paste0("Load in FGSEA Results from the Differential Expression Data")),
                 helpText("A bar plot of FGSEA NES for top pathways based on the slider below"),
                 sliderInput("threshold3", label = "Select the adjusted p-value to determine how many pathway are on the plot", min = -.3, max = .3, value = 0.01),
                 helpText("For the sortable table tab here are the options:"),
                 sliderInput("threshold4", label = "Select the adjusted p-value to filter the table", min = -.3, max = .3, value = 0.01),
                 radioButtons("button6", "Choose which pathways you want to see", c("All", "positive", "negative")),
                 radioButtons("button7", "Choose the column you want to sort the table by", c("pathway","pval", "padj", "log2err", "NES","size")),
                 downloadButton("down", label = "Download Filtered Table"),
                 helpText("For the scatterplot tab here are the options:"),
                 sliderInput("threshold5", label = "Select the adjusted p-value to filter the data", min = -.3, max = .3, value = 0.01)),
                 
                mainPanel(
                tabsetPanel(
                  tabPanel("Barplot", plotOutput("barplot")),
                  tabPanel("Pathway Table", tableOutput("pathtable")),
                  tabPanel("Scatterplot", plotOutput("NEScatter")))
                 ),
               ),
),
)

#load in all the four csvs that are necessary for the different tabs 
server <- function(input, output, session) {
  load_data1 <- reactive({
    if(is.null(input$csv3))
      return(NULL)
    df <- read.csv(input$csv3$datapath, header = TRUE, sep = ",")
    return(df)
  })
  
  load_data2 <- reactive({
    if(is.null(input$csv2))
      return(NULL)
    df <- read.csv(input$csv2$datapath, header = TRUE, sep = ",")
    return(df)
  })
  
  load_data3 <- reactive({
    if(is.null(input$csv4))
      return(NULL)
    df <- read.csv(input$csv4$datapath, header = TRUE, sep = ",")
    return(df)
  })
  
  load_data4 <- reactive({
    if(is.null(input$csv))
      return(NULL)
    df <- read.csv(input$csv$datapath, header = TRUE, sep = ",")
    return(df)
  })
  
  #tab 1 Sample Matrix 
  summary_exp <- function(data){
    if(is.null(input$csv))
      return(NULL)
    samples <- select(data, -Sample_title,-Geo_accession)
    names <- colnames(samples)
    samp <- drop_na(data)
    #tried using apply and ifelse but it wasnt working properly so then I did it semi manually to get the summary
    meanRIN <- mean(samp$RIN)
    meanagedeath <-round(mean(samp$Age_of_death),3)
    meanPMI <- round(mean(samp$PMI),3)
    meanSeq_reads <- round(mean(samp$Seq_reads),3)
    meansamplecount <- round(mean(samp$Sample_channel_count),3)
    meanDiag <- c("Neurologically Normal,Huntington's Disease")
    meanorg <- c("Homo sapiens")
    meansource <- c("Brain")
    meantis <- c("Prefrontal Cortex")
    meantype <- c("SRA")
    samplesum <- data.frame(names)
    totalmeans <- c(meantype,meansamplecount,meansource,meanorg,meantis,meanDiag,meanPMI,meanagedeath,meanRIN,meanSeq_reads)
    #determined this by using class for each of the column names 
    types <- c("Character", "Integer", "Character", "Character","Character","Character","Numeric","Integer","Numeric","Integer")
    samplesum <-cbind(samplesum,totalmeans)
    samplesummary <- cbind(samplesum,types)
    names(samplesummary) <- c("Column Name", "Mean or Distinct Type", "Type" )
    return(samplesummary)
    
  }
  
  sortdata <- function(data,col_name){
    if(is.null(input$csv))
      return(NULL)
    #col_name<-as.vector(col_name)
    col_name <- data[,c(col_name)] #specify a column by name 
    sorted <- data[order(col_name),]
    return(sorted)
  }
  
  histogram<- function(data, x_name, color1,slider){
    if(is.null(input$csv))
      return(NULL)
    col_name <- data[,c(x_name)]
    graph <- ggplot(data, aes(x = col_name)) +
      geom_histogram(bins = slider, col = color1) +
      xlab(x_name) +
      ylab("Counts") +
      ggtitle(paste0("Histogram of ", x_name , " vs Count"))
    #scale_color_manual(values = c(color1))
    return(graph)
  }
  
  density <- function(data, x_name,color1,color2){
    if(is.null(input$csv))
      return(NULL)
    col_name <- data[,c(x_name)]
    graph <- ggplot(data, aes(x = col_name)) +
      geom_density(fill=color1) +
      geom_vline(aes(xintercept=mean(col_name)), color = color2, linetype = "dashed", size = 1) +
      xlab(x_name) +
      ylab("Density") +
      ggtitle(paste0("Density Plot of ", x_name))
    return(graph)
  }
  
  #tab 2 Counts Matrix Functions 
  summarized <- function(data,threshold1, threshold2) {
    if(is.null(input$csv2))
      return(NULL)
    row.names(data) <- data$X
    data <- select(data, -X)
    data <- log(data) #log the data
    data_var <-  mutate(data, variance = apply(data,1,var)) #get the variance of the data
    data_ordered<-data_var[order(data_var$variance,decreasing = T),] #order so that you can get the top results for the percentile filtering
    thres <- threshold1/(100*length(data)) #get the threshold based on the length of the dataset 
    passed <- data_ordered[data_ordered$variance > thres,] #figure out which ones pass the threshold 
    passed <- select(passed, -variance) #drop the variance column so it doesnt mess things up with extra columns 
    data_new <- select(data_ordered, -variance)
    #filter 2 is include genes with at least X samples that are non-zero
    #so count when the value is greater than 0 and when its over the threshold keep it if not get rid of the row 
    thres2<- threshold2/(100*length(data))
    data_nonzero <- mutate(data_new, zeros = rowSums(data_var == 0))
    passed2 <- data_new[data_new$zeros > thres2,]
    total <- rbind(passed, passed2) %>% distinct() #combine the two filters to get the total and remove the duplicates
    num_samples <- ncol(data)
    num_genes <-nrow(data)
    #number and percent of genes passing current filter 
    num_passing <- nrow(total)
    percent_passing <- round((num_passing/num_genes)*100,3)
    num_failing <-num_genes - num_passing 
    percent_failing <- round((100 - percent_passing),3)
    text <- (paste0("The number of samples is ", num_samples , " and the number of genes is ", num_genes, ".", "The number of genes passing the current filter is ", num_passing, " and the percent passing is ", percent_passing,"%.", "The number of genes failing the current filter is ", num_failing , " and the percent failing is ", percent_failing, "%."))
    return(text)
  }
  
  scatter1 <- function(data, threshold1){
    if(is.null(input$csv2))
      return(NULL)
    row.names(data) <- data$X
    data <- select(data, -X)
    data <- log10(data) #follow the same process as above this is the scatter plot for the first filter with varinace
    data_var <-  mutate(data, variance = apply(data,1,var), median = apply(data,1, median))
    data_ordered<-data_var[order(data_var$variance,decreasing = T),]
    thres <- threshold1/(100*length(data))
    passed <- data_ordered[data_ordered$variance > thres,] 
    labelleddata <- mutate(data_ordered, passf1 = ifelse((row.names(data) %in% row.names(passed)), "TRUE", "FALSE")) #assign true and false to get the labels for the colors on what passed on not 
    graph <- ggplot(labelleddata, aes(x=median, y = variance, color = passf1)) +
      geom_point()+
      xlab("Log Median Count") +
      ylab("Log Variance") +
      ggtitle("Scatter Plot of Median Count vs Variance ") 
    #had alot of issues with this plot not sure if it is correct but I tried my hardest, got different result when i logged the data in the plot instead of the beggining and log scaled y 
    return(graph)
  }

  scatter2 <- function(data, threshold2, color){
    if(is.null(input$csv2))
      return(NULL)
    row.names(data) <- data$X
    data <- select(data, -X) #this is the second scatter plot for the second filter non zeros 
    labelleddata <- mutate(data, median = apply(data,1,median), zeros = rowSums(data == 0), passf2 = (ifelse(zeros < threshold2/(100*length(data)), "TRUE", "FALSE")))
    graph <- ggplot(labelleddata, aes(x=log(median), y = zeros , color = passf2)) +
      geom_point()+
      xlab("Median Count") +
      ylab("Number of Zeros") +
      ggtitle("Scatter Plot of Median Count vs Number of Zeros ")
    return(graph)
  }
  
  #add in the heatmap function 
  heatmaplot <- function(data,threshold1, threshold2){
    if(is.null(input$csv2))
      return(NULL)
    row.names(data) <- data$X #same filtering as before 
    data <- select(data, -X)
    data <- log(data)
    data_var <-  mutate(data, variance = apply(data,1,var))
    data_ordered<-data_var[order(data_var$variance,decreasing = T),]
    thres <- threshold1/(100*length(data))
    passed <- data_ordered[data_ordered$variance > thres,]
    passed <- select(passed, -variance)
    data_new <- select(data_ordered, -variance)
    thres2<- threshold2/(100*length(data))
    data_nonzero <- mutate(data_new, zeros = rowSums(data_var == 0))
    passed2 <- data_nonzero[data_nonzero$zeros > thres2,]
    passed2 <- select(passed2, -zeros)
    total <- rbind(passed, passed2) %>% distinct()
    data1<- as.matrix(total)
    heat<- heatmap(data1, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))  
    legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))
    return(heat)
  }
  
  #pca stuff
  plot_pca <- function(data, x_name, y_name) {
    if(is.null(input$csv2))
      return(NULL)
    print(x_name)
    row.names(data) <- data$X
    data <- select(data, -X)
    transposed <- t(data) #transpose the data 
    pca <- prcomp(transposed)
    sd<- pca$sdev
    variance <- (sd)^2/sum((sd)^2) 
    pcs <- colnames(pca$rotation)
    cum<- cumsum(variance)
    pca_tibble <- tibble(variance_explained = variance, principal_components = as.factor(pcs), cumulative = cum)
    pca <- as.data.frame(pca$x)
    colx <- pca[,c(x_name)]
    coly <- pca[,c(y_name)]
    #subset out to only have the number of the PC given so that it can be used to put the variance on the plot
    #turn from string into an integer using strtoi 
    pc_num_x <-strtoi(substr(x_name,3,3))
    pc_num_y <- strtoi(substr(y_name,3,3))
    pcaplot <- ggplot(pca, aes(x=colx, y = coly)) +
      geom_point() +
      ggtitle("Scatter Plot of PCA Plot Projections") +
      xlab(paste0(x_name, " ", round(variance[pc_num_x]*100), "% variance")) +
      ylab(paste0(y_name, " ", round(variance[pc_num_y]*100), "% variance"))
    #need to add in the varinace for each 
    return(pcaplot)
  }
  
  #tab3 DE functions 
  draw_table1 <- function(dataf,col_name) {
    if(is.null(input$csv3))
      return(NULL)
      dataf %>% mutate(pvalue = formatC(.$pvalue, digits = 4, format = 'e'),
               padj = formatC(.$padj, digits = 4, format = 'e')) 
    col_name <- dataf[,c(col_name)] #specify a column by name 
    sorted <- dataf[order(col_name),]
      return(sorted[1:100,])
  }
  #functions from homework 7
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    if(is.null(input$csv3))
      return(NULL)
    plot<- ggplot(dataf,aes(x = !!sym(x_name), y = -log10(!!sym(y_name)))) +
      geom_point(aes(colour = cut(-log(padj), c(-Inf, -(slider), Inf)))) +
      xlab(x_name) +
      ylab(y_name) +
      ggtitle("Volcano Plot for DESeq Differential Expression Results") +
      scale_color_manual(values = c(color1, color2))
    return(plot)
  }
  
  draw_table <- function(dataf, slider) {
    if(is.null(input$csv3))
      return(NULL)
    filt_dataf <- filter(dataf, padj< 10^slider) %>%  #dataf[!(dataf$padj > 10^slider),]
      mutate(pvalue = formatC(.$pvalue, digits = 4, format = 'e'),
             padj = formatC(.$padj, digits = 4, format = 'e'))
    return(filt_dataf)
  }
  
  #tab4 fgsea functions
  barplot <- function(data, threshold){
    if(is.null(input$csv4))
      return(NULL)
    filtered <- data[data$padj > threshold,]
    bargraph <- ggplot(filtered, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj < 0.05)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="Hallmark pathways NES from GSEA") + 
      theme_minimal()
    return(bargraph)
  }
  
  
 #table for pathways gsea
  second <- function(data, threshold, but, col_name){
    if(is.null(input$csv4))
      return(NULL)
    filtered <- data[data$padj > threshold,]
    if(but == "positive"){
      fil <- filter(filtered, NES > 0)
      col_name <- fil[,c(col_name)] #specify a column by name 
      sorted <- fil[order(col_name),]
      return(sorted)
    }
    else if(but == "negative"){
      fil <- filter(filtered, NES < 0)
      col_name <- fil[,c(col_name)] #specify a column by name 
      sorted <- fil[order(col_name),]
      return(sorted)
    }
    else{
      col_name <- filtered[,c(col_name)] #specify a column by name 
      sorted <- filtered[order(col_name),]
      return(sorted)
    }
  }

  
  #fgsea scatter plot
  third <- function(data, threshold){
    if(is.null(input$csv4))
      return(NULL)
    passed <- data[data$padj > threshold, "pathway"]
    newdat <- mutate(data, pass = ifelse(data$pathway %in% passed, "TRUE", "FALSE"))
    scat <- ggplot(newdat,aes(NES, -log10(padj), color = pass)) +
      geom_point()+
      scale_color_manual(values = c("#00AFBB", "grey")) +
      xlab("NES") +
      ylab("-log10(padj)") +
      ggtitle("Scatter Plot of NES vs -log10(padj) ")
    return(scat)
  }

  
  #tab 1 sample matrix tab 
  output$summary <- renderTable(summary_exp(load_data4()))
  output$table1 <- renderTable(sortdata(load_data4(), input$button1))
  output$histogram <- renderPlot(histogram(load_data4(), input$button10, input$color1, input$bins))
  output$density <- renderPlot(density(load_data4(), input$button11, input$color1, input$color2))

  #tab 2 counts tab 
  output$text <- renderText(summarized(load_data2(), input$threshold1, input$threshold2))
  output$scatter1 <- renderPlot(scatter1(load_data2(), input$threshold1))
  output$scatter2 <- renderPlot(scatter2(load_data2(), input$threshold2))
  output$heatmap1 <- renderPlot(heatmaplot(load_data2(), input$threshold1, input$threshold2))
  output$PCscatter <- renderPlot(plot_pca(load_data2(), input$button4, input$button5))
  
  #tab 3 DE tab 
  output$volcano <-renderPlot(volcano_plot(load_data1(),input$button, input$but, input$range, input$col, input$color))
  output$table2 <- renderTable(draw_table(load_data1(), input$range))
  output$tabledata <- renderTable(draw_table1(load_data1(), input$button8))
 
  #tab 4 GSEA tab 
  output$barplot <- renderPlot(barplot(load_data3(),input$threshold3))
  output$pathtable <- renderTable(second(load_data3(),input$threshold4, input$button6, input$button7))
  output$NEScatter <- renderPlot(third(load_data3(),input$threshold5))
  output$down <- downloadHandler(filename = function() {
  
    paste("data-", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    write.csv(second(load_data3(),input$threshold4, input$button6,input$button7),file)
  }
  )
  
}


shinyApp(ui = ui, server = server)
