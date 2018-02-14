library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(plyr)
library(reshape2)
library(limma)
library(biomaRt)
library(clusterProfiler)
library(GeneAnswers)
library(pathview)
shinyServer(function(input, output) {
  
# -------------- The module for Shiny first part ------------------------ #
  
  # Function for filtering the data
  SubsetShinyBasic <- function(data, pubid = "", gpl = "", gse = "", year = "--") {
    if (pubid != "") {
      datapub <- data[grep(pubid, data$PubmedId), ]
    } else {
      datapub <- data
    }
    if (gpl != "") {
      datagpl <- data[grep(gpl, data$GPL), ]
    } else {
      datagpl <- data
    }
    if (gse != "") {
      datagse <- data[grep(gse, data$GSE), ]
    } else {
      datagse <- data
    }
    if (year != "--") {
      datatime <- data[grep(year, data$Year), ]
    } else {
      datatime <- data
    }
    return(join_all(list(datapub, datagpl, datagse, datatime), by = colnames(bcbasic), type = "inner"))
  }
  
  # Function for show the message of every research
  BasicDataMessage <- reactive({
    SubsetShinyBasic(data = bc, pubid = input$pubmedidbasic, gpl = input$gplnumberbasic,
                     gse = input$gsenumberbasic, year = input$yearbasic)
  })
  
  # Process data for table showing
  BasicDataTable <- reactive({
    SubsetShinyBasic(data = bcbasic, pubid = input$pubmedidbasic, gpl = input$gplnumberbasic,
                     gse = input$gsenumberbasic, year = input$yearbasic)
  })
  output$df <- renderDataTable(BasicDataTable(), rownames = F,
                               options = list(pageLength = 5, scrollY = "190px", scrollX = "500px"))
  
  # Process data for message showing
  BasicDataMessage <- reactive({
    SubsetShinyBasic(data = bc, pubid = input$pubmedidbasic, gpl = input$gplnumberbasic,
                     gse = input$gsenumberbasic, year = input$yearbasic)
  })
  output$GSESummary <- renderText(ifelse(length(unique(BasicDataMessage()$GSE_Summary)) != 1, "Please don't choose more than one research", unique(BasicDataMessage()$GSE_Summary)))
  output$GPLTitle <- renderText(ifelse(length(unique(BasicDataMessage()$GSE_Summary)) != 1, "Please don't choose more than one research", unique(BasicDataMessage()$GPL_Title)))
  
  
# ----------------------------- The module for Shiny Second Part ---------------------------------- #
  # Clinical Data Process
  SubsetShinyClinic <- function(data, gpl = "", gse = "", gsm = "", colchoose) {
    if (gpl != "") {
      datagpl <- data[grep(gpl, data$GPL), ]
    } else {
      datagpl <- data
    }
    if (gse != "") {
      datagse <- data[grep(gse, data$GSE), ]
    } else {
      datagse <- data
    }
    if (gsm != "") {
      datagsm <- data[grep(gsm, data$GSM), ]
    } else {
      datagsm <- data
    }
    df <- join_all(list(datagpl, datagse, datagsm), by = colnames(bcclinic), type = "inner")
    return(df[, colchoose])
  }
  
  # Process data for table showing
  ClinicalData <- reactive({
    SubsetShinyClinic(data = bcclinic, gpl = input$gplnumberclinic, gse = input$gsenumberclinic, gsm = input$gsmnumberclinic, colchoose = input$cliniccol)
  })
  output$Cdf <- renderDataTable(ClinicalData(), rownames = F,
                                options = list(pageLength = 5, scrollX = "500px", scrollY = "190px"))
  
  # ---------------------------------- The module for Shiny Third Part ---------------------------- #
  
  # Expression Data Process
  ExpressionMatrix <- function(gpl, gse, n) {
    MatrixNum <- paste(gpl, gse, "express", "genetype", sep = "_")
    data <- eval(parse(text = MatrixNum))
    data$AVGExp <- apply(data[, 2:(ncol(data) - 1)], 1, mean)
    data$AVGExp <- signif(data$AVGExp, 5)
    dataTopX <- ddply(data, "GeneType", function(x) x[order(-x$AVGExp)[1:n],])
    return(list("Coding" = dataTopX[which(dataTopX$GeneType == "Coding"), ], 
                "Noncoding" = dataTopX[which(dataTopX$GeneType == "Noncoding"), ],
                "Matrix" = data,
                "Summary" = dataTopX))
  }
  
  # Process data for table
  MatrixData <- reactive({
    ExpressionMatrix(gpl = input$gplnumberexp, gse = input$gsenumberexp, n = input$TopXtoShow)
  })
  output$Edf <- renderDataTable(MatrixData()$Summary[,c("GeneSymbol", "GeneType", "AVGExp")], rownames = F,
                                options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  output$EdfMatrix <- renderDataTable(MatrixData()$Matrix, rownames = F,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  #Plot heatmap
  output$HeatmapCoding <- renderPlot({
    dataCoding <- melt(data = MatrixData()$Coding, id.vars = c("GeneSymbol", "GeneType", "AVGExp"))
    dataNoncoding <- melt(data = MatrixData()$Noncoding, id.vars = c("GeneSymbol", "GeneType", "AVGExp"))
    scaledata <- scale(c(dataCoding$value, dataNoncoding$value))
    dataCoding$value <- scaledata[1:nrow(dataCoding)]
    dataNoncoding$value <- scaledata[(nrow(dataCoding) + 1):length(scaledata)]
    maxExp <- max(dataCoding$value, dataNoncoding$value)
    minExp <- min(dataCoding$value, dataNoncoding$value)
    pCoding <- ggplot(dataCoding, aes(x = variable, y = GeneSymbol, fill = value)) + labs(x = NULL, y = NULL) + geom_tile(color = "white", size = 0.1)
    pCoding <- pCoding + scale_fill_gradient(limits = c(minExp - 1, maxExp + 1), low = "green", high = "red", name = "Expression")
    pCoding <- pCoding + theme(axis.ticks = element_blank(), legend.key.height = unit(1, "cm"), legend.key.width = unit(0.5, "cm"), panel.background = element_blank(), panel.grid = element_blank(), axis.text.x = element_blank(), legend.title = element_text(size = 9), legend.position = "right")
    pCoding
  })
  output$HeatmapNoncoding <- renderPlot({
    dataCoding <- melt(data = MatrixData()$Coding, id.vars = c("GeneSymbol", "GeneType", "AVGExp"))
    dataNoncoding <- melt(data = MatrixData()$Noncoding, id.vars = c("GeneSymbol", "GeneType", "AVGExp"))
    scaledata <- scale(c(dataCoding$value, dataNoncoding$value))
    dataCoding$value <- scaledata[1:nrow(dataCoding)]
    dataNoncoding$value <- scaledata[(nrow(dataCoding) + 1):length(scaledata)]
    maxExp <- max(dataCoding$value, dataNoncoding$value)
    minExp <- min(dataCoding$value, dataNoncoding$value)
    pNoncoding <- ggplot(dataNoncoding, aes(x = variable, y = GeneSymbol, fill = value)) + labs(x = NULL, y = NULL) + geom_tile(color = "white", size = 0.1)
    pNoncoding <- pNoncoding + scale_fill_gradient(limits = c(minExp - 1, maxExp + 1), low = "green", high = "red", name = "Expression")
    pNoncoding <- pNoncoding + theme(axis.ticks = element_blank(), legend.key.height = unit(1, "cm"), legend.key.width = unit(0.5, "cm"), panel.background = element_blank(), panel.grid = element_blank(), axis.text.x = element_blank(), legend.title = element_text(size = 9), legend.position = "right")
    pNoncoding
  })
  
  # An interactive UI for multiple gene correlation
  output$COselectUIM <- renderUI({
    data <- MatrixData()$Coding
    genes <- as.character(data$GeneSymbol)
    selectInput("COselectM", label = "Select Coding genes", choices = as.list(setNames(c("", genes), c("", genes))), selected = "", multiple = T)
  })
  
  output$NONselectUIM <- renderUI({
    data <- MatrixData()$Noncoding
    genes <- as.character(data$GeneSymbol)
    selectInput("NONselectM", label = "Select Noncoding genes", choices = as.list(setNames(c("", genes), c("", genes))), selected = "", multiple = T)
  })
  
  # For single gene
  output$COselectUIS <- renderUI({
    data <- MatrixData()$Coding
    genes <- as.character(data$GeneSymbol)
    selectInput("COselectS", label = "Select Coding gene", choices = as.list(setNames(c("", genes), c("", genes))), selected = "")
  })
  
  output$NONselectUIS <- renderUI({
    data <- MatrixData()$Noncoding
    genes <- as.character(data$GeneSymbol)
    selectInput("NONselectS", label = "Select Noncoding gene", choices = as.list(setNames(c("", genes), c("", genes))), selected = "")
  })
  
  CorreDataM <- function(gpl, gse, noncoding, coding) {
    if (noncoding == "" | coding == "") {
      return(data.frame())
    } else {
      MatrixNum <- paste(gpl, gse, "express", "genetype", sep = "_")
      data <- eval(parse(text = MatrixNum))
      non <- data[data$GeneSymbol %in% noncoding, 2:(ncol(data) - 1)]
      co <- data[data$GeneSymbol %in% coding, 2:(ncol(data) - 1)]
      matrix <- rbind(non, co)
      row.names(matrix) <- c(noncoding, coding)
      matrix <- t(matrix)
      corre <- cor(matrix)
      corretest <- corre[row.names(corre) %in% noncoding, colnames(corre) %in% coding]
      return(corre)
    }
  }
  
  CorreDataMMap <- reactive({
    CorreDataM(gpl = input$gplnumberexp, gse = input$gsenumberexp, noncoding = input$NONselectM, coding = input$COselectM)
  })
  
  output$corremap <- renderPlot({
    if (is.data.frame(CorreDataMMap())) {
      p <- ggplot(CorreDataMMap())
      p
    } else {
      data <- CorreDataMMap()
      pheatmap::pheatmap(data)
    }
    
  })
  
  CorreDataS <- function(gpl, gse, noncoding, coding) {
    if (noncoding == "" | coding == "") {
      return(data.frame())
    } else {
      sample <- eval(parse(text = paste(gpl, gse, "express_genetype", sep = "_")))
      non <- sample[which(sample$GeneSymbol == noncoding), 2:(ncol(sample) - 1)]
      co <- sample[which(sample$GeneSymbol == coding), 2:(ncol(sample) - 1)]
      data <- as.matrix(t(rbind(non, co)))
      colnames(data) <- c(noncoding, coding)
      return(data)
    }
  }
  
  CorreDataSplot <- reactive({
    CorreDataS(gpl = input$gplnumberexp, gse = input$gsenumberexp, noncoding = input$NONselectS, coding = input$COselectS)
  })
  
  output$correplot <- renderPlot({
    data <- CorreDataSplot()
    if (nrow(data) == 0) {
      p <- ggplot()
      p
    } else {
      plot(data, pch = 16, cex = 1)
      abline(lm(data[,2] ~ data[,1]))
      r <- cor(data)[colnames(data)[1], colnames(data)[2]]
      title(main = paste("Correlation Between", colnames(data)[1], "and", colnames(data)[2], sep = " "), sub = paste0("Correlation Coefficient: ", round(r, 5)))
    }
  })
  
  
  
  # ---------------------------------- The module for Shiny Forth Part ---------------------------- #
  # CVC, CVN, NVN means Canver versus Cancer, Cancer versus Normal and Normal versus Normal respectively
  # Functions in these three parts have same structure but differences in details
  
  # CVC Data Process return table data and heatmap data
  CVCDifferentialExpression <- function(sample1, sample2, n, ExpType, genetype) {
    
    if (sample1 == "--" | sample2 == "--") {
      return(list("datatable" = data.frame(), "HMdata" = data.frame()))
    }
    
    Sample1 <- eval(parse(text = paste(sample1, "_express_genetype", sep = "")))
    Sample2 <- eval(parse(text = paste(sample2, "_express_genetype", sep = "")))
    common.gene <- intersect(Sample1[, "GeneSymbol"], Sample2[, "GeneSymbol"])
    Sample1 <- Sample1[Sample1$GeneSymbol %in% common.gene, ]
    Sample2 <- Sample2[Sample2$GeneSymbol %in% common.gene, ]
    Sample <- cbind(Sample1[order(Sample1$GeneSymbol), c(1, ncol(Sample1))], Sample1[order(Sample1$GeneSymbol), 2:(ncol(Sample1) - 1)], Sample2[order(Sample2$GeneSymbol), 2:(ncol(Sample2) - 1)])
    SampleMatrix <- Sample[, 3:ncol(Sample)]
    row.names(SampleMatrix) <- Sample[, 1]
    SampleCategory <- rbind(subset(bc, bc$GPL_GSE == sample1, select = c(GSM, GPL_GSE)), subset(bc, bc$GPL_GSE == sample2, select = c(GSM, GPL_GSE)))
    
    # Process to follow the guide of using limma
    SampleCategoryFactor <- factor(SampleCategory[, "GPL_GSE"])
    design <- model.matrix(~-1+SampleCategoryFactor)
    contrast.matrix <- makeContrasts(contrasts = paste(colnames(design), collapse = "-"), levels = design)
    fit <- lmFit(SampleMatrix, design)
    fit1 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit1)
    dif <- topTable(fit2, coef = paste(colnames(design), collapse = "-"), n = nrow(fit2), lfc = log2(2.5))
    dif <- dif[dif[, "P.Value"] < 0.01, ]
    dif <- dif[order(row.names(dif)), ]
    dif$GeneType <- Sample$GeneType[Sample$GeneSymbol %in% row.names(dif)]
    dif <- cbind(GeneSymbol = row.names(dif), GeneType = dif[, ncol(dif)], dif[, 1:(ncol(dif) - 1)])
    
    if (ExpType == "High Expression") {
      difTopX <- ddply(dif[dif$logFC > 0, ], "GeneType", function(x) {x[order(-x$logFC)[1:n],]})
      difTopX <- na.omit(difTopX)
      HMSample <- Sample[Sample$GeneSymbol %in% difTopX$GeneSymbol, ]
      HMSample <- subset(HMSample, HMSample$GeneType == genetype)
      HMSample <- melt(data = HMSample, id.vars = c("GeneSymbol", "GeneType"))
      HMSample$Category <- rep(SampleCategory[, 2], each = length(unique(HMSample$GeneSymbol)))
      return(list("datatable" = difTopX, "HMdata" = HMSample))
    }
    
    if (ExpType == "Low Expression") {
      difTopX <- ddply(dif[dif$logFC < 0, ], "GeneType", function(x) {x[order(-abs(x$logFC))[1:n],]})
      difTopX <- na.omit(difTopX)
      HMSample <- Sample[Sample$GeneSymbol %in% difTopX$GeneSymbol, ]
      HMSample <- subset(HMSample, HMSample$GeneType == genetype)
      HMSample <- melt(data = HMSample, id.vars = c("GeneSymbol", "GeneType"))
      HMSample$Category <- rep(SampleCategory[, 2], each = length(unique(HMSample$GeneSymbol)))
      return(list("datatable" = difTopX, "HMdata" = HMSample))
    }
    
  }
  
  # For the DE table and heatmap
  CVCDifferentialExpressionData <- reactive({
    CVCDifferentialExpression(sample1 = input$CVC1_1, sample2 = input$CVC1_2, n = input$GeneNumber, ExpType = input$ExpType, genetype = input$GeneType)
  })
  
  output$CVCDETable <- renderDataTable(CVCDifferentialExpressionData()$datatable, rownames = F, 
                                       options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$CVCDEHeatmap <- renderPlot({
    HMData <- CVCDifferentialExpressionData()$HMdata
    if (nrow(HMData) == 0) {
      p <- ggplot(HMData)
      p
    } else {
      pos1x <- length(unique(HMData$variable[which(HMData$Category == unique(HMData$Category)[1])]))/2
      pos2x <- length(unique(HMData$variable))/2 + pos1x
      label <- unlist(strsplit(unique(HMData$Category), "_"))[c(2, 4)]
      p <- ggplot(HMData, aes(x = variable, y = GeneSymbol, fill = value)) + geom_tile() + labs(x = NULL, y = NULL) + geom_vline(aes(xintercept = 2 * pos1x), color = "white")
      p <- p + scale_fill_gradient(low = "green", high = "red", name = "Expression") + theme(axis.ticks = element_blank(), legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"), panel.background = element_blank(), panel.grid = element_blank(), axis.text.x = element_blank(), legend.title = element_text(size = 9), legend.position = "bottom")
      p + annotate("text", x = c(pos1x, pos2x), y = 1, label = label, size = 3)
    }
  })
  
  # build an renderUI to make the choices in selectInput change along with the sample
  output$selectUI <- renderUI({
    Expdata <- CVCDifferentialExpressionData()$datatable
    genes <- as.character(Expdata$GeneSymbol)
    selectInput("boxplotGeneInput", label = "Specific Gene", choices = as.list(setNames(c("", genes), c("", genes))), selected = "")
  })
  
  # boxplot data process
  CVCBoxplot <- function(sample1, sample2, geneinterest = "") {
    
    if (sample1 == "--" | sample2 == "--" | geneinterest == "") {
      return(data.frame())
    }
    
    Sample1 <- eval(parse(text = paste(sample1, "_express_genetype", sep = "")))
    Sample2 <- eval(parse(text = paste(sample2, "_express_genetype", sep = "")))
    common.gene <- intersect(Sample1[, "GeneSymbol"], Sample2[, "GeneSymbol"])
    Sample1 <- Sample1[Sample1$GeneSymbol %in% common.gene, ]
    Sample2 <- Sample2[Sample2$GeneSymbol %in% common.gene, ]
    Sample <- cbind(Sample1[order(Sample1$GeneSymbol), c(1, ncol(Sample1))], Sample1[order(Sample1$GeneSymbol), 2:(ncol(Sample1) - 1)], Sample2[order(Sample2$GeneSymbol), 2:(ncol(Sample2) - 1)])
    SampleCategory <- rbind(subset(bc, bc$GPL_GSE == sample1, select = c(GSM, GPL_GSE)), subset(bc, bc$GPL_GSE == sample2, select = c(GSM, GPL_GSE)))
    
    boxplotdata <- melt(Sample[which(Sample$GeneSymbol == geneinterest), ], id.vars = c("GeneSymbol", "GeneType"))
    boxplotdata$Category <- SampleCategory$GPL_GSE
    
    return(boxplotdata)
  }
  
  CVCBoxplotData <- reactive({
    CVCBoxplot(sample1 = input$CVC1_1, sample2 = input$CVC1_2, geneinterest = input$boxplotGeneInput)
  })
  
  # Draw boxplot
  output$CVCDEBoxplot <- renderPlot({
    Boxplotdata <- CVCBoxplotData()
    if (nrow(Boxplotdata) == 0) {
      pbox <- ggplot(Boxplotdata)
      pbox
    } else {
      pbox <- ggplot(Boxplotdata) + geom_boxplot(aes(x = Category, y = value), width = 0.3) + theme_bw() + ggtitle(label = unique(Boxplotdata$GeneSymbol))
      pbox
    }
  })
  
  # GO Analyze
  GOplotdata <- reactive({
    dependata <- CVCDifferentialExpression(sample1 = input$CVC1_1, sample2 = input$CVC1_2, n = input$GeneNumber, ExpType = input$ExpType, genetype = input$GeneType)$datatable
    if (nrow(dependata) == 0) {
      return(list("CC" = data.frame(),
                  "BP" = data.frame(),
                  "MF" = data.frame()))
    } else {
      symbol2entrez <- getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "hgnc_symbol", values = dependata$GeneSymbol, mart = mart)
      symbol2entrez <- symbol2entrez[match(dependata$GeneSymbol, symbol2entrez$hgnc_symbol), ]
      dependata$EntrezID <- symbol2entrez$entrezgene
      dependata <- na.omit(dependata)
      return(list("CC" = enrichGO(gene = dependata$EntrezID, ont = "CC", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "BP" = enrichGO(gene = dependata$EntrezID, ont = "BP", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "MF" = enrichGO(gene = dependata$EntrezID, ont = "MF", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "depend" = dependata))
    }
  })
  
  output$GOcnetplotcc <- renderPlot({
    data <- GOplotdata()$depend
    if (is.data.frame(GOplotdata()$CC) | nrow(summary(GOplotdata()$CC)) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(GOplotdata()$CC, categorySize = "pvalue", foldChange = data$logFC, fixed = T)
    }
  })
  
  output$GOtablecc <- renderDataTable(summary(GOplotdata()$CC), rownames = F,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$GOcnetplotmf <- renderPlot({
    data <- GOplotdata()$depend
    if (is.data.frame(GOplotdata()$MF) | nrow(summary(GOplotdata()$MF)) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(GOplotdata()$MF, categorySize = "pvalue", foldChange = data$logFC, fixed = T)
    }
  })
  
  output$GOtablemf <- renderDataTable(summary(GOplotdata()$MF), rownames = F,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$GOcnetplotbp <- renderPlot({
    data <- GOplotdata()$depend
    if (is.data.frame(GOplotdata()$BP) | nrow(summary(GOplotdata()$BP)) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(GOplotdata()$BP, categorySize = "pvalue", foldChange = data$logFC, fixed = T)
    }
  })
  
  output$GOtablebp <- renderDataTable(summary(GOplotdata()$BP), rownames = F,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  #KEGG Analyze
  KEGGdata <- reactive({
    Sample1 <- eval(parse(text = paste(input$CVC1_1, "_express_genetype", sep = "")))
    Sample2 <- eval(parse(text = paste(input$CVC1_2, "_express_genetype", sep = "")))
    common.gene <- intersect(Sample1[, "GeneSymbol"], Sample2[, "GeneSymbol"])
    Sample1 <- Sample1[Sample1$GeneSymbol %in% common.gene, ]
    Sample2 <- Sample2[Sample2$GeneSymbol %in% common.gene, ]
    Sample <- cbind(Sample1[order(Sample1$GeneSymbol), c(1, ncol(Sample1))], Sample1[order(Sample1$GeneSymbol), 2:(ncol(Sample1) - 1)], Sample2[order(Sample2$GeneSymbol), 2:(ncol(Sample2) - 1)])
    
    geneInput <- CVCDifferentialExpression(sample1 = input$CVC1_1, sample2 = input$CVC1_2, n = input$GeneNumber, ExpType = input$ExpType, genetype = input$GeneType)$datatable
    symbol2entrez <- getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "hgnc_symbol", values = geneInput$GeneSymbol, mart = mart)
    symbol2entrez <- symbol2entrez[match(geneInput$GeneSymbol, symbol2entrez$hgnc_symbol), ]
    geneInput$EntrezID <- symbol2entrez$entrezgene
    geneInput <- na.omit(geneInput)
    Expr <- Sample[, 3:ncol(Sample)]
    row.names(Expr) <- Sample$GeneSymbol
    Expr <- Expr[row.names(Expr) %in% geneInput$GeneSymbol, ]
    geneInput <- geneInput[, c("EntrezID", "logFC", "P.Value")]
    y <- geneAnswersBuilder(geneInput, "org.Hs.eg.db", categoryType = "KEGG", testType = "hyperG", pvalueT = 0.1, geneExpressionProfile = Expr, verbose = F)
    yy <- geneAnswersReadable(y, verbose = F)
    
    if (nrow(y@enrichmentInfo) == 0) {
      return(list("table" = data.frame(), "plot" = data.frame()))
    } else {
      return(list("table" = y@enrichmentInfo, "plot" = yy))
    }
  
  })
  
  output$KEGGtable <- renderDataTable(KEGGdata()$table, rownames = T,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$KEGGmap <- renderPlot({
    if (is.data.frame(KEGGdata()$plot)) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      geneAnswersConceptNet(KEGGdata()$plot, colorValueColumn = "logFC", centroidSize = "pvalue")
    }
    
  })
  
  
  
  # CVN Data Process return table data and heatmap data
  CVNDifferentialExpression <- function(sample, n, ExpType, genetype) {
    
    if (sample == "--") {
      return(list("datatable" = data.frame(), "HMdata" = data.frame()))
    }
    
    Sample <- eval(parse(text = paste(sample, "_express_genetype", sep = "")))
    Sample <- Sample[order(Sample$GeneSymbol), ]
    SampleMatrix <- Sample[, 2:(ncol(Sample) - 1)]
    row.names(SampleMatrix) <- Sample[, 1]
    SampleCategory <- subset(bc, bc$GPL_GSE == sample, select = c(GSM, Tumor.OR.Normal))
    SampleCategoryFactor <- factor(SampleCategory[, "Tumor.OR.Normal"])
    design <- model.matrix(~-1+SampleCategoryFactor)
    contrast.matrix <- makeContrasts(contrasts = paste(colnames(design), collapse = "-"), levels = design)
    fit <- lmFit(SampleMatrix, design)
    fit1 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit1)
    dif <- topTable(fit2, coef = paste(colnames(design), collapse = "-"), n = nrow(fit2), lfc = log2(2.5))
    dif <- dif[dif[, "P.Value"] < 0.05, ]
    dif <- dif[order(row.names(dif)), ]
    dif$GeneType <- Sample$GeneType[Sample$GeneSymbol %in% row.names(dif)]
    dif <- cbind(GeneSymbol = row.names(dif), GeneType = dif[, ncol(dif)], dif[, 1:(ncol(dif) - 1)])
    
    if (ExpType == "High Expression") {
      difTopX <- ddply(dif[dif$logFC > 0, ], "GeneType", function(x) {x[order(-x$logFC)[1:n],]})
      difTopX <- na.omit(difTopX)
      HMSample <- Sample[Sample$GeneSymbol %in% difTopX$GeneSymbol, ]
      HMSample <- subset(HMSample, HMSample$GeneType == genetype)
      HMSample <- melt(data = HMSample, id.vars = c("GeneSymbol", "GeneType"))
      HMSample$Category <- rep(SampleCategory[, 2], each = length(unique(HMSample$GeneSymbol)))
      return(list("datatable" = difTopX, "HMdata" = HMSample))
    }
    
    if (ExpType == "Low Expression") {
      difTopX <- ddply(dif[dif$logFC < 0, ], "GeneType", function(x) {x[order(-x$logFC)[1:n],]})
      difTopX <- na.omit(difTopX)
      HMSample <- Sample[Sample$GeneSymbol %in% difTopX$GeneSymbol, ]
      HMSample <- subset(HMSample, HMSample$GeneType == genetype)
      HMSample <- melt(data = HMSample, id.vars = c("GeneSymbol", "GeneType"))
      HMSample$Category <- rep(SampleCategory[, 2], each = length(unique(HMSample$GeneSymbol)))
      return(list("datatable" = difTopX, "HMdata" = HMSample))
    }
  }
  
  CVNDifferentialExpressionData <- reactive({
    CVNDifferentialExpression(sample = input$CVN, n = input$CVNGeneNumber, ExpType = input$CVNExpType, genetype = input$CVNGeneType)
  })
  
  output$CVNDETable <- renderDataTable(CVNDifferentialExpressionData()$datatable, rownames = F, 
                                       options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$CVNDEHeatmap <- renderPlot({
    HMData <- CVNDifferentialExpressionData()$HMdata
    if (nrow(HMData) == 0) {
      p <- ggplot(HMData)
      p
    } else {
      pos1x <- length(unique(HMData$variable[which(HMData$Category == unique(HMData$Category)[1])]))/2
      pos2x <- length(unique(HMData$variable))/2 + pos1x
      label <- unique(HMData$Category)
      p <- ggplot(HMData, aes(x = variable, y = GeneSymbol, fill = value)) + geom_tile() + labs(x = NULL, y = NULL) + geom_vline(aes(xintercept = (2 * pos1x) + 0.5), color = "white")
      p <- p + scale_fill_gradient(low = "green", high = "red", name = "Expression") + theme(axis.ticks = element_blank(), legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"), panel.background = element_blank(), panel.grid = element_blank(), axis.text.x = element_blank(), legend.title = element_text(size = 9), legend.position = "bottom")
      p + annotate("text", x = c(pos1x, pos2x), y = 1, label = label, size = 3)
    }
  })
  
  output$CVNselectUI <- renderUI({
    Expdata <- CVNDifferentialExpressionData()$datatable
    genes <- as.character(Expdata$GeneSymbol)
    selectInput("CVNboxplotGeneInput", label = "Specific Gene", choices = as.list(setNames(genes, genes)), selected = "")
  })
  
  CVNBoxplot <- function(sample, geneinterest) {
    
    if (sample == "--" | geneinterest == "") {
      return(data.frame())
    }
    
    Sample <- eval(parse(text = paste(sample, "_express_genetype", sep = "")))
    Sample <- Sample[order(Sample$GeneSymbol), ]
    SampleMatrix <- Sample[, 2:(ncol(Sample) - 1)]
    row.names(SampleMatrix) <- Sample[, 1]
    SampleCategory <- subset(bc, bc$GPL_GSE == sample, select = c(GSM, Tumor.OR.Normal))
    
    boxplotdata <- melt(Sample[which(Sample$GeneSymbol == geneinterest), ], id.vars = c("GeneSymbol", "GeneType"))
    boxplotdata$Category <- SampleCategory$Tumor.OR.Normal
    
    return(boxplotdata)
  }
  
  CVNBoxplotData <- reactive({
    CVNBoxplot(sample = input$CVN, geneinterest = input$CVNboxplotGeneInput)
  })
  
  output$CVNDEBoxplot <- renderPlot({
    Boxplotdata <- CVNBoxplotData()
    if (nrow(Boxplotdata) == 0) {
      pbox <- ggplot(Boxplotdata)
      pbox
    } else {
      pbox <- ggplot(Boxplotdata) + geom_boxplot(aes(x = Category, y = value), width = 0.3) + theme_bw() + ggtitle(label = unique(Boxplotdata$GeneSymbol))
      pbox
    }
  })
  
  CVNGOplotdata <- reactive({
    dependata <- CVNDifferentialExpression(sample = input$CVN, n = input$CVNGeneNumber, ExpType = input$CVNExpType, genetype = input$CVNGeneType)$datatable
    if (nrow(dependata) == 0) {
      return(list("CC" = data.frame(),
                  "BP" = data.frame(),
                  "MF" = data.frame()))
    } else {
      symbol2entrez <- getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "hgnc_symbol", values = dependata$GeneSymbol, mart = mart)
      symbol2entrez <- symbol2entrez[match(dependata$GeneSymbol, symbol2entrez$hgnc_symbol), ]
      dependata$EntrezID <- symbol2entrez$entrezgene
      dependata <- na.omit(dependata)
      return(list("CC" = enrichGO(gene = dependata$EntrezID, ont = "CC", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "BP" = enrichGO(gene = dependata$EntrezID, ont = "BP", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "MF" = enrichGO(gene = dependata$EntrezID, ont = "MF", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "depend" = dependata))
    }
  })
  
  output$CVNGOcnetplotcc <- renderPlot({
    data <- CVNGOplotdata()$depend
    if (is.data.frame(CVNGOplotdata()$CC) | nrow(summary(CVNGOplotdata()$CC)) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(CVNGOplotdata()$CC, categorySize = "pvalue", fixed = T)
    }
  })
  
  output$CVNGOtablecc <- renderDataTable(summary(CVNGOplotdata()$CC), rownames = F,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$CVNGOcnetplotmf <- renderPlot({
    data <- CVNGOplotdata()$depend
    if (is.data.frame(CVNGOplotdata()$MF) | nrow(summary(CVNGOplotdata()$MF)) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(CVNGOplotdata()$MF, categorySize = "pvalue", fixed = T)
    }
  })
  
  output$CVNGOtablemf <- renderDataTable(summary(CVNGOplotdata()$MF), rownames = F,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  
  output$CVNGOcnetplotbp <- renderPlot({
    data <- CVNGOplotdata()$depend
    if (is.data.frame(CVNGOplotdata()$BP) | nrow(summary(CVNGOplotdata()$BP)) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(CVNGOplotdata()$BP, categorySize = "pvalue", fixed = T)
    }
  })
  
  output$CVNGOtablebp <- renderDataTable(summary(CVNGOplotdata()$BP), rownames = F,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  CVNKEGGdata <- reactive({
    
    Sample <- eval(parse(text = paste(sample = input$CVN, "_express_genetype", sep = "")))
    Sample <- Sample[order(Sample$GeneSymbol), ]
    
    geneInput <- CVNDifferentialExpression(sample = input$CVN, n = input$CVNGeneNumber, ExpType = input$CVNExpType, genetype = input$CVNGeneType)$datatable
    symbol2entrez <- getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "hgnc_symbol", values = geneInput$GeneSymbol, mart = mart)
    symbol2entrez <- symbol2entrez[match(geneInput$GeneSymbol, symbol2entrez$hgnc_symbol), ]
    geneInput$EntrezID <- symbol2entrez$entrezgene
    geneInput <- na.omit(geneInput)
    Expr <- Sample[, 3:ncol(Sample)]
    row.names(Expr) <- Sample$GeneSymbol
    Expr <- Expr[row.names(Expr) %in% geneInput$GeneSymbol, ]
    geneInput <- geneInput[, c("EntrezID", "logFC", "P.Value")]
    y <- geneAnswersBuilder(geneInput, "org.Hs.eg.db", categoryType = "KEGG", testType = "hyperG", pvalueT = 0.1, geneExpressionProfile = Expr, verbose = F)
    yy <- geneAnswersReadable(y, verbose = F)
    
    if (nrow(y@enrichmentInfo) == 0) {
      return(list("table" = data.frame(), "plot" = data.frame()))
    } else {
      return(list("table" = y@enrichmentInfo, "plot" = yy))
    }
    
  })
  
  output$CVNKEGGtable <- renderDataTable(CVNKEGGdata()$table, rownames = T,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$CVNKEGGmap <- renderPlot({
    if (is.data.frame(CVNKEGGdata()$plot)) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      geneAnswersConceptNet(CVNKEGGdata()$plot, colorValueColumn = "logFC", centroidSize = "pvalue")
    }
  })
  
  #NVN Data Process
  
  NVNDifferentialExpression <- function(sample1, sample2, n, ExpType, genetype) {
    
    if (sample1 == "--" | sample2 == "--") {
      return(list("datatable" = data.frame(), "HMdata" = data.frame()))
    }
    
    Sample1 <- eval(parse(text = paste(sample1, "_express_genetype", sep = "")))
    Sample2 <- eval(parse(text = paste(sample2, "_express_genetype", sep = "")))
    common.gene <- intersect(Sample1[, "GeneSymbol"], Sample2[, "GeneSymbol"])
    Sample1 <- Sample1[Sample1$GeneSymbol %in% common.gene, ]
    Sample2 <- Sample2[Sample2$GeneSymbol %in% common.gene, ]
    Sample <- cbind(Sample1[order(Sample1$GeneSymbol), c(1, ncol(Sample1))], Sample1[order(Sample1$GeneSymbol), 2:(ncol(Sample1) - 1)], Sample2[order(Sample2$GeneSymbol), 2:(ncol(Sample2) - 1)])
    SampleMatrix <- Sample[, 3:ncol(Sample)]
    row.names(SampleMatrix) <- Sample[, 1]
    SampleCategory <- rbind(subset(bc, bc$GPL_GSE == sample1, select = c(GSM, GPL_GSE)), subset(bc, bc$GPL_GSE == sample2, select = c(GSM, GPL_GSE)))
    SampleCategoryFactor <- factor(SampleCategory[, "GPL_GSE"])
    design <- model.matrix(~-1+SampleCategoryFactor)
    contrast.matrix <- makeContrasts(contrasts = paste(colnames(design), collapse = "-"), levels = design)
    fit <- lmFit(SampleMatrix, design)
    fit1 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit1)
    dif <- topTable(fit2, coef = paste(colnames(design), collapse = "-"), n = nrow(fit2), lfc = log2(2.5))
    dif <- dif[dif[, "P.Value"] < 0.01, ]
    dif <- dif[order(row.names(dif)), ]
    dif$GeneType <- Sample$GeneType[Sample$GeneSymbol %in% row.names(dif)]
    dif <- cbind(GeneSymbol = row.names(dif), GeneType = dif[, ncol(dif)], dif[, 1:(ncol(dif) - 1)])
    
    if (ExpType == "High Expression") {
      difTopX <- ddply(dif[dif$logFC > 0, ], "GeneType", function(x) {x[order(-x$logFC)[1:n],]})
      difTopX <- na.omit(difTopX)
      HMSample <- Sample[Sample$GeneSymbol %in% difTopX$GeneSymbol, ]
      HMSample <- subset(HMSample, HMSample$GeneType == genetype)
      HMSample <- melt(data = HMSample, id.vars = c("GeneSymbol", "GeneType"))
      HMSample$Category <- rep(SampleCategory[, 2], each = length(unique(HMSample$GeneSymbol)))
      return(list("datatable" = difTopX, "HMdata" = HMSample))
    }
    
    if (ExpType == "Low Expression") {
      difTopX <- ddply(dif[dif$logFC < 0, ], "GeneType", function(x) {x[order(-abs(x$logFC))[1:n],]})
      difTopX <- na.omit(difTopX)
      HMSample <- Sample[Sample$GeneSymbol %in% difTopX$GeneSymbol, ]
      HMSample <- subset(HMSample, HMSample$GeneType == genetype)
      HMSample <- melt(data = HMSample, id.vars = c("GeneSymbol", "GeneType"))
      HMSample$Category <- rep(SampleCategory[, 2], each = length(unique(HMSample$GeneSymbol)))
      return(list("datatable" = difTopX, "HMdata" = HMSample))
    }
    
  }
  
  # For the DE table and heatmap
  NVNDifferentialExpressionData <- reactive({
    NVNDifferentialExpression(sample1 = input$NVN1_1, sample2 = input$NVN1_2, n = input$NVNGeneNumber, ExpType = input$NVNExpType, genetype = input$NVNGeneType)
  })
  
  output$NVNDETable <- renderDataTable(NVNDifferentialExpressionData()$datatable, rownames = F, 
                                       options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$NVNDEHeatmap <- renderPlot({
    HMData <- NVNDifferentialExpressionData()$HMdata
    if (nrow(HMData) == 0) {
      p <- ggplot(HMData)
      p
    } else {
      pos1x <- length(unique(HMData$variable[which(HMData$Category == unique(HMData$Category)[1])]))/2
      pos2x <- length(unique(HMData$variable))/2 + pos1x
      label <- unlist(strsplit(unique(HMData$Category), "_"))[c(2, 4)]
      p <- ggplot(HMData, aes(x = variable, y = GeneSymbol, fill = value)) + geom_tile() + labs(x = NULL, y = NULL) + geom_vline(aes(xintercept = 2 * pos1x), color = "white")
      p <- p + scale_fill_gradient(low = "green", high = "red", name = "Expression") + theme(axis.ticks = element_blank(), legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"), panel.background = element_blank(), panel.grid = element_blank(), axis.text.x = element_blank(), legend.title = element_text(size = 9), legend.position = "bottom")
      p + annotate("text", x = c(pos1x, pos2x), y = 1, label = label, size = 3)
    }
  })
  
  # build an renderUI to make the choices in selectInput change along with the sample
  output$NVNselectUI <- renderUI({
    Expdata <- NVNDifferentialExpressionData()$datatable
    genes <- as.character(Expdata$GeneSymbol)
    selectInput("NVNboxplotGeneInput", label = "Specific Gene", choices = as.list(setNames(genes, genes)), selected = "")
  })
  
  # boxplot data process
  NVNBoxplot <- function(sample1, sample2, geneinterest) {
    
    if (sample1 == "--" | sample2 == "--" | geneinterest == "") {
      return(data.frame())
    }
    
    Sample1 <- eval(parse(text = paste(sample1, "_express_genetype", sep = "")))
    Sample2 <- eval(parse(text = paste(sample2, "_express_genetype", sep = "")))
    common.gene <- intersect(Sample1[, "GeneSymbol"], Sample2[, "GeneSymbol"])
    Sample1 <- Sample1[Sample1$GeneSymbol %in% common.gene, ]
    Sample2 <- Sample2[Sample2$GeneSymbol %in% common.gene, ]
    Sample <- cbind(Sample1[order(Sample1$GeneSymbol), c(1, ncol(Sample1))], Sample1[order(Sample1$GeneSymbol), 2:(ncol(Sample1) - 1)], Sample2[order(Sample2$GeneSymbol), 2:(ncol(Sample2) - 1)])
    SampleCategory <- rbind(subset(bc, bc$GPL_GSE == sample1, select = c(GSM, GPL_GSE)), subset(bc, bc$GPL_GSE == sample2, select = c(GSM, GPL_GSE)))
    
    boxplotdata <- melt(Sample[which(Sample$GeneSymbol == geneinterest), ], id.vars = c("GeneSymbol", "GeneType"))
    boxplotdata$Category <- SampleCategory$GPL_GSE
    
    return(boxplotdata)
  }
  
  NVNBoxplotData <- reactive({
    NVNBoxplot(sample1 = input$NVN1_1, sample2 = input$NVN1_2, geneinterest = input$NVNboxplotGeneInput)
  })
  
  # Draw boxplot
  output$NVNDEBoxplot <- renderPlot({
    Boxplotdata <- NVNBoxplotData()
    if (nrow(Boxplotdata) == 0) {
      pbox <- ggplot(Boxplotdata)
      pbox
    } else {
      pbox <- ggplot(Boxplotdata) + geom_boxplot(aes(x = Category, y = value), width = 0.3) + theme_bw() + ggtitle(label = unique(Boxplotdata$GeneSymbol))
      pbox
    }
  })
  
  # GO Analyze
  NVNGOplotdata <- reactive({
    dependata <- NVNDifferentialExpression(sample1 = input$NVN1_1, sample2 = input$NVN1_2, n = input$NVNGeneNumber, ExpType = input$NVNExpType, genetype = input$NVNGeneType)$datatable
    if (nrow(dependata) == 0) {
      return(list("CC" = data.frame(),
                  "BP" = data.frame(),
                  "MF" = data.frame()))
    } else {
      symbol2entrez <- getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "hgnc_symbol", values = dependata$GeneSymbol, mart = mart)
      symbol2entrez <- symbol2entrez[match(dependata$GeneSymbol, symbol2entrez$hgnc_symbol), ]
      dependata$EntrezID <- symbol2entrez$entrezgene
      dependata <- na.omit(dependata)
      return(list("CC" = enrichGO(gene = dependata$EntrezID, ont = "CC", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "BP" = enrichGO(gene = dependata$EntrezID, ont = "BP", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "MF" = enrichGO(gene = dependata$EntrezID, ont = "MF", pvalueCutoff = 0.01, readable= TRUE, OrgDb = "org.Hs.eg.db"),
                  "depend" = dependata))
    }
  })
  
  output$NVNGOcnetplotcc <- renderPlot({
    data <- NVNGOplotdata()$depend
    if (is.data.frame(NVNGOplotdata()$CC) | nrow(summary(NVNGOplotdata()$CC)) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(NVNGOplotdata()$CC, categorySize = "pvalue", foldChange = data$logFC, fixed = T)
    }
  })
  
  output$NVNGOtablecc <- renderDataTable(summary(NVNGOplotdata()$CC), rownames = F,
                                         options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$NVNGOcnetplotmf <- renderPlot({
    data <- NVNGOplotdata()$depend
    if (is.data.frame(NVNGOplotdata()$MF) | nrow(summary(NVNGOplotdata()$MF)) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(NVNGOplotdata()$MF, categorySize = "pvalue", foldChange = data$logFC, fixed = T)
    }
  })
  
  output$NVNGOtablemf <- renderDataTable(summary(NVNGOplotdata()$MF), rownames = F,
                                         options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$NVNGOcnetplotbp <- renderPlot({
    data <- NVNGOplotdata()$depend
    if (is.data.frame(NVNGOplotdata()$BP) || nrow(NVNGOplotdata()$BP@result) == 0) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      cnetplot(NVNGOplotdata()$BP, categorySize = "pvalue", foldChange = data$logFC, fixed = T)
    }
  })
  
  output$NVNGOtablebp <- renderDataTable(summary(NVNGOplotdata()$BP), rownames = F,
                                         options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  #KEGG Analyze
  NVNKEGGdata <- reactive({
    Sample1 <- eval(parse(text = paste(input$NVN1_1, "_express_genetype", sep = "")))
    Sample2 <- eval(parse(text = paste(input$NVN1_2, "_express_genetype", sep = "")))
    common.gene <- intersect(Sample1[, "GeneSymbol"], Sample2[, "GeneSymbol"])
    Sample1 <- Sample1[Sample1$GeneSymbol %in% common.gene, ]
    Sample2 <- Sample2[Sample2$GeneSymbol %in% common.gene, ]
    Sample <- cbind(Sample1[order(Sample1$GeneSymbol), c(1, ncol(Sample1))], Sample1[order(Sample1$GeneSymbol), 2:(ncol(Sample1) - 1)], Sample2[order(Sample2$GeneSymbol), 2:(ncol(Sample2) - 1)])
    
    geneInput <- NVNDifferentialExpression(sample1 = input$NVN1_1, sample2 = input$NVN1_2, n = input$NVNGeneNumber, ExpType = input$NVNExpType, genetype = input$NVNGeneType)$datatable
    symbol2entrez <- getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "hgnc_symbol", values = geneInput$GeneSymbol, mart = mart)
    symbol2entrez <- symbol2entrez[match(geneInput$GeneSymbol, symbol2entrez$hgnc_symbol), ]
    geneInput$EntrezID <- symbol2entrez$entrezgene
    geneInput <- na.omit(geneInput)
    Expr <- Sample[, 3:ncol(Sample)]
    row.names(Expr) <- Sample$GeneSymbol
    Expr <- Expr[row.names(Expr) %in% geneInput$GeneSymbol, ]
    geneInput <- geneInput[, c("EntrezID", "logFC", "P.Value")]
    y <- geneAnswersBuilder(geneInput, "org.Hs.eg.db", categoryType = "KEGG", testType = "hyperG", pvalueT = 0.1, geneExpressionProfile = Expr, verbose = F)
    yy <- geneAnswersReadable(y, verbose = F)
    
    if (nrow(y@enrichmentInfo) == 0) {
      return(list("table" = data.frame(), "plot" = data.frame()))
    } else {
      return(list("table" = y@enrichmentInfo, "plot" = yy))
    }
    
  })
  
  output$NVNKEGGtable <- renderDataTable(NVNKEGGdata()$table, rownames = T,
                                      options = list(pageLength = 10, scrollX = "500px", scrollY = "190px"))
  
  output$NVNKEGGmap <- renderPlot({
    if (is.data.frame(NVNKEGGdata()$plot)) {
      p <- ggplot(data.frame()) + annotate("text", label = "No gene can be found included in this function", x = 0, y = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + labs(x = NULL, y = NULL)
      p
    } else {
      geneAnswersConceptNet(NVNKEGGdata()$plot, colorValueColumn = "logFC", centroidSize = "pvalue")
    }
    
  })
})
