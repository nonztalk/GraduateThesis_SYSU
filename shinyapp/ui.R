library(shinydashboard)
library(shiny)
library(DT)

header <- dashboardHeader(title = "LncBRCA", titleWidth = 250)

sidebar <- dashboardSidebar(
  width = 250,
  sidebarMenu(
    menuItem(tags$p("Basic Research Message Table", style = "font-size: 120%;"), 
             icon = icon("bars"), tabName = "basic"),
    menuItem(tags$p("Clinical Data", style = "font-size: 120%;"), 
             icon = icon("hospital-o"), tabName = "clinic"),
    menuItem(tags$p("Expression Matrix", style = "font-size: 120%;"), icon = icon("th"), tabName = "ExpData"),
    menuItem(tags$p("Different Expression", style = "font-size: 120%;"), icon = icon("th"), 
             menuSubItem("Cancer vs. Cancer", tabName = "CVC"),
             menuSubItem("Cancer vs. Normal", tabName = "CVN"),
             menuSubItem("Normal vs. Normal", tabName = "NVN"))
    
  )
)

body <- dashboardBody(
  tabItems(
# ----------------------- Basic Message Table ---------------------- #
    # The filter
    tabItem(
      "basic",
      box(
        width = 12,
        title = "Filter",
        status = "primary",
        solidHeader = T,
        collapsible = T,
        column(3,
               textInput("pubmedidbasic", label = "PubMed ID", value = "")
        ),
        column(3,
               textInput("gplnumberbasic", label = "GPL Number", value = "")
        ),
        column(3,
               textInput("gsenumberbasic", label = "GSE Number", value = "")
        ),
        column(3,
               selectInput("yearbasic", label = "Research Year", choices = c("--", 2003:2016), selected = "--")
        ),
        br(),
        submitButton("Confirm")
      ),
      # Show Data
      tabBox(
        width = 12,
        title = "BRCA Research Message",
        tabPanel("Table", 
                 dataTableOutput("df"),
                 downloadButton('downloadbasic', 'Download')),
        tabPanel("Summary", 
                 h3("Research Summary"),
                 textOutput("GSESummary"),
                 h3("Microarray Message"),
                 textOutput("GPLTitle"))
      )
    ),
# ----------------------- Clinical Data Table ---------------------- #
    # Filter box
    tabItem(
      "clinic",
      fluidRow(
        box(
          width = 12,
          title = "Filter",
          status = "primary",
          solidHeader = T,
          collapsible = T,
          fluidRow(
            column(4,
                  textInput("gplnumberclinic", label = "GPL Number", value = "")
            ),
            column(4,
                  textInput("gsenumberclinic", label = "GSE Number", value = "")
            ),
            column(4,
                  textInput("gsmnumberclinic", label = "GSM Number", value = "")
            ),
            column(6, 
                   selectInput("cliniccol", label = "Variables to show", choices = as.list(setNames(colnames(bcclinic), colnames(bcclinic))), selected = colnames(bcclinic)[c(1:3,17)], multiple = T))
          ),
          
          submitButton("Confirm")
        )
      ),
      # Data Exhibition
      fluidRow(
        box(
          title = "Clinical Data",
          width = 12,
          status = "primary",
          dataTableOutput("Cdf"),
          downloadButton('downloadclinic', 'Download')
        )
      )
    ),
# ----------------------- Expression Data Table ---------------------- #
    # Filter
    tabItem(
      "ExpData",
      fluidRow(
        box(
          title = "Series Choose",
          width = 12,
          status = "primary",
          solidHeader = T,
          collapsible = T,
          column(4,
               textInput("gplnumberexp", label = "GPL Number", value = "GPL8300")
          ),
          column(4,
               textInput("gsenumberexp", label = "GSE Number", value = "GSE1299")
          ),
          column(4,
               numericInput("TopXtoShow", label = "Gene to show", value = 10)
          ),
          h4("Correlation Map Control (Multiple Genes)"),
          column(6, uiOutput("COselectUIM")),
          column(6, uiOutput("NONselectUIM")),
          h4("Correlation Plot Control (Single Genes)"),
          column(6, uiOutput("COselectUIS")),
          column(6, uiOutput("NONselectUIS")),
          submitButton("Confirm")
        )
      ),
      
      # Show Noncoding Gene Table
      fluidRow(
        tabBox(
          title = "RMA Normalized Expression",
          width = 12,
          tabPanel("Expression Summary", dataTableOutput("Edf")),
          tabPanel("Expression Matrix", dataTableOutput("EdfMatrix"), downloadButton('downloadmatrix', 'Download'))
        )
      ),
      fluidRow(
        box(
          title = "Expression Heatmap (Coding Gene)",
          width = 6,
          status = "primary",
          plotOutput("HeatmapCoding")
        ),
        box(
          title = "Expression Heatmap (Noncoding Gene)",
          width = 6,
          status = "primary",
          plotOutput("HeatmapNoncoding")
        )
      ),
      fluidRow(
        tabBox(title = "Correlation between Coding and Noncoding Gene",
               width = 12,
               tabPanel("Correlatiom Map (Multiple genes)",
                        
                        plotOutput("corremap")),
               tabPanel("Correlation Figure (Single gene)",
                        
                        plotOutput("correplot"))
        )
      )
      
    ),

# ----------------------- Differential Expression ---------------------- #
# Cancer vs. Cancer
    tabItem(
      "CVC",
      fluidRow(
        box(
          title = "Settings",
          width = 12,
          status = "primary",
          solidHeader = T,
          collapsible = T,
          h4("Sample Select"),
          column(6, 
                 selectInput("CVC1_1", label = "First Sample", choices = as.list(setNames(c("--", tumor), c("--", tumor))), selected = "--")),
          column(6,
                 selectInput("CVC1_2", label = "Second Sample", choices = as.list(setNames(c("--", tumor), c("--", tumor))), selected = "--")),
          
          
          h4("Differential Expression Result Table Control"),
          column(6,
                 selectInput("ExpType", label = "Expression Type", choices = list("High Expression" = "High Expression", "Low Expression" = "Low Expression"), selected = "High Expression")),
          column(6,
                 numericInput("GeneNumber", label = "Top x to show", value = 10)),
          column(6,
                 selectInput("GeneType", label = "Heatmap GeneType", choices = list("Coding" = "Coding", "Noncoding" = "Noncoding"), selected = "Noncoding")),
          column(6,
                 uiOutput("selectUI")),
          
          submitButton("Confirm")
        )
      ),
      fluidRow(
        box(
          title = "Differential Expression Result Table",
          width = 12,
          status = "primary",
          dataTableOutput("CVCDETable")
        )
      ),
      fluidRow(
        box(
          title = "Differential Expression Result Heatmap",
          width = 6,
          status = "primary",
          plotOutput("CVCDEHeatmap")
        ),
        box(
          title = "Differential Expression Result Gene",
          width = 6,
          status = "primary",
          plotOutput("CVCDEBoxplot")
        )
      ),
      fluidRow(
        tabBox(
          title = "GO Analysis",
          width = 12,
          tabPanel("Cellular Component", 
                   dataTableOutput("GOtablecc"),
                   plotOutput("GOcnetplotcc")),
          tabPanel("Molecular Function", 
                   dataTableOutput("GOtablemf"),
                   plotOutput("GOcnetplotmf")),
          tabPanel("Biological Process", 
                   dataTableOutput("GOtablebp"),
                   plotOutput("GOcnetplotbp"))
        )
      ),
      fluidRow(
        tabBox(
          title = "KEGG Analyze",
          width = 12,
          tabPanel("Result table", dataTableOutput("KEGGtable")),
          tabPanel("Result map", plotOutput("KEGGmap"))
        )
      )
      
      
    ),
# Cancer vs. Normal
    tabItem(
      "CVN",
      fluidRow(
        box(
          title = "Settings",
          width = 12,
          status = "primary",
          solidHeader = T,
          collapsible = T,
          h4("Sample Select"),
          column(12, 
                 selectInput("CVN", label = "Sample", choices = as.list(setNames(c("--", tumornormal), c("--", tumornormal))), selected = "--")),
         
                 
          h4("Differential Expression Result Table Control"),
          column(6,
                 selectInput("CVNExpType", label = "Expression Type", choices = list("High Expression" = "High Expression", "Low Expression" = "Low Expression"), selected = "High Expression")),
          column(6,
                 numericInput("CVNGeneNumber", label = "Top x to show", value = 10)),
          column(6,
                 selectInput("CVNGeneType", label = "Heatmap GeneType", choices = list("Coding" = "Coding", "Noncoding" = "Noncoding"), selected = "Noncoding")),
          column(6,
                 uiOutput("CVNselectUI")),
          
          submitButton("Confirm")
        )
      ),
      fluidRow(
        box(
          title = "Differential Expression Result Table",
          width = 12,
          status = "primary",
          dataTableOutput("CVNDETable")
        )
      ),
      fluidRow(
        box(
          title = "Differential Expression Result Heatmap",
          width = 6,
          status = "primary",
          plotOutput("CVNDEHeatmap")
        ),
        box(
          title = "Differential Expression Result Gene",
          width = 6,
          status = "primary",
          plotOutput("CVNDEBoxplot")
        )
      ),
      fluidRow(
        tabBox(
          title = "GO Analysis",
          width = 12,
          tabPanel("Cellular Component", 
                   dataTableOutput("CVNGOtablecc"),
                   plotOutput("CVNGOcnetplotcc")),
          tabPanel("Molecular Function", 
                   dataTableOutput("CVNGOtablemf"),
                   plotOutput("CVNGOcnetplotmf")),
          tabPanel("Biological Process", 
                   dataTableOutput("CVNGOtablebp"),
                   plotOutput("CVNGOcnetplotbp"))
        )
      ),
      fluidRow(
        tabBox(
          title = "KEGG Analyze",
          width = 12,
          tabPanel("Result table", dataTableOutput("CVNKEGGtable")),
          tabPanel("Result map", plotOutput("CVNKEGGmap"))
        )
      )
      
    ),
# Normal vs. Normal
    tabItem(
      "NVN",
      fluidRow(
        box(
          title = "Settings",
          width = 12,
          status = "primary",
          solidHeader = T,
          collapsible = T,
          h4("Sample Select"),
          column(6, 
                 selectInput("NVN1_1", label = "First Sample", choices = as.list(setNames(c("--", normal), c("--", normal))), selected = "--")),
          column(6,
                 selectInput("NVN1_2", label = "Second Sample", choices = as.list(setNames(c("--", normal), c("--", normal))), selected = "--")),
          
          
          h4("Differential Expression Result Table Control"),
          column(6,
                 selectInput("NVNExpType", label = "Expression Type", choices = list("High Expression" = "High Expression", "Low Expression" = "Low Expression"), selected = "High Expression")),
          column(6,
                 numericInput("NVNGeneNumber", label = "Top x to show", value = 10)),
          column(6,
                 selectInput("NVNGeneType", label = "Heatmap GeneType", choices = list("Coding" = "Coding", "Noncoding" = "Noncoding"), selected = "Noncoding")),
          column(6,
                 uiOutput("NVNselectUI")),
          
          submitButton("Confirm")
        )
      ),
      fluidRow(
        box(
          title = "Differential Expression Result Table",
          width = 12,
          status = "primary",
          dataTableOutput("NVNDETable")
        )
      ),
      fluidRow(
        box(
          title = "Differential Expression Result Heatmap",
          width = 6,
          status = "primary",
          plotOutput("NVNDEHeatmap")
        ),
        box(
          title = "Differential Expression Result Gene",
          width = 6,
          status = "primary",
          plotOutput("NVNDEBoxplot")
        )
      ),
      fluidRow(
        tabBox(
          title = "GO Analysis",
          width = 12,
          tabPanel("Cellular Component", 
                   dataTableOutput("NVNGOtablecc"),
                   plotOutput("NVNGOcnetplotcc")),
          tabPanel("Molecular Function", 
                   dataTableOutput("NVNGOtablemf"),
                   plotOutput("NVNGOcnetplotmf")),
          tabPanel("Biological Process", 
                   dataTableOutput("NVNGOtablebp"),
                   plotOutput("NVNGOcnetplotbp"))
        )
      ),
      fluidRow(
        tabBox(
          title = "KEGG Analyze",
          width = 12,
          tabPanel("Result table", dataTableOutput("NVNKEGGtable")),
          tabPanel("Result map", plotOutput("NVNKEGGmap"))
        )
      )
    )
    
  )
)

dashboardPage(
  skin = "blue",
  header,
  sidebar,
  body
)
