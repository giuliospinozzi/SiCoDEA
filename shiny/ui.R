library(shiny)
library(shinyjs)

shinyUI(fluidPage(
    
    navbarPage("Drug Plates analysis",
               tabPanel("Single drug", 
                        fluidRow(
                            column(4,
                                   wellPanel(
                                       fluidRow(
                                           column(12,align="center",fileInput('file_s1','Choose CSV File for dose:', 
                                                                              accept=c(
                                                                                  "text/csv",
                                                                                  "text/comma-separated-values,text/plain",
                                                                                  ".csv")
                                           )),
                                           column(6,radioButtons('dec_s1', 'Decimal separator', 
                                                                 c("point"='.',"comma"=','),
                                                                 inline = F)),
                                           column(6,radioButtons('sep_s1', 'Column separator', 
                                                                 c("semicolon"=';',"comma"=',',"tab"="\t"),
                                                                 inline = F))
                                       ))),
                            column(4,
                                   wellPanel(
                                       fluidRow(
                                           column(12,align="center",
                                                  selectInput("norm_s",strong("Choose normalization method:"),
                                                              c("Max signal" = "max",
                                                                "Zero drug concentration" = "zero"))),
                                           column(12,align="center",
                                                  selectizeInput("drug_list", "Choose drug:", choices=NULL, selected=NULL)))),
                                       wellPanel(
                                           fluidRow(
                                           column(12,align="center",
                                                  downloadButton('downloadAll_s','Download Report',class="butt"),
                                                  tags$head(tags$style(".butt{background-color:#e9f4f7;} .butt{color: #337ab7;}")))
                                           
                                       ))),
                            column(4,
                                   wellPanel(
                                       fluidRow(
                                           column(12,align="center",fileInput('file_s2','Choose CSV File for response:', 
                                                                              accept=c(
                                                                                  "text/csv",
                                                                                  "text/comma-separated-values,text/plain",
                                                                                  ".csv")
                                           )),
                                           column(6,radioButtons('dec_s2', 'Decimal separator', 
                                                                 c("point"='.',"comma"=','),
                                                                 inline = F)),
                                           column(6,radioButtons('sep_s2', 'Column separator', 
                                                                 c("semicolon"=';',"comma"=',',"tab"="\t"),
                                                                 inline = F))
                                       ))),
                            column(8,plotOutput("drug_s1")),
                            column(4,align="center",br(),br(),br(),tableOutput("table_s1"),br(),
                                                  downloadButton('downloadAll_s1','Download Plot',class="butt"),
                                                  tags$head(tags$style(".butt{background-color:#e9f4f7;} .butt{color: #337ab7;}")))
                        )
               ),
               tabPanel("Drug comparison", 
                        sidebarLayout(
                            sidebarPanel(
                                shinyjs::useShinyjs(),
                                fluidRow(
                                column(12,align="center",fileInput('dose1','Choose CSV File for dose:', 
                                                                   accept=c(
                                                                       "text/csv",
                                                                       "text/comma-separated-values,text/plain",
                                                                       ".csv")
                                )),
                                column(6,radioButtons('dec_ss0', 'Decimal separator', 
                                                      c("point"='.',"comma"=','),
                                                      inline = F)),
                                column(6,radioButtons('sep_ss0', 'Column separator', 
                                                      c("semicolon"=';',"comma"=',',"tab"="\t"),
                                                      inline = F))
                                
                                ),
                                hr(),
                                fluidRow(column(6, align="center",
                                                checkboxGroupInput("cell_show", "Cell lines to show:",
                                                                   c("Cell Line 1" = "cell_sh1",
                                                                     "Cell Line 2" = "cell_sh2",
                                                                     "Cell Line 3" = "cell_sh3",
                                                                     "Cell Line 4" = "cell_sh4")),
                                #numericInput("cell_num", strong("Number of cell lines (up to 4):"), 1,min=1,max=4)
                                ),
                                column(6,radioButtons('rep_no', 'Plot representation', 
                                                                     c("Replicates"='rep',"Mean"='mea'),
                                                                     inline = F))
                                ),
                                #selectInput("cell_num",strong("Number of cell lines:"),c("1","2","3","4")),
                                conditionalPanel(condition="input.cell_show.includes('cell_sh1')",
                                                 
                                hr(),
                                fluidRow(
                                    column(12,align="center",textInput('name1',strong("Cell Line 1 name:"),'Cell Line 1')),
                                    column(12,align="center",fileInput('cell1','Choose CSV file for Cell Line 1:', 
                                                                       accept=c(
                                                                           "text/csv",
                                                                           "text/comma-separated-values,text/plain",
                                                                           ".csv")
                                    )),
                                    column(6,radioButtons('dec_ss1', 'Decimal separator', 
                                                          c("point"='.',"comma"=','),
                                                          inline = F)),
                                    column(6,radioButtons('sep_ss1', 'Column separator', 
                                                          c("semicolon"=';',"comma"=',',"tab"="\t"),
                                                          inline = F))
                                )
                                ),
                                conditionalPanel(condition="input.cell_show.includes('cell_sh2')",
                                                 hr(),
                                                 fluidRow(
                                                     column(12,align="center",textInput('name2',strong("Cell Line 2 name:"),'Cell Line 2')),
                                                     column(12,align="center",fileInput('cell2','Choose CSV file for Cell Line 2:', 
                                                                                        accept=c(
                                                                                            "text/csv",
                                                                                            "text/comma-separated-values,text/plain",
                                                                                            ".csv")
                                                     )),
                                                     column(6,radioButtons('dec_ss2', 'Decimal separator', 
                                                                           c("point"='.',"comma"=','),
                                                                           inline = F)),
                                                     column(6,radioButtons('sep_ss2', 'Column separator', 
                                                                           c("semicolon"=';',"comma"=',',"tab"="\t"),
                                                                           inline = F))
                                                     
                                                 )
                                ),
                                conditionalPanel(condition="input.cell_show.includes('cell_sh3')",
                                                 hr(),
                                                 fluidRow(
                                                     column(12,align="center",textInput('name3',strong("Cell Line 3 name:"),'Cell Line 3')),
                                                     column(12,align="center",fileInput('cell3','Choose CSV file for Cell Line 3:', 
                                                                                        accept=c(
                                                                                            "text/csv",
                                                                                            "text/comma-separated-values,text/plain",
                                                                                            ".csv")
                                                     )),
                                                     column(6,radioButtons('dec_ss3', 'Decimal separator', 
                                                                           c("point"='.',"comma"=','),
                                                                           inline = F)),
                                                     column(6,radioButtons('sep_ss3', 'Column separator', 
                                                                           c("semicolon"=';',"comma"=',',"tab"="\t"),
                                                                           inline = F))
                                                     
                                                 )
                                ),
                                conditionalPanel(condition="input.cell_show.includes('cell_sh4')",
                                                 hr(),
                                                 fluidRow(
                                                     column(12,align="center",textInput('name4',strong("Cell Line 4 name:"),'Cell Line 4')),
                                                     column(12,align="center",fileInput('cell4','Choose CSV file for Cell Line 4:', 
                                                                                        accept=c(
                                                                                            "text/csv",
                                                                                            "text/comma-separated-values,text/plain",
                                                                                            ".csv")
                                                     )),
                                                     column(6,radioButtons('dec_ss4', 'Decimal separator', 
                                                                           c("point"='.',"comma"=','),
                                                                           inline = F)),
                                                     column(6,radioButtons('sep_ss4', 'Column separator', 
                                                                           c("semicolon"=';',"comma"=',',"tab"="\t"),
                                                                           inline = F))
                                                     
                                                 )
                                )

                            ),
                            mainPanel(
                                wellPanel(
                                fluidRow(
                                column(4,align="center",
                                       selectInput("norm_ss",strong("Choose normalization method:"),
                                                   c("Max signal" = "max",
                                                     "Zero drug concentration" = "zero"))),
                                column(4,align="center",
                                       selectizeInput("drug_list_s", "Choose drug:", choices=NULL, selected=NULL)),
                                column(4,align="center",
                                           selectInput("drug_ss",strong("Choose model to use:"),
                                                       c("Log-logistic"="log_log",
                                                         "Log-logistic [0]"="log_log_0",
                                                         "Log-logistic [01]"="log_log_01", 
                                                         "Log-logistic [01] v2"="log_log_01_v2",
                                                         "Log-logistic [1]"="log_log_1",
                                                         "Median-effect"="med_eff"))
                                           ))
                            ),
                            fluidRow(
                                column(12,align="center",plotOutput("drug_cell")),
                                column(12,align="center",downloadButton('do','Download Plot'))),
                            #        tags$head(tags$style(".butt{background-color:#e9f4f7;} .butt{color: #337ab7;}")))
                            br(),
#                            wellPanel(
                            fluidRow(
                                # column(8,align="center",
                                #             textInput('out_dir',strong("Project Name (no spaces):"),paste0("Project_",Sys.Date()))),
                                column(12, align="center", br(),
                                            downloadButton('downloadAll_ss1','Download all plots',class="butt"),
                                            tags$head(tags$style(".butt{background-color:#e9f4f7;} .butt{color: #337ab7;}"))))
                            # actionButton("report","Generate report")))
                            
#                            )
                            )
                        )),
               tabPanel("Drug combination",fluidRow(
                   column(2),
                   column(4,
               wellPanel(
                         fluidRow(
                             column(12,align="center",fileInput('file1','Choose CSV File:', 
                             accept=c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                   )),
                   column(6,radioButtons('dec', 'Decimal separator', 
                                         c("point"='.',"comma"=','),
                                inline = F)),
                   column(6,radioButtons('sep', 'Column separator', 
                                         c("semicolon"=';',"comma"=',',"tab"="\t"),
                                inline = F))
                   ))),
        column(4,
               wellPanel(
                   fluidRow(
                   column(12,align="center",
                         selectInput("norm",strong("Choose normalization method:"),
                                     c("Max signal" = "max",
                                       "Zero drug concentration" = "zero"))),
                   column(9,textInput('name1',strong("Drug 1 name:"),"")),
                   column(3,br(),actionButton("go1", "Change")),
                   column(9,textInput('name2',strong("Drug 2 name:"),"")),
                   column(3,br(),actionButton("go2", "Change"))
               ))),
        column(12,align="center",
               downloadButton('downloadAll','Download Report',class="butt"),
               tags$head(tags$style(".butt{background-color:#e9f4f7;} .butt{color: #337ab7;}"))),
        column(12,align="center",
               br(),
               tableOutput("table")),
        column(12,
               tabsetPanel(
                   tabPanel("Single drug",
                            fluidRow(
                                column(6,
                                       wellPanel(align="center",
                                           plotOutput("drug1"),
                                           br(),
                                           downloadButton('downloadPlot1','Download Plot')
                                           )
                                       ),
                                column(6,
                                       wellPanel(align="center",
                                                 plotOutput("drug2"),
                                                 br(),
                                                 downloadButton('downloadPlot2','Download Plot')
                                       )
                                       ),
                                column(3),
                                column(6,
                                       wellPanel(align="center",
                                           plotOutput("r2_plot"),
                                           br(),
                                           downloadButton('downloadPlot','Download Plot')
                                       ))),
                   ),
                   tabPanel("Combination Index",
                            fluidRow(
                                column(2),
                                column(4,align="center",
                                       wellPanel(
                                       selectInput("drug",strong("Choose model to use:"),
                                        c("Log-logistic",
                                          "Log-logistic [0]",
                                          "Log-logistic [01]", 
                                          "Log-logistic [01] v2",
                                          "Log-logistic [1]",
                                          "Median-effect")),
                                       uiOutput("formula"),
                                       uiOutput("formula1"))),
                                column(4,wellPanel(tableOutput("table_r2"))),
                                column(2),
                                br()),
                            fluidRow(
                            column(6,
                                   wellPanel(align="center",
                                             plotOutput("plot1a"),
                                             br(),
                                             downloadButton('downloadPlot1a','Download Plot')
                                   )),
                            column(6,
                                   wellPanel(align="center",
                                             plotOutput("plot1b"),
                                             br(),
                                             downloadButton('downloadPlot1b','Download Plot')
                                   ))
                            ))
               ))
    )
               ))
))