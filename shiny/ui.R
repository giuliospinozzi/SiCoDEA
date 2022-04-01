library(shiny)
library(shinyjs)

shinyUI(fluidPage(
  
  navbarPage("Drug Plates analysis",
             tabPanel("Single drug", 
                      fluidRow(
                        column(3,
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
                                                         inline = F)),
                                   br(),
                                   column(12,align="center",downloadButton('ex1','Download Example'))
                                 ))),
                        column(6,
                               wellPanel(
                                 fluidRow(
                                   column(6,align="center",
                                          selectInput("norm_s",strong("Choose normalization method:"),
                                                      c("Max signal" = "max",
                                                        "Zero drug concentration" = "zero",
                                                        "No normalization" = "none"))),
                                   column(6,align="center",
                                          numericInput("out_pv_s",strong("Choose outlier p-value 
                                                       (0 is for no outliers,
                                                       1 remove one replicate per dose concentration):"),
                                                       0,min=0, max=1,step = 0.01)),
                                   column(6,align="center",
                                          selectInput("in_via_s",strong("Inhibition/Viability:"),
                                                      c("Inhibition" = "in",
                                                        "Viability" = "via"))),
                                   column(6,align="center",
                                          selectizeInput("drug_list", "Choose drug:", choices=NULL, selected=NULL)))),
                               wellPanel(
                                 fluidRow(
                                   column(12,align="center",
                                          downloadButton('downloadAll_s','Download Report',class="butt"),
                                          tags$head(tags$style(".butt{background-color:#e9f4f7;} .butt{color: #337ab7;}")))
                                   
                                 ))),
                        column(3,
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
                                                         inline = F)),
                                   br(),
                                   column(12,align="center",downloadButton('ex2','Download Example'))
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
                              column(6,align="center",
                                     selectInput("norm_ss",strong("Choose normalization method:"),
                                                 c("Max signal" = "max",
                                                   "Zero drug concentration" = "zero",
                                                   "No normalization" = "none"))),
                              column(6,align="center",
                                     numericInput("out_pv_ss",strong("Choose outlier p-value 
                                                       (0 is for no outliers,
                                                       1 remove one replicate per dose concentration):"),
                                                  0,min=0, max=1,step = 0.01))),
                            fluidRow(
                              column(4,align="center",
                                     selectizeInput("drug_list_s", "Choose drug:", choices=NULL, selected=NULL)),
                              column(4,align="center",
                                     selectInput("in_via_s1",strong("Inhibition/Viability:"),
                                                 c("Inhibition" = "in",
                                                   "Viability" = "via"))),
                              column(4,align="center",
                                     selectInput("drug_ss",strong("Choose model to use:"),
                                                 c("Log-logistic"="log_log",
                                                   "Log-logistic [0]"="log_log_0",
                                                   "Log-logistic [01]"="log_log_01",
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
               column(4,
                      wellPanel(
                        fluidRow(
                          column(4,align="center",
                                 selectInput("mat_tab",strong("Input file:"),
                                             c("Matrix" = "mat",
                                               "Table" = "tab"))),
                          column(8,align="center",fileInput('file1','Choose CSV File:', 
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
                                                inline = F)),
                          br(),
                          column(6,align="center",downloadButton('ex3','Matrix Example')),
                          column(6,align="center",downloadButton('ex4','Table Example'))
                        ))),
               column(8,
                      wellPanel(
                        fluidRow(
                          column(4,align="center",
                                 selectInput("norm",strong("Choose normalization method:"),
                                             c("Max signal" = "max",
                                               "Zero drug concentration" = "zero",
                                               "No normalization" = "none"))),
                          column(4,align="center",
                                 selectInput("in_via",strong("Inhibition/Viability:"),
                                             c("Inhibition" = "in",
                                               "Viability" = "via"))),
                          column(4,align="center",
                                 numericInput("out_pv",strong("Choose outlier p-value
                                                       (0 is for no outliers,
                                                       1 remove one replicate per dose concentration):"),
                                              0,min=0, max=1,step = 0.01)),
                          column(4,textInput('name_drug1',strong("Drug 1 name:"),NULL)),
                          column(4,textInput('name_drug2',strong("Drug 2 name:"),NULL))
                        ))),
               column(12,align="center",
                      downloadButton('downloadAll','Download Report',class="butt"),
                      tags$head(tags$style(".butt{background-color:#e9f4f7;} .butt{color: #337ab7;}"))),
               column(2),
               column(8,align="center",
                      br(),
                      plotOutput("table")),
               column(12,align="center",downloadButton('downloadPlot_tab','Download Plot')),
               br(),
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
                                            selectInput("ci_model",strong("Choose combination index model to use:"),
                                                        c("Loewe additivity model"="Loewe",
                                                          "Response additivity model"="Response Additivity",
                                                          "Highest single agent model (HSA)"="HSA", 
                                                          "Bliss independence model"="Bliss",
                                                          "Zero Interaction Potency model (ZIP)"="ZIP")),
                                            uiOutput("formula0"),
                                            conditionalPanel(condition="input.ci_model == 'Loewe' || input.ci_model == 'ZIP'",
                                                             selectInput("drug",strong("Choose drug-response curve model to use:"),
                                                                         c("Log-logistic",
                                                                           "Log-logistic [0]",
                                                                           "Log-logistic [01]", 
                                                                           "Log-logistic [1]",
                                                                           "Median-effect")),
                                                             uiOutput("formula"),
                                                             uiOutput("formula1"))
                                          )),
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
             ),
             tabPanel("READ ME",h1("SiCoDEA"),
                      p("SiCoDEA (Single and Combined Drug Effect Analysis) is a simple, fast and 
           complete app for analyzing the effect of individual drugs and their 
           combinations and is available at:",
                        tags$a(href="https://sicodea.shinyapps.io/shiny/",
                               "https://sicodea.shinyapps.io/shiny/"),
                        ". A video tutorial is also available to understand how it works: ",
                        tags$a(href="https://youtu.be/kzzdU83r0_I",
                               "https://youtu.be/kzzdU83r0_I"),
                        ". The app consists of three panels, each of which allows you to carry 
           out a different analysis: a single drug analysis,
           a comparison of different drugs or an analysis of drug combinations."),
                      h2("Single drug analysis"),
                      img(src = "screenshot1.png",height = 650, width = 1000),
                      p("It can be done in the first tab, called “Single drug”. The data relating
           to drug doses and responses must be uploaded to the two side panels. To
           get an idea of the required formatting, you can download the sample files.
           For the dose file you need a matrix in which the first column contains 
           the names of the different drugs under examination and in which each row
           contains the sequence of doses used for the relative drug. For the second
           file, however, a matrix with the response values calculated for the 
           different doses is required. The sequence and position in the response 
           file must be corresponding to the dose values represented in the other 
           file. Both files must be in CSV format, but you can choose whether to 
           use comma, semicolon or tab separators."),
                      p("Once the files have been uploaded, it is then possible to choose 
         different options in the central panel. Normalization can be performed 
         considering the maximum signal of the entire plate as a basis, or 
         considering the average of the replicates obtained for drug concentrations
         equal to zero as the baseline. If, on the other hand, the normalization 
         has already been carried out on the loaded data, you can choose the “no 
         normalization” option. Based on the experiment in question, it is then 
         possible to choose between inhibition and viability to indicate the 
         information provided by the response values loaded. A Grubbs’ test is 
         performed on the replicates to assess the presence of any outliers and 
         the filter applied is based on the calculated p-value. It can be set in 
         the central panel using values ranging from 0 to 1, where in the first 
         case replicas are not eliminated while in the second one is eliminated 
         for each dose concentration. Rather than setting a fixed value, we 
         preferred to leave the possibility of choice in order to evaluate in 
         real time and visually from the plot the actual nature of an outlier. 
         Finally, there is a drop-down menu that allows you to navigate between 
         the different drugs present in the input file."),
                      p("The plot below shows, for each drug, the individual values at the 
         different doses as points, while the curves represent the expected trend 
         of the five different models examined. For each, the value of the IC50 
         is also shown with a dashed vertical line. The next table summarizes 
         the value of each IC50, as well as the value of",
                        HTML(paste("R", tags$sup(2), sep = "")),"for each model, so 
         as to be able to evaluate the most suitable for the data being analyzed. 
         The closer to 1 the value of",HTML(paste("R", tags$sup(2), sep = "")),
                        "the better the model for the data under consideration."),
                      p("Finally, there is a “Download plot” button that allows you to download
         the single plot that is displayed at that moment, or the “Download report”
         button, which generates a report containing all the plots for the various
         drugs."),
                      h2("Comparison of different drugs"),
                      img(src = "screenshot2.png",height = 600, width = 1000),
                      p("The second tab concerns the comparison between drugs. Again, there is a 
           side panel for uploading files. The formatting required for doses and 
           responses is identical to that of the previous tab. In this case, however,
           it is possible to load up to four files for the response values, while 
           the dose file remains the same; for example, these may be four different
           cell lines to which the same sequence of doses of a drug is applied to 
           evaluate its effectiveness."),
                      p("It is possible to choose which and how many cell lines to load and show
           by selecting the different ticks and with each new tick selected a new 
           panel will appear to load the file. You can also choose whether to 
           represent replicates as distinct points in the plot or as mean and 
           standard deviation."),
                      p("In this case also, as in the previous one, the names of the drugs are 
           those reported in the input file, while the names of the cell lines can
           be assigned in the specific spaces."),
                      p("The options in the central panel are the same as those seen in the 
           previous tab. However, we also find the option relating to the model 
           to be used. In this case, in fact, in the plot the different cell lines
           are compared and not the different models and therefore a specific model
           must be chosen to be represented also, possibly, on the basis of the",
                        HTML(paste("R", tags$sup(2), sep = "")),"calculated in the previous tab."),
                      h2("Analysis of drug combinations"),
                      img(src = "screenshot3.png",height = 600, width = 1000),
                      p("Finally, in the third tab we have the analysis of drug combinations. 
           It is possible to choose in which format to load the input file, whether 
           as a table, showing a column with the doses of a drug, a column with the
           doses of the other drug and a column with the recorded values, or as a
           matrix, reporting the doses of a drug in the first column, those of the 
           other drug in the first row and the values recorded in the central cells.
           There must also be doses of one drug and the other equal to zero in order
           to calculate the baseline value and build the curve of the individual 
           drugs. For clarity, there is a downloadable example for each format."),
                      p("The options for the normalization method, for the choice between 
           viability and inhibition, for the p-value of the outliers are the same 
           as already seen. In addition, it is also possible to change the name of
           the drugs directly from the interface."),
                      p("The first plot shows a matrix with the average of the replicates present
           in the source file for each dose combination and there is the usual 
           button to download it. Below then we have two tabs, one for the plots 
           relating to the individual drugs and the other for the combination."),
                      img(src = "screenshot4.png",height = 600, width = 1000),
                      p("In the first tab we have the plots for the dose-response curve of the
           two drugs under examination and the representation of all five models. 
           This is the same plot also seen in the “Single Drug” tab. In this case 
           also, the",HTML(paste("R", tags$sup(2), sep = "")),"value is calculated
           for each of the five models and in both drugs. The results are shown in 
           the last plot, in which there are the points of intersection of the",
                        HTML(paste("R", tags$sup(2), sep = "")),"for one drug and the other. 
           The points that are most located at the top right are those with the most
           suitable model for the data under examination. This way we can choose 
           the right model for our experiment."),
                      img(src = "screenshot5.png",height = 650, width = 1000),
                      p("In the second tab it is possible to choose the model for the combination
           index from five different options:",em("Response Additivity"),"model,",
                        em("Highest Single Agent"),"(HSA) model,",em("Bliss Independence"),
                        "model,",em("Loewe Additivity"),"model and",em("Zero Interaction Potency"),
                        "(ZIP) model. By selecting the different models from the drop-down menu,
           the corresponding formula will appear in the panel. In the case of the 
           Loewe Additivity model and ZIP model, since they use an Effect-Based 
           Strategy, it is also necessary to choose the model for the dose-response
           curve."),
                      p("The next table shows the two",HTML(paste("R", tags$sup(2), sep = "")),
                        "values for the five models, the same ones also present in the previous
           plot. The two plots below show the values of the combination index for 
           each combination of doses under examination based on the models chosen;
           blue represents synergy, red antagonism and black additivity. These are
           two different representations of the same results: in the first case 
           there are dots with dimensions proportional to the strength of synergy 
           or antagonism, while in the second case we have a matrix with colors of
           intensity proportional to the values of combination index reported within
           the cells. You can then browse through the various models to see how the
           results vary for the same data. Finally, you can download the single 
           plots or the report with all the plots.")
             )
  )
))