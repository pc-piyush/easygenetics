#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)

header<-dashboardHeader(title = HTML(paste(icon('dna'),'Easy Genetics')) ,
                        dropdownMenuOutput("messageMenu")
)


body <- dashboardBody(
  tabItems(

    # First tab content
    tabItem(tabName = "Simulate_Genotypes",
            titlePanel("Simulating Genotypes at a SNP over Multiple Generations"),
            fluidRow(
     box(
            numericInput("initial.n", "Initial Sample Size", 500),
            numericInput("Rate", "Mutation Rate", 0.001),
            actionButton("go", "Run")
     ), box(
            numericInput("Growth", "Average Number of Offspring Per Generation", 2.15),
            numericInput("Generations", "Number of Generations", 50),
            numericInput("seed2", "Enter Seed", NULL),
            checkboxInput("setseed", "Set Seed", F)
     )
            ),
     fluidRow(
     box(
     plotOutput(outputId = "plot1")
     ),

     box(
       plotOutput(outputId = "plot2")
     )

    )
    ),

    tabItem(tabName = "Haseman-Elston_Regression",
            # Application title
            titlePanel("Haseman-Elston Regression"),

            # Sidebar with a slider input for number of bins
            sidebarLayout(
              sidebarPanel(

                actionButton("runhe", "Run", width='100%'),

                selectInput("plothe", "Type of Plot",
                            choices = c("All", "Histogram of IBD", "IBD vs Differences in Sibling Phenotypes","Phenotype- Sibling 1 vs Sibling 2","Phenotype IBD=0 Sibling 1 vs Sibling 2", "Phenotype IBD=1 Sibling 1 vs Sibling 2", "Phenotype IBD=2 Sibling 1 vs Sibling 2")),

                numericInput("NumFamhe",
                             "Enter Number of Families",
                             value = 100),
                #num.families,var.total,var.genetic,var.env,maf,seed,row,column
                numericInput("TotVarhe",
                             "Total Variance",
                             value = 4),

                numericInput("GenVarhe",
                             "Genetic Variance",
                             value = 2),

                numericInput("EnvVarhe",
                             "Environmental Variance",
                             value = 1),

                numericInput("mafhe",
                             "Minor Allele Frequency",
                             min=0, max=0.5,
                             value = 0.5),

                numericInput("seedhe",
                             "Enter Seed",
                             NULL),

                checkboxInput("setseedhe", "Set Seed", F)
              ),

              # Show a plot of the generated distribution
              mainPanel(
                plotOutput("HE")
              )

            )

            ),

tabItem(tabName = "Breeding_Values",
        titlePanel("Genotypic Values, Breeding Values, and Dominance Deviation Table"),
        #sidebarLayout(
        #  sidebarPanel(
        fluidRow(

        box(
          title = "Observations of Phenotype Means", width=4,
        numericInput("A1A1", "A1A1", 140),
        numericInput("A1A2", "A1A2", 127),
        numericInput("A2A2", "A2A2", 79)

        ),


        box(width=4,
            numericInput("p", "Frequency of Allele A1", min=0, max=1, 0.3),
          numericInput("q", "Frequency of Allele A2", min=0, max=1, 0.7)

        ),

        box(width=4,
          title = "Other Estimates",
          #textOutput("q"),
          #textOutput("a"),
          tableOutput("estimates")
        )

        ),
        #mainPanel(
        fluidRow(
          column(12,
        box(width = NULL,
          tableOutput("bvplot"))
        )
        )
          ),

tabItem(tabName = "Multiple_QTL_Mapping",
        # Application title
        titlePanel("Multiple QTL Mapping by Haley-Knott Regression"),

        # Sidebar with a slider input for number of bins
        sidebarLayout(
          sidebarPanel(
            actionButton("runqtl", "Run", width = "100%"),
            selectInput("plotqtl", "Type of Function",
                        choices = c("Genetic Map",
                                    "Single QTL",
                                    "Two QTLs",
                                    "Heat Map-Two QTLs Comparison",
                                    "Two QTLs Comparison",
                                    "Composite Interval Mapping")),

            numericInput("len",
                         "Length of Chromosome",
                         value = 100),
            #num.families,var.total,var.genetic,var.env,maf,seed,row,column
            numericInput("nmark",
                         "Number of Markers",
                         value = 10),

            numericInput("ind",
                         "Number of Individuals",
                         value = 100),


            selectInput("space",
                        "Equal Spacing",
                        choices = c("TRUE","FALSE")),

            selectInput("type",
                        "Type of Cross",
                        choices = c("bc","f2")),


            numericInput("pos1",
                         "QTL 1- Position",
                         value = 10),

            numericInput("eff1",
                         "QTL 1-Effect",
                         value = 5),

            numericInput("pos2",
                         "QTL 2- Position",
                         value = 20),

            numericInput("eff2",
                         "QTL 2-Effect",
                         value = -5),

            numericInput("seedqtl",
                         "Enter Seed",
                         NULL),

            checkboxInput("setseedqtl", "Set Seed", F)


          ),

          # Show a plot of the generated distribution
          mainPanel(
            plotOutput("mQTLP"),
            verbatimTextOutput("mQTLT")
          )

        )

        ),

tabItem(tabName = "Interval_Mapping",
        titlePanel("Interval Mapping"),

        # Sidebar with a slider input for number of bins
        sidebarLayout(
          sidebarPanel(
            actionButton("runim", "Run", width='100%'),
            selectInput("plotim", "Function",
                        choices = c("Properties of Cross",
                                    "Plot",
                                    "LOD Thresholds",
                                    "LOD Intervals")),


            numericInput("lenim",
                         "Length of Chromosome",
                         value = 100),
            #num.families,var.total,var.genetic,var.env,maf,seed,row,column
            numericInput("nmarkim",
                         "Number of Markers",
                         value = 10),

            numericInput("indim",
                         "Number of Individuals",
                         value = 100),


            selectInput("spaceim",
                        "Equal Spacing",
                        choices = c("TRUE","FALSE")),

            selectInput("typeim",
                        "Type of Cross",
                        choices = c("bc","f2")),


            selectInput("qtlim",
                        "No. of QTLs",
                        choices = c("1","2")),

            selectInput("metim",
                        "Interval Mapping Method",
                        choices = c("EM","HK","EHK","All")),

            numericInput("permim",
                         "No. of Permutations",
                         value = 0),

            numericInput("pos1im",
                         "QTL 1- Position",
                         value = 10),

            numericInput("eff1im",
                         "QTL 1-Effect",
                         value = 5),

            numericInput("pos2im",
                         "QTL 2- Position",
                         value = 20),

            numericInput("eff2im",
                         "QTL 2-Effect",
                         value = -5),

            numericInput("seedim",
                         "Enter Seed",
                         NULL),

            checkboxInput("setseedim", "Set Seed", F)



          ),

          # Show a plot of the generated distribution
          mainPanel(
            plotOutput("IntMapP"),
            verbatimTextOutput("IntMapT")
          )

        )

        ),

tabItem(tabName = "Covariance_Estimation",
        #fluidPage(
        titlePanel("Simulating Genotype of Parent by Child- Estimating Covariance"),

        # Sidebar with a slider input for number of bins
        sidebarLayout(
          sidebarPanel(

            actionButton("runcov", "Run", width='100%'),

            numericInput("int",
                         "Enter n",
                         value = 10),

            numericInput("seedcov",
                        "Enter Seed",
                        NULL),

            checkboxInput("setseedcov", "Set Seed", F)
          ),

         #Show a plot of the generated distribution
         mainPanel(
          verbatimTextOutput("PCCov")
        )
        )
        )
)


        )
#)
#)



sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Simulate_Genotypes", tabName = "Simulate_Genotypes"),
    menuItem("Breeding_Values", tabName = "Breeding_Values"),
    menuItem("Interval_Mapping", tabName = "Interval_Mapping"),
    menuItem("Multiple_QTL_Mapping", tabName = "Multiple_QTL_Mapping"),
    menuItem("Covariance_Estimation", tabName = "Covariance_Estimation"),
    menuItem("Haseman-Elston_Regression", tabName = "Haseman-Elston_Regression"),
    id = "tab" #Candice
  )
)

ui <- dashboardPage(skin = "black",
              header,
              sidebar,
            footer = dashboardFooter(
                #left = "By Divad Nojnarg",
              left = tagList(
                "Created by ",
                tags$a(
                  href = "https://csammons7.wordpress.com",
                  target = "_blank",
                  "Candice Sammons"
                ),
                "and ",
                tags$a(
                  href = "https://pc-piyush.github.io",
                  target = "_blank",
                  "Piyush Chaudhari"
                ),
              ),
                right = "© 2018"
              ),
              body
)
