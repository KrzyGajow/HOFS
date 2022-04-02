runShinyHOFS <- function( launch.browser = T ){

  app <- shinyApp(

    #### UI ####
    ui = navbarPage( "HOFS",

      #### Description ####
      tabPanel( "Description",

        fluidPage(
          h2("Higher Order Mutual Information Approximation for Feature Selection"),

          br(),

          h3("Description"),
          "A higher order Mutual Information (MI) based approximation technique. Instead of producing a single list of features,
           method produces a ranked collection of feature subsets that maximizes MI,
           giving better comprehension (feature ranking) as to which features work best together when selected,
           due to their underlying interdependent structure.",
          br(),

          h3("Author"),

          a(href="http://krzysztof_gajowniczek.users.sggw.pl/", "Krzysztof Gajowniczek, PhD"),

          h3("Uploading Files"),

          p(strong("data:"), "Data.frame in which to interpret the parameters Target and Attributes.
          Accepted extensions: text/csv, text/comma-separated-values,text/plain, .csv, .arff"),

          h3("Parameters"),

          p(strong("data:"), "Table with potential features."),
          p(strong("label:"), "Target variable."),
          p(strong("C:"), "Predefined constant C, assessing whether create new subset."),
          p(strong("Nfeatures:"), "Number of relevant features to be selected."),
          p(strong("samplePct:"), "Percentage of observations that are sampled to speed up the computation."),
          p(strong("seed:"), "Seed for Pseudo-Random Number Generator."),

        )
      ),

      #### Run Selection ####
      tabPanel("Run Selection",

       shinyjs::useShinyjs(),

        tabsetPanel(

          tabPanel( "Uploading Files",

            sidebarLayout(

              sidebarPanel(

                fileInput("fileS", "Choose Table", multiple = TRUE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".arff")),

                hr(),

                checkboxInput("headerS", "Header", TRUE),

                radioButtons("sepS", "Separator", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = ","),

                radioButtons("quoteS", "Quote", choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"), selected = '"')

              ),

              mainPanel(

                tableOutput("InputFileS")

              )

            )

          ),

          tabPanel( "Parameters",

            sidebarLayout(

              sidebarPanel(

                actionButton( "startS", "Start Selection"),

                hr(),

                textOutput("wrongTargetS"),

                selectInput("Y_nameS", "Target", ""),

                selectInput("X_namesS", "Attributes", "", multiple = TRUE),

                numericInput("typeS", "Number of features", value = 10, 1, Inf, 1),

                numericInput("cS", "C", value = 0.2, 0, 1, 0.1 ),

                numericInput("pctS", "Percentage of observations", value = 0.1, 0.01, 1, 0.1 ),

                numericInput("seedS", "Seed fo PRNG", value = 666, 1, Inf, 1)

              ),

              mainPanel(

                tabPanel( "Tree",

                  # br(),

                  hr(),

                  verbatimTextOutput("plotS")


                )

              )

            )

          )

        )

      )

    ),

    #### Server ####
    server <- function( input, output, session ) {

      #### Init ####
      options( width = 10000 )

      sapply( c("OkTargetS","OkfileS","TarLevelsS"), function(x){ assign(x, T, envir = .GlobalEnv) } )

      #### Input files ####
      input_fileS <- reactive({

        req( input$fileS )

        OkfileS <<- T

        if( length( grep("arff",input$fileS$name) ) ){

          dfS <- read.arff( input$fileS$datapath )

        }else{

          dfS <- read.table( input$fileS$datapath, header = input$headerS, sep = input$sepS, quote = input$quoteS, check.names = F, stringsAsFactors = T )

        }

        if( any( sapply( dfS, anyNA ) ) ){

          OkfileS <<- F
          dfS <- data.frame( " " = "*** Missing values are not allowed ***", check.names = F )

        }

        return( dfS )

      })

      #### Dynamic names ####
      observeEvent( input_fileS(), ignoreInit = T, {

        input_namesS <- colnames( input_fileS() )

        updateSelectInput( session, "Y_nameS", choices = input_namesS, selected = "" )
        updateSelectInput( session, "X_namesS", choices = input_namesS )

      })

      #### Files render ####
      output$InputFileS <- renderTable( input_fileS() )

      #### Run Selection ####
      observeEvent( input$Y_nameS, ignoreInit = T, {

        req( OkfileS == T & input$Y_nameS != "" )

        if( is.numeric( input_fileS()[,input$Y_nameS] ) ){

          OkTargetS <<- F
          output$wrongTargetS <- renderText( "Target should be a factor" )

        }else{

          OkTargetS <<- T
          output$wrongTargetS <- renderText("")

        }

      })

      observeEvent( input$startS, {

        req( input$Y_nameS, input$X_namesS, OkTargetS == T,
             input$pctS > 0 & input$pctS <= 1, input$cS > 0 & input$cS <= 1 )

        ModelS <- HOFS:::HOFS_Shiny( data = input_fileS()[, input$X_namesS ],
                                     label = input_fileS()[, input$Y_nameS, drop = F ],
                                     C = input$cS,
                                     samplePct = input$pctS,
                                     Nfeatures = input$typeS,
                                     verbose = F,
                                     seed = input$seedS )

        output$plotS <- renderPrint({

          ModelS

        })

      })

    }

  )

  #### RunApp ####
  runApp( app, launch.browser = launch.browser )

}
