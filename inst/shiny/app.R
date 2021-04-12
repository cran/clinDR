#######################################################
#' clinDR shiny application
#'
#' Scope of use: 
#'  Planning of dose-response studies. The application 
#'  provides the user with basic settings only for clinDR
#'  emaxsim and emaxsimB simulations. However, code is 
#'  presented to allow the user to recreate results outside
#'  of the application. This code can be "tweaked" to use
#'  any more advanced settings required.
#'
#' The application provides a wrapper aI round the clinDR
#' package written by Neal Thomas (Pfizer) and published 
#' on CRAN v1.8 (https://cran.r-project.org/web/packages/clinDR/index.html)
#'
#' Author / Code maintainer: Mike K Smith
#' Date: 19 May 2020
#' Version: 1.4
#'
#' Revisions:
#' v1.0 03 December 2018 - Initial version (MKS)
#' v1.1 02 January 2019 - Add markdown reporting download functionality
#' v1.2 20 Feburary 2020 - Update to clinDR v2.2.1
#' v1.3 22 April 2020 - Refactor
#' v1.4 19 May 2020 - Updates following review by Neal Thomas:
#'                    Add sample variability to inputs visualisation.
#'                    Add ability to specify mean levels directly (rather than
#'                    via parameters)
#' v1.5 02 September 2020 - Incorporate input from reviewers re: UI and from Neal 
#'                     Thomas re: Bayesian priors. Test with clinDR 2.3
#' v1.6 17 November 2020 - Handle user feedback if prior settings not reviewed
#' v1.7 22 March 2021 - Minor changes to handling variable selection in
#'                    fit.quantiles and associated glue code for reproducibility.
#' v1.8 31 March 2021 - Revising prior to CRAN submission
#'
#' Pre-requisites:
#'  This code requires the following packages in addition to shiny
#'    * clinDR
#'      - Need to run the clinDR function compileStanModels() before execution.
#'    * glue
#'    * dplyr, tidyr, purrr, tibble, magrittr
#'    * waiter
#'
#########################################################

# tools:::Rd2HTML(utils:::.getHelpFile(help(package = clinDR,
#                                           topic = emaxsim)),
#                 out = "emaxsim.html")
# 
# tools:::Rd2HTML(utils:::.getHelpFile(help(package = clinDR,
#                                           topic = emaxsimB)),
#                 out = "emaxsimB.html")
# tools:::Rd2HTML(utils:::.getHelpFile(help(package = clinDR,
#                                           topic = emaxPrior.control)),
#                 out = "emaxPrior_control.html")

requireNamespace("clinDR", quietly=TRUE)
requireNamespace("ggplot2", quietly=TRUE)
requireNamespace("shiny", quietly=TRUE)

waiting_screen <- tagList(
  waiter::spin_fading_circles(),
  h4("Simulations running...")
)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel(title = "Dose-Response simulation using clinDR"),
  withTags({
    div(class="header", checked=NA,
        p("Simulation of dose response studies using Bayesian or Maximum 
        Likelihood Emax model estimation.  Binary and continuous data 
        simulations are supported.  For binary data, the Emax model is logit 
        transformed.  This application executes the emaxsim or emaxsimB 
        functions in R package clinDR"),
        p(strong("Shiny App maintainer:"), 
          a(href="mailto:mike.k.smith@pfizer.com","Mike K Smith")),
        p(strong("clinDR package developer:"), 
          a(href="mailto:neal.thomas@pfizer.com", "Neal Thomas")),
        br()
    )
  }),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(type="tabs",
                  tabPanel("Simulation Settings", 
                           style = "overflow-y:scroll; max-height: 900px; position:relative;", 
                           h3("Simulation inputs"),
                           h4("Design settings"),
                           textInput("doselev", 
                                     label="Dose levels - separate with commas",
                                     value="0,5,25,50,100"),
                           textInput("n",
                                     label="Number subjects on each dose - separate with 
                           commas",
                                     value="100, 50, 50, 50, 100"),
                           hr(),
                           
                           h4("Simulation model"),
                           checkboxInput("binary", 
                                         label = "Binary outcome", 
                                         value = FALSE),
                           
                           selectInput("inputVals", 
                                       label = "Specify dose response curve", 
                                       choices = list("Using means/proportions for each dose" = 0,
                                                      "Using Emax model parameters" = 1),
                                       selected = 1),
                           uiOutput("inputParameters"),
                           uiOutput("inputMeanlev"),
                           conditionalPanel(condition = "input.binary == 0",
                                            verticalLayout(
                                              numericInput("resSD", 
                                                           label="Residual SD", 
                                                           value=8, step = 0.1)
                                            )
                           ),
                           
                           hr()),
                  tabPanel("Analysis Settings",
                           h3("Analysis settings"),
                           selectInput("modType", 
                                       label="Emax model fitted",
                                       choices=list("Hyperbolic (3-parameter) Emax"=3,
                                                    "Sigmoidal (4-parameter) Emax"=4),
                                       selected = 4),
                           checkboxInput("bayes", 
                                         label = "Use Bayesian methods", 
                                         value = TRUE),
                           uiOutput("inputBayesPriors")
                  )
      )
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  ### Displaying what the Emax dose-response looks like for 
                  ### a given set of inputs. Helps the user refine input values.
                  tabPanel("Inputs check",
                           h3("Check your inputs provided in the 'Simulation 
                           Settings' tab on the left using the graph and table below"),
                           strong("Once you are happy with the Simulation Settings
                           click on the 'Run & Summary' tab above to run 
                                  simulations and view outputs."),
                           br(),
                           br(),
                           p("Visual check of the specified inputs.
          Plot of response versus dose. If using binary responses / proportions 
          please ensure the 'Binary outcome' box is checked to see accurate 
          description of the Emax dose-response for this type of data.
                             The plot also optionally shows sampling variability
                             (approximate 95% CIs) of the mean values for the 
                             given inputs."),
                           
                           checkboxInput("showSampMean", 
                                         label = "Show sampling variability", 
                                         value = FALSE),
                           br(),
                           plotOutput("inputShow"),
                           br(),
                           p("For the given inputs, the table below shows the 
                             predicted outcome at each dose."),
                           tableOutput("predictionTable")),
                  tabPanel("Run & Summary",
                           h4("Check Analysis Settings in the left-hand tab
                           then run simulations by clicking the run button below"),
                           waiter::use_waiter(),
                           fluidPage(
                             h3("Simulation settings"),
                             fluidRow(
                               column(width=4,
                                      br(),
                                      actionButton("Run", "Run" )),
                               column(width=4,
                                      numericInput("seed", 
                                                   label="Simulation seed", 
                                                   value= floor(
                                                     runif(1,0, 99999))
                                                   )),
                               column(width=4,##NT change 
                                      uiOutput("nsimsIO")),
                               column(width=4, 
                                      uiOutput("nprocIO"))
                             ),
                             p("Check with the system administrator before 
                             requesting more than 1 processor on a multi-user 
                             system.  The number of simulations should be 
                             divisible by the number of processors for 
                            computational efficiency"),
                             hr(),
                             ### Default output from clinDR. Refer to {clinDR}
                             ### package for more information.
                             fluidRow(column(
                               width=12,
                               p("For pairwise comparisons, the 
                               'most favorable pairwise comparison' means the 
                               dose with the best difference versus placebo is 
                               compared to the population mean response for the 
                               selected dose, thus the target value for 
                               coverage, bias, and RMSE changes depending on the 
                               selected dose."),
                               verbatimTextOutput("result")))
                           )),
                  tabPanel("Visual Summary",
                           ### Visual summary of simulation results. 
                           h2("Histograms and correlations between parameter estimates."),
                           p("Based on estimates or posterior medians. One value 
                           per simulation. For maximum likelihood methods i.e. 
                           not Bayesian analysis this shows only estimates for 
                           the nominated Emax model i.e. either 3 or 4 parameter
                             Emax."),
                           p(strong("This plot takes a moment to render. Please 
                                    be patient.")),
                           br(),
                           plotOutput("histplot", width = "100%"),
                           br(),
                           tableOutput("fitType"),
                           br(),
                           h2("Q-Q plot of comparisons"),
                           uiOutput("QQPlotText"),
                           plotOutput("QQplot", width = "95%")
                  ),
                  tabPanel("Individual Simulation Results",
                           ### plot and table of results from individual 
                           ### simulation replicates.
                           uiOutput("indivResults")
                  ),
                  tabPanel("Code",
                           ### code for reproducibility. Users should copy and
                           ### paste this into an R session to recreate results.
                           p("The code below can be used to reproduce the 
                             results provided through this application. 
                             It can also be used for further refinement of 
                             simulation settings and conditions."),
                           p("If using Bayesian analysis methods, remember to
                             click on the 'Analysis Settings' tab on the left
                             to specify priors for analysis. If you do not 
                             do so, the prior definition code will not be shown
                             in the code chunk below."),
                           br(),
                           # p("Alternatively, click on the 'Download report' 
                           # button below to export an HTML of the output"),
                           # downloadButton("report", "Download report"),
                           verbatimTextOutput("code")
                  ),
                  tabPanel("HELP",
                           h3("Help for key {clinDR} package functions"),
                           p("Below you'll find the manual pages for key
                             {clinDR} functions `emaxsim`, `emaxsimB` and
                             `prior.control`."),
                           fluidPage( 
                             style = "overflow-y:scroll; max-height: 800px; position:relative;",
                             tabsetPanel(
                               tabPanel("emaxsim", includeHTML("emaxsim.html")),
                               tabPanel("emaxsimB", includeHTML("emaxsimB.html")),
                               tabPanel("emaxPrior.control", includeHTML("emaxPrior_control.html"))
                             ))
                  ),
                  tabPanel("About",
                           includeMarkdown("clindr_Shiny_application_documentation.md")
                  )
                  
      )
    )
  )
)

server <- function(input, output) {
  
  w <- waiter::Waiter$new()
  
  ShinyPlotTheme <- theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))
  
  ######################################################################
  ##
  ## DESIGN INPUTS
  ##
  ######################################################################
  
  doselev <- reactive({
    as.numeric(unlist(strsplit(input$doselev,",")))
  })
  
  n <- reactive({
    as.numeric(unlist(strsplit(input$n,",")))
  })
  
  Ndose<- reactive({
    length(doselev())
  })
  

  doseSeq <- reactive({seq(min(doselev()), max(doselev()),length=100)})
  
  ######################################################################
  ##
  ## Simulation inputs
  ##
  ######################################################################
  
  ## POPULATION ESTIMATES FOR GENERATING MEAN VALUES PER DOSE
  ##  User specifies E0, target effect at specified dose
  ##   and clinDR function solveEmax calculates the Emax for simulation.
  
  ## For binary outcomes: User specifies E0, target effect on logit scale
  ##   pop() will use solveEmax to calculate Emax.
  
  ## targetDose defaults to maximum value of doses provided in design
  ## but the user can select other doses. These are then used to derive Emax
  ## using the clinDR function solveEmax.
  
  targetDose <- reactive({
    return(max(doselev()))
  })
  
  output$inputParameters <- renderUI({
    if(input$binary){
      conditionalPanel(condition = "input.inputVals == 1",
                       verticalLayout(
                         h4("Population parameters"),
                         numericInput("e0", 
                                      label="E0 (logit scale)", 
                                      value=0, step = 0.1),
                         numericInput("ed50", 
                                      label="ED50", 
                                      value=15, step = 1),
                         numericInput("target", 
                                      label="Target effect at specified dose 
                                      (logit scale)", 
                                      value=1.5, step = 0.1),
                         numericInput("targetDose", 
                                      label="Dose achieving specified target effect", 
                                      value = targetDose()),
                         checkboxInput("pboadj", 
                                       label="Is target effect placebo adjusted?", 
                                       value = TRUE),
                         numericInput("lambda", 
                                      label="Lambda", 
                                      value=1, step = 0.1)
                       )
      )
    } else {
      conditionalPanel(condition = "input.inputVals == 1",
                       verticalLayout(
                         h4("Population parameters"),
                         numericInput("e0", 
                                      label="E0", 
                                      value=2, step = 0.1),
                         numericInput("ed50", 
                                      label="ED50", 
                                      value=15, step = 1),
                         numericInput("target", 
                                      label="Target effect at specified dose", 
                                      value=6, step = 0.1),
                         numericInput("targetDose", 
                                      label="Dose achieving specified target effect", 
                                      value = targetDose()),
                         checkboxInput("pboadj", 
                                       label="Is target effect placebo adjusted?", 
                                       value = TRUE),
                         numericInput("lambda", 
                                      label="Lambda", 
                                      value=1, step = 0.1)
                       )
      )}
  })
  
  output$inputMeanlev <- renderUI({
    ### Defaults for meanlev correspond to values for Emax model given 
    ### default input parameter values above.
    if(input$binary){
      conditionalPanel(condition = "input.inputVals == 0",
                       verticalLayout(
                         textInput("meanlev", 
                                   label="Population response rate (proportions) 
                                   at each dose",
                                   value="0.5,0.61,0.75,0.79,0.82")
                       )
      )
    } else {
      conditionalPanel(condition = "input.inputVals == 0",
                       verticalLayout(
                         textInput("meanlev", 
                                   label="Population mean at each dose", 
                                   value="2.0, 3.72, 6.31, 7.31, 8.0")
                       )
      )
    }
  })
  
  pop <- reactive({
    if(input$inputVals == 1){
      validate(
        ### Checks that input values are valid.
        need(input$ed50 & input$ed50 > 0, "ED50 must be a positive number"),
        need(input$e0, "E0 must be supplied"),
        need(input$target, "Target effect at specified dose must be supplied")
      )
      led50 <- log(input$ed50)
      lambda <- input$lambda
      E0 <- input$e0
      Emax <- clinDR::solveEmax(target = input$target, 
                        dose = input$targetDose,
                        led50 = led50,
                        lambda = lambda,
                        e0 = E0,
                        pboadj = input$pboadj)
      c(led50=led50, lambda=lambda, emax=Emax, e0=E0)
    } else {
      return(NULL)
    }
  })
  
  ## For given inputs, calculate population mean effect at each dose
  
  meanlev <- reactive({
    validate(
      need(length(doselev())>2, "Number of doses must be greater than 2 for emax
           model fit"),
      need(length(doselev()) == length(n()), 
           "Must supply a sample size for each dose arm. Too many doses or too 
           many n values supplied.")
    )
    if(input$inputVals==1){
      if(!input$binary){
        emaxfun(dose = doselev(), parm = pop())
      } else {
        plogis(emaxfun(dose = doselev(), parm = pop()))
      }
    } else {
      values <- as.numeric(unlist(strsplit(input$meanlev,",")))
      validate(
        need(length(doselev()) == length(values),
             "Need a mean response level for each dose. Too many doses or too 
             many mean responses supplied.")
      )
      if(input$binary){
        validate(
          need((min(values)>0 & max(values)<1), "All response inputs must be 
               between 0, 1.")
        )
      }
      values
    }
  })
  
  ## collect inputs for simulation
  
  gen <- reactive({
    if(!input$binary){
      clinDR::FixedMean(n = n(), 
                doselev = doselev(), 
                meanlev = meanlev(),
                parm = pop(),
                resSD = input$resSD) 
    } else {
      clinDR::FixedMean(n = n(), 
                doselev = doselev(), 
                meanlev = meanlev(),
                parm = pop(),
                binary = TRUE)
    }
  })
  
  ## For the inputs checking tab, we calculate model predictions across a 
  ## wider dose-range 0 - maximum dose specified in design.
  
  predictions <- reactive({
    if(input$inputVals==1){
      if(!input$binary){
        clinDR::emaxfun(dose = doseSeq(), parm = pop())
      } else {
        plogis(clinDR::emaxfun(dose = doseSeq(), parm = pop()))
      }
    } else {
      meanlev()
    }
  })
  
  ## Calculate lower and upper CI of sample mean variability.
  ## For binary, this is based on quantiles of observed proportions from a 
  ## binomial draw using n per dose from simulation design.
  
  loCImean <- reactive({
    if(!input$binary){
      meanlev() - 1.96*(input$resSD/sqrt(n()))
    } else {
      sampleBinom <- purrr::map2(n(), meanlev(),
                          .f = function(.x, .y){rbinom(n=10000, .x, .y)/.x}) %>%
        purrr::map_dbl(quantile,probs=0.025)
    }
  })
  
  hiCImean <- reactive({
    if(!input$binary){
      meanlev() + 1.96*(input$resSD/sqrt(n()))
    } else {
      sampleBinom <- purrr::map2(n(), meanlev(),
                          .f = function(.x, .y){rbinom(n=10000, .x, .y)/.x}) %>%
        purrr::map_dbl(quantile,probs=0.975)
    }
  })
  
  ## Plot the dose-response and associated sampling variability
  
  output$inputShow <- renderPlot({
    validate(
      need(length(doselev())>2, "Number of doses must be greater than 2 for emax
           model fit"),
      need(length(doselev()) == length(n()), 
           "Must supply a sample size for each dose arm. Too many doses or too 
           many n values supplied.")
    )
    if(input$inputVals==1){
      plot1 <- ggplot() +
        geom_line(mapping = aes(x = doseSeq(), y = predictions())) + 
        labs(x = "Dose", y = "Response")  +
        ShinyPlotTheme
    } else {
      plot1 <- ggplot() +
        geom_line(mapping = aes(x = doselev(), y = meanlev())) + 
        labs(x = "Dose", y = "Response") +
        ShinyPlotTheme
    }
    
    if(!input$showSampMean){
      print(plot1)
    } else {
      plot1 +
        geom_ribbon(mapping = aes(x = doselev(),
                                  ymin = loCImean(),
                                  ymax = hiCImean()),
                    fill = "blue", alpha=0.2)
    }
  })
  
  ## Create a table of population means for each of the doses in the input design.
  
  output$predictionTable <- renderTable({
    data.frame(Dose= doselev(),
               Response = meanlev())
  })
  
  ## Collect together inputs for simulation.
  
  parameters <- reactive({
    doselev <- doselev()
    n <- n()
    Ndose<-length(doselev())
    ### population parameters for simulation
    pop <- pop() 
    meanlev <- meanlev()
    gen <- gen()      
  })
  
  ######################################################################
  ##
  ## Bayesian analysis prior specification
  ##
  ######################################################################  
  
  ## Because some prior values are calculated based on other inputs, we
  ## need to perform the calculation in the server() function and then
  ## pass the UI back into the ui() function.
  
  ## Neal Thomas has provided input to reasonable settings of Bayesian priors.
  ## The user is free to select their own values.
  
  prior.p50 <- reactive({
    doses <- doselev()
    nzDoses <- doses[doses!=0]
    round(nzDoses[1] + (nzDoses[2] - nzDoses[1])/2, 1)
  })
  
  prior.difTargetmu <- reactive({
    tDose <- doselev() == input$dTarget
    difTargetmu <- ifelse(input$binary,
                          round(qlogis(meanlev()[tDose],4)),
                          meanlev()[tDose]
    )
    return(difTargetmu)
  })
  
  ## Neal Thomas suggests prior for scale parameters should be 5*resSD
  prior.difTargetsca <- reactive({
    ifelse(input$binary, 4, round(10*input$resSD,1))
  })
  
  prior.epsca <- reactive({
    ifelse(input$binary, 4, round(10*input$resSD,1))
  })
  
  prior.epmu <- reactive({
    ifelse(input$binary,round(qlogis(meanlev()[1]),4), meanlev()[1])
  })
  
  prior.dTarget <- reactive({
    if(input$inputVals==1){
      input$targetDose
    } else {
      targetDose()
    }
  })
  
  output$inputBayesPriors <- renderUI({
    conditionalPanel(condition = "input.bayes == 1",
                     if(!input$binary){
                       verticalLayout(
                         h4("Prior settings"),
                         numericInput("epmu",
                                      label="Prior Mean E0", 
                                      value=prior.epmu()),
                         ## Need to specify UI components in server function
                         ## as these rely on inputs for default value calculation
                         numericInput("epsca", 
                                      label="Prior Scale parameter for E0", 
                                      value=prior.epsca()),
                         numericInput("dTarget", 
                                      label="Target dose for priors", 
                                      value= prior.dTarget()),
                         numericInput("difTargetmu", 
                                      label="Prior Mean for effect at dTarget dose", 
                                      value=prior.difTargetmu()),
                         numericInput("difTargetsca", 
                                      label="Scale parameter for effect at dTarget dose", 
                                      value=prior.difTargetsca()),
                         numericInput("p50",
                                      label="Projected ED50",
                                      value = prior.p50()),
                         numericInput("sigmalow", 
                                      label="Lower bound for residual SD", 
                                      value=input$resSD/10),
                         numericInput("sigmaup", 
                                      label="Upper bound for residual SD", 
                                      value=input$resSD*10)
                       )
                     } else {
                       verticalLayout(
                         h4("Prior settings"),
                         numericInput("epmu",
                                      label="Prior Mean for E0 (logit)", 
                                      value=prior.epmu()),
                         ## Need to specify UI components in server function
                         ## as these rely on inputs for default value calculation
                         numericInput("epsca", 
                                      label="Prior Scale parameter for E0", 
                                      value=prior.epsca()),
                         numericInput("dTarget", 
                                      label="Target dose for priors", 
                                      value= input$targetDose),
                         numericInput("difTargetmu", 
                                      label="Prior Mean for effect at target dose", 
                                      value=input$target),
                         numericInput("difTargetsca", 
                                      label="Scale parameter for effect at target dose", 
                                      value=prior.difTargetsca()),
                         numericInput("p50",label="Projected ED50",
                                      value = prior.p50())
                       )})
  })
  
  
  prior <- reactive({
    if(!input$binary){
      prior <- clinDR::emaxPrior.control(epmu = input$epmu,
                                 epsca = input$epsca,
                                 difTargetmu = input$difTargetmu,
                                 difTargetsca = input$difTargetsca,
                                 dTarget = input$dTarget,
                                 p50 = input$p50,
                                 sigmalow = input$sigmalow,
                                 sigmaup = input$sigmaup,
                                 parmDF = 5)
    } else {
      prior <- clinDR::emaxPrior.control(epmu = input$epmu,
                                 epsca = input$epsca,
                                 difTargetmu = input$difTargetmu,
                                 difTargetsca = input$difTargetsca,
                                 dTarget = input$dTarget,
                                 p50 = input$p50,
                                 parmDF = 5,
                                 binary=TRUE)
    }
    return(prior)
  })
  
  ######################################################################
  ##
  ## Run clinDR emaxsim or emaxsimB
  ##
  ## NOTE: Because of reactivity rules, the simulation DOES NOT run
  ## until the user clicks through to the output tabs (reactive values
  ## are not calculated until needed)
  ##
  ######################################################################  
  
  ### add number of processors for parallel computation
    maxprocs <- parallel::detectCores()
    np <- ifelse(.Platform$OS.type == 'windows', maxprocs, 1)
    

  output$nprocIO <- renderUI({
    maxprocs <- parallel::detectCores()
    labtxt <- paste0("Number of processors [1-", maxprocs, "]")
    sliderInput(
      "nprocs",
      label = labtxt,
      value = np,
      min = 1,
      max = maxprocs,
      step = ifelse(np<=8, 1, 2)
    )
  })

  nsims <- reactive({
     req(input$nprocs)
     if (input$nprocs<3) nsims <- 20
     else if (input$nprocs>2 & input$nprocs<5) nsims <- 24
     else nsims <- 5 * input$nprocs
    return(nsims)
  })
  
  output$nsimsIO <- renderUI({
    numericInput(
      "nsims",
      label = "Number of simulations",
      value = nsims(),
      min = 1,
      max = 1000
    )
  })
  
  simulation <- eventReactive(input$Run,{
    req(gen())
    set.seed(input$seed)
    if(!input$bayes){
      waiter::waiter_show(html = waiting_screen, color = "grey")
      
      if(!input$binary){
        results <- clinDR::emaxsim(nsim = input$nsims,
                           genObj = gen(),
                           modType=as.numeric(input$modType),
                           nproc = input$nprocs)
      } else {
        results <- clinDR::emaxsim(nsim = input$nsims,
                           genObj = gen(),
                           modType=as.numeric(input$modType),
                           nproc = input$nprocs,
                           binary=TRUE)                       
      }
    } else {
      req(prior())
      
      rstan::rstan_options(auto_write = TRUE)
      # options(mc.cores = parallel::detectCores())
      # Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
      
      mcmc<-clinDR::mcmc.control(chains=1,warmup=500,iter=5000,
                         propInit=0.15,adapt_delta = 0.95)
      

        showNotification("You MUST review Bayesian priors in the 'Analysis 
        Settings' tab BEFORE running the simulation with estimation
        using Bayesian analysis method",
                         type = "error",
                         duration = NULL)
        
        req(prior(), cancelOutput = TRUE)

        waiter::waiter_show(html = waiting_screen, color = "grey")
      
      results <- clinDR::emaxsimB(nsim = input$nsims,
                          genObj = gen(),
                          prior = prior(),
                          modType=as.numeric(input$modType),
                          mcmc=mcmc,
                          check=FALSE,
                          nproc = input$nprocs,
                          binary=(input$binary))
    }
    waiter::waiter_hide()
    return(results)
  })
  
  ######################################################################
  ##
  ## Simulation summaries aggregated across simulated trials
  ##
  ######################################################################  
  
  ## estimates holds the parameter estimates from the model fits
  
  estimates <- reactive({
    req(simulation())
    
    if(simulation()$modType==3 & class(simulation()) =="emaxsim"){
      est <- simulation()$est3
    }
    if(simulation()$modType==4 & class(simulation()) =="emaxsim"){
      est <- simulation()$est4
    }
    if(class(simulation())=="emaxsimB"){est <- simulation()$est}
    
    est <- est %>%
      tibble::as_tibble() %>%
      dplyr::mutate(ED50 = exp(led50)) %>%
      dplyr::select(-led50) %>%
      dplyr::rename(E0 = e0,
             Emax = emax)
    
    if(simulation()$binary){
      est 
    } else {
      est %>%
        dplyr::mutate('Residual_SD' = simulation()$residSD)
    }
  })
  
  ## Standard clinDR print.emaxsim/emaxsimB output
  
  output$result <- renderPrint({
    req(simulation())
    summary(simulation())
  })
  
  ## Standard clinDR plot.emaxsim/emaxsimB output
  
  output$QQPlotText <- renderUI({
    if(input$bayes == FALSE){
      tagList(
        p("A Q-Q plot of the dose response estimate of the 
                             mean at the highest dose minus its population value 
                             divided by the standard error of the estimator 
                             (computed using the delta method)"),
        p("Results based on alternative model fits 
                           i.e. where fitted model is not as specified in 
                           the Analysis Settings are displayed in red")
      )
    } else {
      tagList(
        p("A Q-Q plot of the posterior median at the highest dose minus its 
        population value divided by the posterior standard deviation")
      )
    }
  })
  
  output$QQplot <- renderPlot({
    # req(simulation())
    if(input$bayes){ clinDR::plot.emaxsimB(simulation())
    }else clinDR::plot.emaxsim(simulation())
  })
  
  ## Table of fit types for ML. Does not apply to Bayesian analysis.
  
  output$fitType <- renderTable({
    ## Only calculate for non-Bayesian estimation.
    ## Bayesian approach does not use alternative model estimation
    if(class(simulation()) =="emaxsimB") return()
    if(class(simulation()) =="emaxsim"){
      nsims <- nrow(simulation()$results$predpop)
      
      fitType <- factor(simulation()$fitType,
                        levels = c("4","3","L","LL","E"))
      dplyr::bind_rows(summary(fitType)) %>%
        tidyr::pivot_longer(cols = tidyselect::everything(),
                     names_to = "Var1",
                     values_to = "number") %>%
        dplyr::mutate(Var1 = dplyr::recode(Var1,
                             "3" = "3-parameter Emax",
                             "4" = "4-parameter sigmoid Emax",
                             "L" = "Linear",
                             "LL" = "Log Linear",
                             "E" = "Exponential")) %>%
        dplyr::rename("Fit Type" = Var1)
    }
  })
  
  ## Scatter plot matrix of parameter estimates. 
  ## For ML this is only applies for the chosen type of model 
  ##   i.e. if 4-parameter Emax model is chosen, then only displays estimates
  ##   from 4-parameter Emax models.
  ## For Bayesian analysis / MCMC plots the posterior median parameter estimates 
  ##   for each simulation
  
  output$histplot <- renderPlot({
    req(estimates())
    graphics::pairs(estimates())
  })
  

  ######################################################################
  ##
  ## Individual simulated trial results
  ##
  ###################################################################### 
  
  ## Plot individual trial results with associated interval estimates 
  ##   of predictions at each dose.
  
  plot.indivResult <- function(x, plotDif){
    if(!x$binary) {
      ylower <- min(x$fitpredv-2*x$sepredv)
      yupper <- max(x$fitpredv+2*x$sepredv)
    } else {
      ylower <- 0
      yupper <- 1
    }
    if(input$bayes){pout<-clinDR::plot.emaxsimBobj(x[input$sim],
         ylim = c(ylower, yupper),
         plotDif = plotDif)
    } else pout<-clinDR::plot.emaxsimobj(x[input$sim],
          ylim = c(ylower, yupper),
          plotDif = plotDif)
   return(pout)
  }
  output$plotIndiv <- renderPlot({
    req(simulation())
    plot.indivResult(simulation(), 
                     input$plotDif)
  })
  
  ## Prepare table of model parameter estimates for individual trial
  ## For binary, these are returned on the logit scale so need to be 
  ##   back-transformed to the proportion scale (0,1)
  ## 
  ## For Bayesian analysis, we're using median predictions for parameters
  
  print.indivPars <- function(x){
    fitdifv <- x$fitpred - x$fitpred[1]
    names(fitdifv) <- NULL
    fitdifP <- x$predpop - x$predpop[1]
    names(fitdifP) <- NULL
    negC <- x$negC
    bigC <- x$bigC
    pval <- round(x$pVal, 3)
    resid <- NA
    if(!input$bayes){
      est3 <- x$est3
      est4 <- x$est4
      if (x$fitType == "4") {
        est4[1] <- exp(x$est4[1])
        est4 <- round(est4,3)
        names(est4) <- NULL
      }
      if (x$fitType == "3") {
        est3[1] <- exp(x$est3[1])
        est3 <- round(est3,3)
        names(est3) <- NULL
      }
      if(!input$binary)resid <- round(x$residSD,3)
      noFit <- (x$fitType == "4" & any(is.na(est4)) | 
                  (x$fitType == "3" & any(is.na(est3))))
      
      
      zval <- round((fitdifv[x$idmax] - fitdifP[x$idmax])/x$sedif[x$idmax], 3)
      names(zval) <- NULL
      
      if (x$fitType== "4") {
        indivPars <- data.frame(fitType = x$fitType, 
                                noFit = noFit, 
                                negC = negC, 
                                bigC = bigC, 
                                ED50 = est4[1],
                                lambda = est4[2],
                                Emax = est4[3],
                                E0 = est4[4],
                                resid = resid,
                                zval = zval, 
                                pval = pval)
      } else {
        indivPars <- data.frame(fitType = x$fitType, 
                                noFit = noFit, 
                                negC = negC, 
                                bigC = bigC, 
                                ED50 = est3[1],
                                Emax = est3[2],
                                E0 = est3[3],
                                resid = resid,
                                zval = zval, 
                                pval = pval)
      }
    }
    
    if(input$bayes){
      if (x$modType == "4") {
        est4 <- summary(x$bfit$estanfit)$summary[c("led50","lambda","emax","e0[1]"),c(4,6,8)]
        est4[1,] <- exp(est4[1,])
       
        est4 <- paste(round(est4[,2],3),
                      " (",round(est4[,1],3),", ",
                      round(est4[,3],3),")",sep="")
      } else {
        est3 <- summary(x$bfit$estanfit)$summary[c("led50","emax","e0[1]"),c(4,6,8)]
        est3[1,] <- exp(est3[1,])
        
        est3 <- paste(round(est3[,2],3),
                      " (",round(est3[,1],3),", ",
                      round(est3[,3],3),")",sep="")   
      }
      resid <- ifelse(input$binary, 
                      NA, 
                      round(summary(x$bfit$estanfit)$summary["sigma[1]",6],3))
      
      if (x$modType== "4") {
        indivPars <- data.frame(fitType = x$modType, 
                                ED50 = est4[1],
                                lambda = est4[2],
                                Emax = est4[3],
                                E0 = est4[4],
                                resid = resid)
      } else {
        indivPars <- data.frame(fitType = x$modType, 
                                ED50 = est3[1],
                                Emax = est3[2],
                                E0 = est3[3],
                                resid = resid)
      }
    }
    
    return(indivPars)
  }
  
  output$printIndivPars <- renderTable({
    req(simulation())
    print.indivPars(simulation()[input$sim])
  })
  
  ## Calculate table of model predictions, differences and uncertainty
  ### estimates at each dose.
  
  indivEst <- function(x, sim ){
    results <- tibble::tibble(mean = x$fitpred[sim,],
                      se = x$sepred[sim,],
                      diff = x$fitpred[sim,] - x$fitpred[sim,1],
                      sediff = x$sedif[sim,])
    if(!input$bayes){
      results <- results %>%
        dplyr::mutate('0.25' = diff - qnorm(0.975)*sediff,
               '0.75' = diff + qnorm(0.975)*sediff)
    } else {
      results <- results %>%
        dplyr::mutate(diff = x$fitdif[sim,])
      lb <- x$lb[,sim,]
      lb <- rbind('0'=rep(0,ncol(lb)), lb)
      ub <- x$ub[,sim,]
      ub <- rbind('0'=rep(0,ncol(ub)), ub)
      ub <- ub[,c(3,2,1)]
      results <- cbind(results, lb, ub)
    }
    # results <- tibble::rownames_to_column(results, var="dose")
    results <- results %>%
      dplyr::mutate(dose = doselev()) %>%
      dplyr::select(dose, tidyselect::everything())
    return(results)
  }
  
  output$printIndivEst <- renderTable({
    req(simulation())
    indivEst(simulation(), input$sim)
  })
  
  ## UI output for individual results tables
  
  output$indivResults <- renderUI({
    MLEtableText <- "Table below presents MLEs"
    BayestableText <- "Table below presents posterior medians, lower 2.5%, upper 97.5%"
    tagList(
      h3("A plot of simulated data and associated model fit."),
      p(strong("This plot takes a moment to render. Please be patient.")),
      p(em("Please select a trial replicate number and press ENTER")),
      numericInput("sim",
                   label="Simulated dataset", 
                   value=1, min = 1, max = input$nsims),
      checkboxInput("plotDif", 
                    label="Plot changes from placebo?", 
                    value = FALSE),
      # validate(
      #   need(input$sim > 0 & input$sim <= input$nsims), 
      #   "Individual results requested must be within the number of simulated 
      #   replicates"
      # ),
      
      plotOutput("plotIndiv"),
      p("Note:  Dashed curve is population, solid curve is estimated"),
      p("Red stars are observed means at each dose"),
      p("Black intervals are dose-group mean Interval estimates"),
      p("Grey intervals are predictive intervals"),
      p(glue::glue("90% interval estimates shown")),
      hr(),
      h3("Treatment estimates for each dose within the nominated study"),
      p("Table below presents fitted means, se's, differences from placebo and 
        associated uncertainty in these differences"),
      tableOutput("printIndivEst"),
      h3("Parameter estimates for the simulated study"),
      p(glue::glue({ifelse(input$bayes,BayestableText, MLEtableText)})),
      tableOutput("printIndivPars"),
      if(!input$bayes){
        tagList(
          p("noFit  : No convergence"),
          p("negC   : Converged with ED50<lower limit"),
          p("bigC   : Converged with ED50>upper limit"),
          p("StdBias: (estimate-population)/SE for highest dose vs PBO"),
          p("P-val  : MCP-Mod P-value for test of no drug effect for highest 
            dose vs placebo")
        )
      }
    )
  })
  
  
  ######################################################################
  ##
  ## Code for reproducibility
  ##  Aim is for colleagues to copy and paste this code into a new
  ##  R session in order to replicate what is seen in the Shiny application.
  ##  At present this uses `glue` functionality to paste inputs to appropriate
  ##  code. 
  ## 
  ## TODO:
  ##   Refactor to use the `shinymeta` package if possible.
  ##
  ###################################################################### 
  
  output$code <- renderText({
    design <- glue::glue("library(clinDR)",
                   "library(tidyverse)",
                   "",
                   "set.seed({input$seed})",
                   "",
                   "# Design aspects",
                   "nsim <- {input$nsims}",
                   "doselev <- c({input$doselev})",
                   "n <- c({input$n})",
                   "Ndose <- length(doselev)",
                   .sep="\n")
    
    param <- glue::glue("# Parameters",
                  "led50 <- log({input$ed50})",
                  "lambda <- {input$lambda}",
                  "e0 <- {input$e0}",
                  "target <- {input$target}",
                  "targetDose <- {input$targetDose}",
                  "emax <- solveEmax(target = target, 
                        dose = targetDose,
                        led50 = led50,
                        lambda = lambda,
                        e0 = e0,
                        pboadj = {input$pboadj})",
                  "",
                  "pop <- c(led50 = led50, lambda = lambda, emax = emax, e0 = e0)",
                  .sep="\n")
    
    inputMeans.cont <- glue::glue("# mean values",
                       "meanlev <- c({input$meanlev})",
                       "resSD <- {input$resSD}",
                       "",
                       "gen <- FixedMean(n, doselev, meanlev, resSD)",
                       .sep="\n")

    inputMeans.bin <- glue::glue("# proportions",
                            "meanlev <- c({input$meanlev})",
                            "gen <- FixedMean(n, doselev, meanlev, binary = TRUE)",
                            .sep="\n")
    
    param.cont <- glue::glue("resSD <- {input$resSD}",
                       "",
                       "meanlev <- emaxfun(doselev, pop)",
                       "gen <- FixedMean(n, doselev, meanlev, resSD, parm = pop)",
                       .sep="\n")
    
    param.bin <- glue::glue("meanlev <- plogis(emaxfun(doselev, pop))",
                      "gen <- FixedMean(n, doselev, meanlev, parm = pop, binary = TRUE)",
                      .sep="\n")
    
    binary <- input$binary
    input2 <- ifelse(binary, inputMeans.bin, inputMeans.cont)
    param2 <- ifelse(binary, 
                     glue::glue(param, param.bin, .sep="\n"), 
                     glue::glue(param, param.cont, .sep="\n"))
    
    simInputs <- ifelse(input$inputVals==1, param2, input2)
    
    ## Simulate for Maximum Likelihood estimation
    emaxsim <- glue::glue("# Simulate outcomes",
                    "D1 <- emaxsim(nsim, gen, modType={input$modType}, nproc={input$nprocs}, binary = {binary})",
                    "",
                    "summary(D1, testalph=0.05)",
                    "plot(D1)",
                    "", .sep="\n")
    
    ## Simulate for Bayesian estimation
    prior.cont <-  glue::glue("# Priors",
                        "prior <- emaxPrior.control(epmu = {input$epmu},",
                        "\tepsca = {input$epsca},",
                        "\tdifTargetmu = {input$difTargetmu},",
                        "\tdifTargetsca = {input$difTargetsca},",
                        "\tdTarget = {input$dTarget},",
                        "\tp50 = {input$p50},",
                        "\tsigmalow = {input$sigmalow},",
                        "\tsigmaup = {input$sigmaup},",
                        "\tparmDF = 5)\n",
                        .sep="\n")
    
    prior.bin <- glue::glue("# Priors",
                      "prior <- emaxPrior.control(epmu = {input$epmu},",
                      "\tepsca = {input$epsca},",
                      "\tdifTargetmu = {input$difTargetmu},",
                      "\tdifTargetsca = {input$difTargetsca},",
                      "\tdTarget = {input$dTarget},",
                      "\tp50 = {input$p50},",
                      "\tparmDF = 5,",
                      "\tbinary = TRUE)\n",
                      .sep="\n")
    
    prior <- ifelse(binary, prior.bin, prior.cont)
    
    mcmc <- glue::glue("# Stan and MCMC settings",
                 "rstan::rstan_options(auto_write = TRUE)",
                 "mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,",
                 "\tpropInit=0.15,adapt_delta = 0.95)",
                 .sep="\n")
    
    ##NT changed nproc
    emaxsimB <- glue::glue("# Simulate outcomes",
                     "D1 <- emaxsimB(nsim, gen, prior,",
                     "\tmodType={input$modType},mcmc=mcmc,",
                     "\tcheck=FALSE, nproc={input$nprocs}, binary = {binary})",
                     "summary(D1)",
                     "plot(D1)", .sep="\n")
    
    if(input$modType=="3" ){ est <- "D1$est3" }
    if(input$modType=="4" ){ est <- "D1$est4" }
    if(input$bayes) { est <- "D1$est" }
    
    bayes <- ifelse(input$bayes,
                    glue::glue(prior, mcmc, .sep="\n"),
                    "")
    emaxsim <- ifelse(input$bayes,
                      emaxsimB,
                      emaxsim)
    plot1 <- glue::glue("# Additional plots of simulation results",
                  "library(tidyverse,quietly = TRUE)",
                  " ",
                  "{est} %>%",
                  "as_tibble() %>%",
                  "mutate(ED50 = exp(led50)) %>%",
                  "select(-led50) %>%",
                  .sep="\n")
    plot1.bin <- glue::glue(
                      "pairs(.)",                     
                      .sep="\n")
    plot1.cont <- glue::glue("rename(E0 = e0, ",
                       "\tEmax = emax) %>%",
                       "mutate('Residual SD' = D1$residSD) %>%",
                       "pairs(.)",
                       .sep="\n")
    plot1 <- ifelse(binary,
                    glue::glue(plot1, plot1.bin,  .sep="\n"),
                    glue::glue(plot1, plot1.cont, .sep="\n"))

    glue::glue(design,
         " ",
         simInputs,
         " ",
         bayes,
         " ",
         emaxsim,
         " ",
         plot1,
         .sep="\n")
  })
}


shinyApp(ui, server)