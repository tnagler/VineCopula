#' Shiny app for bivariate copula selection
#'
#' The function starts a shiny app which visualizes copula data and allows to
#' compare it with overlays of density contours or simulated data from different
#' copula families with fitted parameters. Several specifications for the
#' margins are available.
#'
#' @param u1,u2 Data vectors of equal length with values in [0,1].
#' @param familyset Vector of bivariate copula families to select from.
#' The vector has to include at least one bivariate copula
#' family that allows for positive and one that allows for negative dependence.
#' If \code{familyset = NA} (default), selection among all possible families is
#' performed. If a vector of negative numbers is provided, selection among all
#' but \code{abs(familyset)} families is performed. Coding of bivariate copula
#' families: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param rotations If \code{TRUE}, all rotations of the families in
#' \code{familyset} are included (or substracted).
#'
#' @return A \code{\link{BiCop}} object containing the model selected by the
#' user.
#'
#' @author Matthias Killiches, Thomas Nagler
#'
#' @examples
#'# load data
#'data(daxreturns)
#'
#'# find a suitable copula family for the first two stocks
#'\donttest{fit <- BiCopCompare(daxreturns[, 1], daxreturns[, 2])}
#'
BiCopCompare <- function(u1, u2, familyset = NA, rotations = TRUE) {
    if (!requireNamespace("shiny"))
        stop("The 'shiny' package must be installed to run the app.")

    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    remove_nas,
                    check_nobs,
                    check_if_01,
                    prep_familyset,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    ## assign data to global environment
    z1 <- qnorm(u1)
    z2 <- qnorm(u2)

    ## create list of admissible families
    tau <- cor(u1, u2, method = "kendall")
    if (tau > 0) {
        allfamlst <- list("0 - Independence " = 0,
                          "1 - Gaussian " = 1,
                          "2 - Student t " = 2,
                          "3 - Clayton " = 3,
                          "4 - Gumbel " = 4,
                          "5 - Frank " = 5,
                          "6 - Joe " = 6,
                          "7 - BB1 " = 7,
                          "8 - BB6 " = 8,
                          "9 - BB7 " = 9,
                          "10 - BB8 " = 10,
                          "13 -  Clayton  (180 deg)" = 13,
                          "14 -  Gumbel  (180 deg)" = 14,
                          "16 -  Joe  (180 deg)" = 16,
                          "17 -  BB1  (180 deg)" = 17,
                          "18 -  BB6  (180 deg)" = 18,
                          "19 -  BB7  (180 deg)" = 19,
                          "20 -  BB8  (180 deg)" = 20,
                          "104 - Tawn 1 " = 104,
                          "114 -  Tawn 1  (180 deg)" = 114,
                          "204 - Tawn 2 " = 204,
                          "214 -  Tawn 2  (180 deg)" = 214)
    } else {
        allfamlst <- list("0 - Independence " = 0,
                          "1 - Gaussian " = 1,
                          "2 - Student t " = 2,
                          "5 - Frank " = 5,
                          "23 -  Clayton  (90 deg)" = 23,
                          "24 -  Gumbel  (90 deg)" = 24,
                          "26 -  Joe  (90 deg)" = 26,
                          "27 -  BB1  (90 deg)" = 27,
                          "28 -  BB6  (90 deg)" = 28,
                          "29 -  BB7  (90 deg)" = 29,
                          "30 -  BB8  (90 deg)" = 30,
                          "33 -  Clayton  (270 deg)" = 33,
                          "34 -  Gumbel  (270 deg)" = 34,
                          "36 -  Joe  (270 deg)" = 36,
                          "37 -  BB1  (270 deg)" = 37,
                          "38 -  BB6  (270 deg)" = 38,
                          "39 -  BB7  (270 deg)" = 39,
                          "40 -  BB8  (270 deg)" = 40,
                          "124 -  Tawn 1  (90 deg)" = 124,
                          "134 -  Tawn 1  (270 deg)" = 134,
                          "224 -  Tawn 2  (90 deg)" = 224,
                          "234 -  Tawn 2  (270 deg)" = 234)
    }
    famlst <- allfamlst[unlist(allfamlst) %in% familyset]

    ## gather information about fits of bivariate copulas
    comp <- BiCopEstList(u1 = u1, u2 = u2, familyset = familyset)
    lst <- list(u1 = u1, u2 = u2, comp = comp)

    ## start shiny app
    suppressWarnings(shiny::runApp(list(
        ui = shiny::fluidPage(

            # Application title
            shiny::titlePanel("Compare bivariate copulas"),

            # General settings
            shiny::wellPanel(shiny::fluidRow(
                # Level
                shiny::column(3,
                              shiny::selectInput("margins", "Margins:",
                                                 list("uniform" = "unif",
                                                      "normal" = "normal",
                                                      "exponential" = "exp",
                                                      "flipped exponential" = "flexp"),
                                                 "uniform")),
                # Display mode
                shiny::column(3,
                              shiny::selectInput("dispmod", "Display mode:",
                                                 list("simulated data" = "simdata",
                                                      "contours" = "contours",
                                                      "contours and simulated data" = "both"),
                                                 "simulated data")
                ),
                # Sample size
                shiny::column(3,
                              shiny::conditionalPanel(condition = "input.dispmod != 'contours'",
                                                      shiny::sliderInput("nsim", "Sample size:",
                                                                         min=0, max=20000, value=5000, step=1000),
                                                      width = "80%")
                              #),shiny::column(1, shiny::textOutput("")
                ),
                # Select & close
                shiny::column(3,
                              shiny::radioButtons("radio", label = "Select family:",
                                                  choices = list("A" = "1", "B" = "2",
                                                                 "C" = "3", "D" = "4"),
                                                  selected = 1, inline = TRUE),
                              shiny::actionButton("close", "select & close",
                                                  width="100%", class="btn-primary btn-lg")
                )
            )
            ),
            shiny::wellPanel(
                shiny::fluidRow(
                    # specifications
                    shiny::column(3, shiny::selectInput("fam1", "Family A:", famlst, 1)),
                    shiny::column(3, shiny::selectInput("fam2", "Family B:", famlst, 1)),
                    shiny::column(3, shiny::selectInput("fam3", "Family C:", famlst, 1)),
                    shiny::column(3, shiny::selectInput("fam4", "Family D:", famlst, 1)))
                ,
                shiny::fluidRow(
                    # plots
                    shiny::column(3, shiny::plotOutput("plot1")),
                    shiny::column(3, shiny::plotOutput("plot2")),
                    shiny::column(3, shiny::plotOutput("plot3")),
                    shiny::column(3, shiny::plotOutput("plot4")))
            ),
            shiny::wellPanel(shiny::fluidRow(
                # AIC/BIC/ll ranking plot
                shiny::column(9, shiny::plotOutput("rankingplot")
                ),

                shiny::column(3,
                              # Selection criterion
                              shiny::selectInput("selcrit",
                                                 "Selection criterion:",
                                                 list("AIC" = "AIC",
                                                      "BIC" = "BIC",
                                                      "log-likelihood" = "logLik"),
                                                 "AIC"),
                              shiny::br(),
                              # Order plots according to selected criterion
                              shiny::actionButton("order",
                                                  "sort",
                                                  width="33%",
                                                  class="btn-primary btn-lg"))
            )
            )
        ),
        server = ## server function
            function(input, output, session) {

                # needed so that 'positive' can be used in ui.R
                # plot 1
                output$plot1 <- shiny::renderPlot({
                    # determine family number and corresponding parameters
                    fam <- input$fam1
                    index <- which(comp$summary$family == fam)
                    par <- comp$models[[index]]$par
                    par2 <- comp$models[[index]]$par2

                    # plot
                    plot_function(index, input, lst)
                })
                # plot 2
                output$plot2 <- shiny::renderPlot({
                    # determine family number and corresponding parameters
                    fam <- input$fam2
                    index <- which(comp$summary$family == fam)
                    par <- comp$models[[index]]$par
                    par2 <- comp$models[[index]]$par2

                    # plot
                    plot_function(index, input, lst)
                })
                # plot 3
                output$plot3 <- shiny::renderPlot({
                    # determine family number and corresponding parameters
                    fam <- input$fam3
                    index <- which(comp$summary$family == fam)
                    par <- comp$models[[index]]$par
                    par2 <- comp$models[[index]]$par2

                    # plot
                    plot_function(index, input, lst)
                })
                # plot 4
                output$plot4 <- shiny::renderPlot({
                    # determine family number and corresponding parameters
                    fam <- input$fam4
                    index <- which(comp$summary$family == fam)
                    par <- comp$models[[index]]$par
                    par2 <- comp$models[[index]]$par2

                    # plot
                    plot_function(index, input, lst)
                })


                # AIC/BIC/ll ranking plot
                output$rankingplot <- shiny::renderPlot({
                    if (input$selcrit == "AIC"){
                        val <- sort(comp$summary$AIC)
                        cop_names <- as.character(comp$summary$family[order(comp$summary$AIC)])
                    } else if (input$selcrit == "BIC"){
                        val <- sort(comp$summary$BIC)
                        cop_names <- as.character(comp$summary$family[order(comp$summary$BIC)])
                    } else {
                        val <- sort(comp$summary$logLik, decreasing = TRUE)
                        cop_names <- as.character(comp$summary$family[order(comp$summary$logLik,
                                                                            decreasing = TRUE)])
                    }
                    barplot(height = val, names.arg = cop_names, border = TUMivory, col = TUMlightblue)
                })


                # select families according to selected ranking criterion
                shiny::observe({
                    input$order

                    # ranking of models according to selection criterion
                    if (input$selcrit == "AIC"){
                        max_ind <- comp$summary$family[order(comp$summary$AIC)]
                    } else if (input$selcrit == "BIC"){
                        max_ind <- comp$summary$family[order(comp$summary$BIC)]
                    } else {
                        max_ind <- comp$summary$family[order(comp$summary$logLik, decreasing = TRUE)]
                    }
                    if (length(max_ind) < 4)
                        max_ind <- c(max_ind, rep(max_ind[length(max_ind)], 4 - length(max_ind)))

                    # overwrite 'fam' and 'rot' for plots 1-4
                    for (i in 1:4){
                        shiny::updateSelectInput(session, paste0("fam", i),
                                                 selected = max_ind[i])
                    }
                })

                # select family and close app
                shiny::observe({
                    if(input$close > 0){
                        fam <- switch(input$radio,
                                      "1" = as.numeric(input$fam1),
                                      "2" = as.numeric(input$fam2),
                                      "3" = as.numeric(input$fam3),
                                      "4" = as.numeric(input$fam4))
                        shiny::stopApp(comp$models[[which(comp$summary$family == fam)]])
                    }
                })
            }
    )
    ))
}


## TUM colors
TUMblue <- rgb(0, 101, 189, maxColorValue=255)
TUMlightblue <- rgb(100, 160, 200, maxColorValue=255)
TUMtrlightblue <- rgb(100, 160, 200, 25, maxColorValue=255)
TUMgreen <- rgb(162, 173, 0, maxColorValue=255)
TUMorange <- rgb(227, 114, 034, maxColorValue=255)
TUMivory <- rgb(218, 215, 203, maxColorValue=255)
tr_gray <- gray(0.25, 0.5)


## functions

# generates plots according to specified settings
plot_function <- function(index, input, lst){
    # plotting input depending on specified settings
    x <- switch(input$margins,
                "unif" = lst$u1,
                "normal" = qnorm(lst$u1),
                "exp" = qexp(lst$u1),
                "flexp" = -qexp(1 - lst$u1))
    y <- switch(input$margins,
                "unif" = lst$u2,
                "normal" = qnorm(lst$u2),
                "exp" = qexp(lst$u2),
                "flexp" = -qexp(1 - lst$u2))
    xlim <- switch(input$margins,
                   "unif" = c(0,1),
                   "normal" = c(-3.5,3.5),
                   "exp" = c(0,10),
                   "flexp" = c(-10,0))
    ylim <- xlim

    plot(NULL, xlab = "", ylab = "", xlim = xlim, ylim = ylim)

    # plot simulated data (if required)
    if (input$dispmod %in% c("simdata", "both")){
        set.seed(7)
        temp <- BiCopSim(N = input$nsim, lst$comp$models[[index]])
        simdata <- switch(input$margins,
                          "unif" = temp,
                          "normal" = qnorm(temp),
                          "exp" = qexp(temp),
                          "flexp" = -qexp(1 - temp))
        points(simdata, pch = 20, col = TUMtrlightblue)
    }

    # plot actual data
    points(x, y, pch = ".", cex = 2, col = TUMblue)

    # plot contour lines (if required)
    if (input$dispmod %in% c("contours", "both")) {
        cont_margin <- switch(input$margins,
                              "unif" = "unif",
                              "normal" = "norm",
                              "exp" = "exp",
                              "flexp" = "flexp")
        suppressWarnings(contour(lst$comp$models[[index]],
                                 margins = cont_margin,
                                 lwd = 1,
                                 drawlabels = FALSE,
                                 add = TRUE))
    }

}
