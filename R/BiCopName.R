#' Bivariate Copula Family Names
#'
#' This function transforms the bivariate copula family number into its
#' character expression and vice versa.
#'
#'
#' @param family Bivariate copula family, either its number or its character
#' expression (see table below).
#' \tabular{rll}{
#' No.        \tab Short name         \tab Long name \cr
#' `0`   \tab `"I"`         \tab `"Independence"` \cr
#' `1`   \tab `"N"`         \tab `"Gaussian"` \cr
#' `2`   \tab `"t"`         \tab `"t"` \cr
#' `3`   \tab `"C"`         \tab `"Clayton"` \cr
#' `4`   \tab `"G"`         \tab `"Gumbel"` \cr
#' `5`   \tab `"F"`         \tab `"Frank"` \cr
#' `6`   \tab `"J"`         \tab `"Joe"` \cr
#' `7`   \tab `"BB1"`       \tab `"BB1"` \cr
#' `8`   \tab `"BB6"`       \tab `"BB6"` \cr
#' `9`   \tab `"BB7"`       \tab `"BB7"` \cr
#' `10`  \tab `"BB8"`       \tab `"Frank-Joe"` \cr
#' `13`  \tab `"SC"`        \tab `"Survival Clayton"` \cr
#' `14`  \tab `"SG"`        \tab `"Survival Gumbel"` \cr
#' `16`  \tab `"SJ"`        \tab `"Survival Joe"` \cr
#' `17`  \tab `"SBB1"`      \tab `"Survival BB1"` \cr
#' `18`  \tab `"SBB6"`      \tab `"Survival BB6"` \cr
#' `19`  \tab `"SBB7"`      \tab `"Survival BB7"` \cr
#' `20`  \tab `"SBB8"`      \tab `"Survival BB8"` \cr
#' `23`  \tab `"C90"`       \tab `"Rotated Clayton 90 degrees"` \cr
#' `24`  \tab `"G90"`       \tab `"Rotated Gumbel 90 degrees"` \cr
#' `26`  \tab `"J90"`       \tab `"Rotated Joe 90 degrees"` \cr
#' `27`  \tab `"BB1_90"`    \tab `"Rotated BB1 90 degrees"` \cr
#' `28`  \tab `"BB6_90"`    \tab `"Rotated BB6 90 degrees"` \cr
#' `29`  \tab `"BB7_90"`    \tab `"Rotated BB7 90 degrees"` \cr
#' `30`  \tab `"BB8_90"`    \tab `"Rotated Frank-Joe 90 degrees"` \cr
#' `33`  \tab `"C270"`      \tab `"Rotated Clayton 270 degrees"` \cr
#' `34`  \tab `"G270"`      \tab `"Rotated Gumbel 270 degrees"` \cr
#' `36`  \tab `"J270"`      \tab `"Rotated Joe 270 degrees"` \cr
#' `37`  \tab `"BB1_270"`   \tab `"Rotated BB1 270 degrees"` \cr
#' `38`  \tab `"BB6_270"`   \tab `"Rotated BB6 270 degrees"` \cr
#' `39`  \tab `"BB7_270"`   \tab `"Rotated BB7 270 degrees"` \cr
#' `40`  \tab `"BB8_270"`   \tab `"Rotated Frank-Joe 270 degrees"` \cr
#' `104` \tab `"Tawn"`      \tab `"Tawn type 1"` \cr
#' `114` \tab `"Tawn180"`   \tab `"Rotated Tawn type 1 180 degrees"` \cr
#' `124` \tab `"Tawn90"`    \tab `"Rotated Tawn type 1 90 degrees"` \cr
#' `134` \tab `"Tawn270"`   \tab `"Rotated Tawn type 1 270 degrees"` \cr
#' `204` \tab `"Tawn2"`     \tab `"Tawn type 2"` \cr
#' `214` \tab `"Tawn2_180"` \tab `"Rotated Tawn type 2 180 degrees"` \cr
#' `224` \tab `"Tawn2_90"`  \tab `"Rotated Tawn type 2 90 degrees"` \cr
#' `234` \tab `"Tawn2_270"` \tab `"Rotated Tawn type 2 270 degrees"` \cr
#' }
#' @param short Logical; if the number of a bivariate copula family is used and
#' `short = TRUE` (default), a short version of the corresponding
#' character expression is returned, otherwise the long version.
#'
#' @return The transformed bivariate copula family (see table above).
#'
#' @author Ulf Schepsmeier
#'
#' @seealso [RVineTreePlot()]
#'
#' @examples
#'
#' ## family number to character expression
#' family <- 1
#' BiCopName(family, short = TRUE)	 # short version
#' BiCopName(family, short = FALSE)	 # long version
#'
#' ## family character expression (short version) to number
#' family <- "C"
#' BiCopName(family)	# as number
#'
#' ## family character expression (long version) to number
#' family <- "Clayton"
#' BiCopName(family)	# as number
#'
#' ## vectors of families
#' BiCopName(1:10)    # as character expression
#' BiCopName(c("Clayton","t","J"))    # as number
#'
BiCopName <- function(family, short = TRUE) {
    stopifnot(is.logical(short))
    sapply(family, fam_name, short = short)
}

fam_name <- function(family, short) {
    fam <- NA
    ## get family name given a number
    if (is.numeric(family)) {
        if (short == TRUE) {
            ## short name
            if (family == 0)
                fam <- "I"
            if (family == 1)
                fam <- "N"
            if (family == 2)
                fam <- "t"
            if (family == 3)
                fam <- "C"
            if (family == 4)
                fam <- "G"
            if (family == 5)
                fam <- "F"
            if (family == 6)
                fam <- "J"
            if (family == 7)
                fam <- "BB1"
            if (family == 8)
                fam <- "BB6"
            if (family == 9)
                fam <- "BB7"
            if (family == 10)
                fam <- "BB8"
            if (family == 13)
                fam <- "SC"
            if (family == 14)
                fam <- "SG"
            if (family == 16)
                fam <- "SJ"
            if (family == 17)
                fam <- "SBB1"
            if (family == 18)
                fam <- "SBB6"
            if (family == 19)
                fam <- "SBB7"
            if (family == 20)
                fam <- "SBB8"
            if (family == 23)
                fam <- "C90"
            if (family == 24)
                fam <- "G90"
            if (family == 26)
                fam <- "J90"
            if (family == 27)
                fam <- "BB1_90"
            if (family == 28)
                fam <- "BB6_90"
            if (family == 29)
                fam <- "BB7_90"
            if (family == 30)
                fam <- "BB8_90"
            if (family == 33)
                fam <- "C270"
            if (family == 34)
                fam <- "G270"
            if (family == 36)
                fam <- "J270"
            if (family == 37)
                fam <- "BB1_270"
            if (family == 38)
                fam <- "BB6_270"
            if (family == 39)
                fam <- "BB7_270"
            if (family == 40)
                fam <- "BB8_270"
            if (family == 100)
                fam <- "NP"
            if (family == 41)
                fam <- "1-par AS"
            if (family == 51)
                fam <- "1-par AS180"
            if (family == 61)
                fam <- "1-par AS90"
            if (family == 71)
                fam <- "1-par AS270"
            if (family == 104)
                fam <- "Tawn"
            if (family == 114)
                fam <- "Tawn180"
            if (family == 124)
                fam <- "Tawn90"
            if (family == 134)
                fam <- "Tawn270"
            if (family == 204)
                fam <- "Tawn2"
            if (family == 214)
                fam <- "Tawn2_180"
            if (family == 224)
                fam <- "Tawn2_90"
            if (family == 234)
                fam <- "Tawn2_270"
        } else {
            ## long names
            if (family == 0)
                fam <- "Independence"
            if (family == 1)
                fam <- "Gaussian"
            if (family == 2)
                fam <- "t"
            if (family == 3)
                fam <- "Clayton"
            if (family == 4)
                fam <- "Gumbel"
            if (family == 5)
                fam <- "Frank"
            if (family == 6)
                fam <- "Joe"
            if (family == 7)
                fam <- "BB1"
            if (family == 8)
                fam <- "BB6"
            if (family == 9)
                fam <- "BB7"
            if (family == 10)
                fam <- "BB8"
            if (family == 13)
                fam <- "Survival Clayton"
            if (family == 14)
                fam <- "Survival Gumbel"
            if (family == 16)
                fam <- "Survival Joe"
            if (family == 17)
                fam <- "Survival BB1"
            if (family == 18)
                fam <- "Survival BB6"
            if (family == 19)
                fam <- "Survival BB7"
            if (family == 20)
                fam <- "Survival BB8"
            if (family == 23)
                fam <- "Rotated Clayton 90 degrees"
            if (family == 24)
                fam <- "Rotated Gumbel 90 degrees"
            if (family == 26)
                fam <- "Rotated Joe 90 degrees"
            if (family == 27)
                fam <- "Rotated BB1 90 degrees"
            if (family == 28)
                fam <- "Rotated BB6 90 degrees"
            if (family == 29)
                fam <- "Rotated BB7 90 degrees"
            if (family == 30)
                fam <- "Rotated BB8 90 degrees"
            if (family == 33)
                fam <- "Rotated Clayton 270 degrees"
            if (family == 34)
                fam <- "Rotated Gumbel 270 degrees"
            if (family == 36)
                fam <- "Rotated Joe 270 degrees"
            if (family == 37)
                fam <- "Rotated BB1 270 degrees"
            if (family == 38)
                fam <- "Rotated BB6 270 degrees"
            if (family == 39)
                fam <- "Rotated BB7 270 degrees"
            if (family == 40)
                fam <- "Rotated BB8 270 degrees"
            if (family == 100)
                fam <- "Nonparametric"  #changed Mathias
            if (family == 41)
                fam <- "1-parametric asymmetric"
            if (family == 51)
                fam <- "Rotated 1-parametric asymmetric 180 degree"
            if (family == 61)
                fam <- "Rotated 1-parametric asymmetric 90 degree"
            if (family == 71)
                fam <- "Rotated 1-parametric asymmetric 270 degree"
            if (family == 104)
                fam <- "Tawn  type 1"
            if (family == 114)
                fam <- "Rotated Tawn type 1 180 degrees"
            if (family == 124)
                fam <- "Rotated Tawn type 1 90 degrees"
            if (family == 134)
                fam <- "Rotated Tawn type 1 270 degrees"
            if (family == 204)
                fam <- "Tawn  type 2"
            if (family == 214)
                fam <- "Rotated Tawn type 2 180 degrees"
            if (family == 224)
                fam <- "Rotated Tawn type 2 90 degrees"
            if (family == 234)
                fam <- "Rotated Tawn type 2 270 degrees"
        }
    } else {
        ## get family number given a name
        if (family == "I" || family == "Independence")
            fam <- 0
        if (family == "N" || family == "Gaussian")
            fam <- 1
        if (family == "t")
            fam <- 2
        if (family == "C" || family == "Clayton")
            fam <- 3
        if (family == "G" || family == "Gumbel")
            fam <- 4
        if (family == "F" || family == "Frank")
            fam <- 5
        if (family == "J" || family == "Joe")
            fam <- 6
        if (family == "BB1" || family == "BB1")
            fam <- 7
        if (family == "BB6" || family == "BB6")
            fam <- 8
        if (family == "BB7" || family == "BB7")
            fam <- 9
        if (family == "SC" || family == "Survival Clayton")
            fam <- 13
        if (family == "SG" || family == "Survival Gumbel")
            fam <- 14
        if (family == "SJ" || family == "Survival Joe")
            fam <- 16
        if (family == "SBB1" || family == "Survival BB1")
            fam <- 17
        if (family == "SBB6" || family == "Survival BB6")
            fam <- 18
        if (family == "SBB7" || family == "Survival BB7")
            fam <- 19
        if (family == "SBB8" || family == "Survival BB8")
            fam <- 20
        if (family == "C90" || family == "Rotated Clayton 90 degrees")
            fam <- 23
        if (family == "G90" || family == "Rotated Gumbel 90 degrees")
            fam <- 24
        if (family == "J90" || family == "Rotated Joe 90 degrees")
            fam <- 26
        if (family == "BB1_90" || family == "Rotated BB1 90 degrees")
            fam <- 27
        if (family == "BB6_90" || family == "Rotated BB6 90 degrees")
            fam <- 28
        if (family == "BB7_90" || family == "Rotated BB7 90 degrees")
            fam <- 29
        if (family == "BB8_90" || family == "Rotated BB8 90 degrees")
            fam <- 30
        if (family == "C270" || family == "Rotated Clayton 270 degrees")
            fam <- 33
        if (family == "G270" || family == "Rotated Gumbel 270 degrees")
            fam <- 34
        if (family == "J270" || family == "Rotated Joe 270 degrees")
            fam <- 36
        if (family == "BB1_270" || family == "Rotated BB1 270 degrees")
            fam <- 37
        if (family == "BB6_270" || family == "Rotated BB6 270 degrees")
            fam <- 38
        if (family == "BB7_270" || family == "Rotated BB7 270 degrees")
            fam <- 39
        if (family == "BB8_270" || family == "Rotated BB8 270 degrees")
            fam <- 40
        if (family == "NP" || family == "Nonparametric")
            fam <- 100
        if (family == "1-par AS" || family == "1-parametric asymmetric")
            fam <- 41
        if (family == "1-par AS180" || family == "Rotated 1-parametric asymmetric 180 degree")
            fam <- 51
        if (family == "1-par AS90" || family == "Rotated 1-parametric asymmetric 90 degree")
            fam <- 61
        if (family == "1-par AS270" || family == "Rotated 1-parametric asymmetric 270 degree")
            fam <- 71
        if (family == "Tawn" || family == "Tawn type 1")
            fam <- 104
        if (family == "Tawn180" || family == "Rotated Tawn type 1 180 degrees")
            fam <- 114
        if (family == "Tawn90" || family == "Rotated Tawn type 1 90 degrees")
            fam <- 124
        if (family == "Tawn270" || family == "Rotated Tawn type 1 270 degrees")
            fam <- 134
        if (family == "Tawn2" || family == "Tawn type 2")
            fam <- 204
        if (family == "Tawn2_180" || family == "Rotated Tawn type 2 180 degrees")
            fam <- 214
        if (family == "Tawn2_90" || family == "Rotated Tawn type 2 90 degrees")
            fam <- 224
        if (family == "Tawn2_270" || family == "Rotated Tawn type 2 270 degrees")
            fam <- 234
    }
    if (is.na(fam))
        stop("Family not implemented.")

    fam
}
