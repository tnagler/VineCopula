#########################################################
# C2RVine    					#
#							#
# Input:						#
# family	vector of copula families		#
# par		vector of the (first) copula parameters	#
# par2		vector of the (second) copula parameters#
#							#
# Output:						#
# RVM		RVM-object				#
#########################################################

C2RVine <- function(order, family, par, par2 = rep(0, length(family))) {
    dd <- length(family)
    if (dd < 1) 
        stop("Length of 'family' has to be positive.")
    if (length(par) != length(par2)) 
        stop("Lengths of 'par' and 'par2' do not match.")
    if (length(par) != dd) 
        stop("Lengths of 'family' and 'par' do not match.")
    
    for (i in 1:dd) {
        if (!(family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40))) 
            stop("Copula family not implemented.")
        # Parameterbereiche abfragen
        if ((family[i] == 1 || family[i] == 2) && abs(par[i]) >= 1) 
            stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
        if (family[i] == 2 && par2[i] <= 2) 
            stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
        if ((family[i] == 3 || family[i] == 13) && par[i] <= 0) 
            stop("The parameter of the Clayton copula has to be positive.")
        if ((family[i] == 4 || family[i] == 14) && par[i] < 1) 
            stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
        if ((family[i] == 6 || family[i] == 16) && par[i] <= 1) 
            stop("The copula parameter of the Joe copula has to be in the interval (1,oo).")
        if (family[i] == 5 && par[i] == 0) 
            stop("The parameter of the Frank copula has to be unequal to 0.")
        if ((family[i] == 7 || family[i] == 17) && par[i] <= 0) 
            stop("The first parameter of the BB1 copula has to be positive.")
        if ((family[i] == 7 || family[i] == 17) && par2[i] < 1) 
            stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
        if ((family[i] == 8 || family[i] == 18) && par[i] <= 0) 
            stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
        if ((family[i] == 8 || family[i] == 18) && par2[i] < 1) 
            stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
        if ((family[i] == 9 || family[i] == 19) && par[i] < 1) 
            stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
        if ((family[i] == 9 || family[i] == 19) && par2[i] <= 0) 
            stop("The second parameter of the BB7 copula has to be positive.")
        if ((family[i] == 10 || family[i] == 20) && par[i] < 1) 
            stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
        if ((family[i] == 10 || family[i] == 20) && (par2[i] <= 0 || par2[i] > 1)) 
            stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
        if ((family[i] == 23 || family[i] == 33) && par[i] >= 0) 
            stop("The parameter of the rotated Clayton copula has to be negative.")
        if ((family[i] == 24 || family[i] == 34) && par[i] > -1) 
            stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
        if ((family[i] == 26 || family[i] == 36) && par[i] >= -1) 
            stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
        if ((family[i] == 27 || family[i] == 37) && par[i] >= 0) 
            stop("The first parameter of the rotated BB1 copula has to be negative.")
        if ((family[i] == 27 || family[i] == 37) && par2[i] > -1) 
            stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 28 || family[i] == 38) && par[i] >= 0) 
            stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 28 || family[i] == 38) && par2[i] > -1) 
            stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 29 || family[i] == 39) && par[i] > -1) 
            stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 29 || family[i] == 39) && par2[i] >= 0) 
            stop("The second parameter of the rotated BB7 copula has to be negative.")
        if ((family[i] == 30 || family[i] == 40) && par[i] > -1) 
            stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 30 || family[i] == 40) && (par2[i] >= 0 || par2[i] < (-1))) 
            stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
    }
    
    d <- (1 + sqrt(1 + 8 * dd))/2
    
    if (length(order) != d) 
        stop("Length of 'order' and dimension of the C-vine do not match.")
    
    Matrix <- matrix(rep(0, d * d), d, d)
    Copula.Params <- matrix(rep(0, d * d), d, d)
    Copula.Params2 <- matrix(rep(0, d * d), d, d)
    Copula.Types <- matrix(rep(0, d * d), d, d)
    
    # structure
    for (i in 1:d) {
        for (j in 1:(d - i + 1)) {
            Matrix[(d - i + 1), j] <- order[i]
        }
    }
    
    # copula properties
    k <- 1
    for (i in 1:(d - 1)) {
        for (j in 1:(d - i)) {
            Copula.Types[(d - i + 1), (d - j - i + 1)] <- family[k]
            Copula.Params[(d - i + 1), (d - j - i + 1)] <- par[k]
            Copula.Params2[(d - i + 1), (d - j - i + 1)] <- par2[k]
            k <- k + 1
        }
    }
    
    RVM <- RVineMatrix(Matrix = Matrix,
                       family = Copula.Types, 
                       par = Copula.Params, 
                       par2 = Copula.Params2)
    
    return(RVM)
}
