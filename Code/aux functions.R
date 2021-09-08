###
###  Load and/or install package
###

loadPackage <- function(package,install=T) {
    package   <- as.character(substitute(package))
    installed <- any(package==.packages(all.available = TRUE))
    if (!installed && install) {
        install.packages(package)
    }
    if(suppressMessages(suppressWarnings(require(package, character.only = TRUE)))) {
        return(invisible(T))
    }
    else {
        warning("\n Failed to install/load: ", paste0("'",package,"'"))
        return(suppressMessages(invisible(require(package, character.only = TRUE))))
    }
}

###
###  Ranges
###

lim2vec <- function(lim) {lim[1L]:lim[2L]}

lim2length <- function(lim) {lim[2L]-lim[1L]+1L}

lim2lim <- function(lim) {
    if (length(lim) == 1) return(c(lim,lim))
    if (lim[1L] > lim[2L]) return(rev(lim))
    else return(lim)
}

###
###
###

logit <- function(p) {log(p/(1-p))}
expit <- function(a) {exp(a)/(exp(a)+1)}

dskewnormal <- function(xvec,d,xi,omega,alpha) { d*(2/omega)*dnorm((xvec-xi)/omega)*pnorm(alpha*(xvec-xi)/omega) }

# Fit cubic polynomial P(x) = a + b(x-x1) + c(x-x1)^2 + d(x-x1)^3 satisfying
# P(x1) = y1, P'(x1) = slope1, P(x2) = y2, P'(x2) = slope2
# It is required that x1 != x2
fitCubicSpline <- function(x1,y1,slope1,x2,y2,slope2) {
  
  delta <- x2-x1
  if (delta == 0) stop('x1 and x2 must differ.')

  a    <- y1
  b    <- slope1
  tmp1 <- (y2-a-b*delta)/delta^2
  tmp2 <- (slope2-b)/delta
  c    <- 3*tmp1 - tmp2
  d    <- (tmp2-2*tmp1)/delta
  
  function(x) { a + b*(x-x1) + c*(x-x1)^2 + d*(x-x1)^3} 
}

###
###  Return propercase of a string
###

properCase <- function(s) sub("(.)", ("\\U\\1"), tolower(s), pe = TRUE)

###
###  Life expectancies
###

life.exp.base <- function(musurf, birthdate, calcdate, maxage = NULL, type = c("period", "cohort")) {
    if (is.null(maxage)) maxage <- musurf$agelim[2]
    type <- match.arg(type)
    
    calcdate.year <- floor(calcdate)  # year in which calculation of life expectancy starts
    todate        <- birthdate+maxage
    todate.year   <- ceiling(todate)-1   # last year we need to visit to calculate life expectancy
    fromage       <- calcdate-birthdate
    fromage.year  <- floor(fromage)
    
    if (calcdate < birthdate)               stop("Input 'calcdate' must be larger than input 'birthdate'.\n")
    if (maxage > (musurf$agelim[2L]+1L))    stop("Maxage not within age span of musurf. Maxage can be at most ",musurf$agelim[2L]+1L,"\n")
    if (calcdate.year < musurf$timelim[1L]) stop("Calcdate, ",calcdate,", not within time range of musurf. Calcdate must be at least ",musurf$timelim[1L],".\n")
    if (fromage.year < musurf$agelim[1L])   stop("Age at calcdate, ",fromage,", not within age range of surf. Age at fromdate must be at least ",musurf$agelim[1L],".\n")
    if ((todate.year > musurf$timelim[2L]) && (type == "cohort")) warning("\n\n  Date at which maxage is attained, ",todate,", is not within time range of surf. Calculation performed on boundary year.\n")
    
    # Argument checking complete - calculate life expectancy
    life <- 0
    if (maxage > fromage) {
        # Find intensity "pieces"
        agevec <- fromage
        muvec <- NULL
    
        ageoffset   <- musurf$agelim[1L]-1L
        yearoffset  <- musurf$timelim[1L]-1L
        curdate     <- calcdate
        fromyearidx <- floor(calcdate)-yearoffset
        
        mu <- musurf$mu
        while (curdate < todate) {
            curdate.year <- floor(curdate)
            curage       <- curdate-birthdate
            curage.year  <- floor(curage)
            # Find next date at which either age or year boundary is crossed
            nextdate <- min(todate, curdate + min(curdate.year+1-curdate, curage.year+1-curage))
            agevec   <- c(agevec, nextdate-birthdate)
            if (type == "cohort") {
                muvec <- c(muvec,mu[curage.year-ageoffset, min(curdate.year-yearoffset, ncol(mu))])
            } 
            else {
                muvec <- c(muvec,mu[curage.year-ageoffset, fromyearidx])
            }    
            curdate <- nextdate
        }
        
        life <- NaN
        if (all(!is.nan(muvec))) {
            # Calculate life expectancy as the integrated survival probability
            # Due to the intensity being piecewise constant the integral can be calculated analytically
            dagevec <- diff(agevec)
            surv1   <- exp(-dagevec*muvec)
            surv    <- c(1, cumprod(surv1[-length(surv1)]))
            life    <- sum(surv*ifelse(muvec > 0, (1-surv1)/muvec, dagevec))
        }
    }
    life
}

surf2lifeexp <- function(surf, groups = NULL, timelim = NULL, agelim = NULL, type = c("period","cohort"), maxage = NULL) {
    if (class(surf)[1L] != "surf") stop("Class attribute '",class(surf),"' is not valid. First argument must be of class 'surf'.\n")
    if (is.null(groups)) groups <- surf$groups
    else {
        if (!all(groups %in% surf$groups)) stop("The following groups are not contained within 'surf': ", paste(groups[!(groups %in% surf$groups)], sep = ", "))
        groups <- surf$groups[surf$groups %in% groups]
    }
    if (is.null(agelim))  agelim  <- surf$agelim
    else {
        agelim  <- lim2lim(agelim)
        if(agelim[1] < surf$agelim[1] | agelim[2] > surf$agelim[2]) stop(paste(c("Ages must be in the interval ", surf$agelim[1]," to ", surf$agelim[2])))
    }
    if (is.null(timelim)) timelim <- surf$timelim
    else {
        timelim <- lim2lim(timelim)
        if(timelim[1] < surf$timelim[1] | timelim[2] > surf$timelim[2]) stop(paste(c("Years must be in the interval ", surf$timelim[1]," to ", surf$timelim[2])))
    }
    type <- match.arg(type)
    lifeexp <- array(NA, dim = c(length(groups), lim2length(agelim), lim2length(timelim)), dimnames = list("Group"=groups, "Age"=lim2vec(agelim), "Year"=lim2vec(timelim)))
    for (group in groups) {
        tmp <- list(agelim=surf$agelim,timelim=surf$timelim,mu=as.matrix(surf$mu[group,,]))
        for (year in lim2vec(timelim)) {
            calcdate <- year + 0.5
            for (age in lim2vec(agelim)) {
                birthdate <- calcdate - age
                lifeexp[group,as.character(age),as.character(year)] <- life.exp.base(tmp, birthdate, calcdate, maxage, type)
            }
        }
    }
    structure(list(groups = groups, agelim = agelim, timelim = timelim, lifeexp = lifeexp), type = type, class = "lifexp")
}

###
###  surf2val / val2surf
###

mu2surf <- function(groups, agelim, timelim, mu) {
    stopifnot(is.array(mu))
    stopifnot(dim(mu)[1] == length(groups))
    stopifnot(dim(mu)[2] == lim2length(agelim))
    stopifnot(dim(mu)[3] == lim2length(timelim))
    dimnames(mu) <- list("Group" = groups, "Age" = lim2vec(agelim), "Year" = lim2vec(timelim))
    structure(list(groups = groups, agelim = agelim, timelim = timelim, mu = mu), class='surf')
}

fer2surf <- function(agelim, timelim, fer) {
    stopifnot(is.matrix(fer))
    stopifnot(nrow(fer) == lim2length(agelim))
    stopifnot(ncol(fer) == lim2length(timelim))
    dimnames(fer) <- list("Age" = lim2vec(agelim), "Year" = lim2vec(timelim))
    structure(list(agelim = agelim, timelim = timelim, fer = fer), class='surf.fer')
}

surf2mu <- function(surf, groups = NULL, agelim = NULL, timelim = NULL) {
    if (class(surf)[1L] != "surf") stop("Class attribute '",class(surf),"' is not valid. First argument must be of class 'surf'.\n")
    if (is.null(groups)) groups <- surf$groups
    else {
        if (!all(groups %in% surf$groups)) stop("The following groups are not contained within 'surf': ", paste(groups[!(groups %in% surf$groups)], sep = ", "))
        groups <- surf$groups[surf$groups %in% groups]
    }
    if (is.null(agelim)) agelim  <- surf$agelim
    else {
        agelim  <- lim2lim(agelim)
        if(agelim[1] < surf$agelim[1] | agelim[2] > surf$agelim[2]) stop(paste(c("Ages must be in the interval ", surf$agelim[1]," to ", surf$agelim[2])))
    }
    if (is.null(timelim)) timelim <- surf$timelim
    else {
        timelim <- lim2lim(timelim)
        if(timelim[1] < surf$timelim[1] | timelim[2] > surf$timelim[2]) stop(paste(c("Years must be in the interval ", surf$timelim[1]," to ", surf$timelim[2])))
    }
    gidx <- surf$groups %in% groups
    aidx <- lim2vec(surf$agelim) %in% lim2vec(agelim)
    tidx <- lim2vec(surf$timelim) %in% lim2vec(timelim)
    surf$mu[gidx, aidx, tidx,drop=F]
}

surf2fer <- function(surf.fer, agelim = NULL, timelim = NULL) {
    if (class(surf.fer)[1L] != "surf.fer") stop("Class attribute '",class(surf.fer),"' is not valid. First argument must be of class 'surf.fer'.\n")
    if (is.null(agelim)) agelim  <- surf.fer$agelim
    else {
        agelim  <- lim2lim(agelim)
        if(agelim[1] < surf.fer$agelim[1] | agelim[2] > surf.fer$agelim[2]) stop(paste(c("Ages must be in the interval ", surf.fer$agelim[1]," to ", surf.fer$agelim[2])))
    }
    if (is.null(timelim)) timelim <- surf.fer$timelim
    else {
        timelim <- lim2lim(timelim)
        if(timelim[1] < surf.fer$timelim[1] | timelim[2] > surf.fer$timelim[2]) stop(paste(c("Years must be in the interval ", surf.fer$timelim[1]," to ", surf.fer$timelim[2])))
    }
    aidx <- lim2vec(surf.fer$agelim) %in% lim2vec(agelim)
    tidx <- lim2vec(surf.fer$timelim) %in% lim2vec(timelim)
    surf.fer$fer[aidx, tidx,drop=F]
}


# -----------------------------------------------------------------------------
#           LC-estimation
# -----------------------------------------------------------------------------

# Function to estimate the Poisson version of Lee-Carter model, cf. Brouhns et. al (2002)
lc.estimate <- function(Omat,Emat, maxiter = 5000L, tol = 1e-10) {
    alpha <- rep(0,len=nrow(Omat))
    beta  <- rep(1/nrow(Omat),len=nrow(Omat))
    kappa <- rep(0,len=ncol(Omat))
    ones  <- rep(1,len=ncol(Omat))
    
    logmu        <- matrix(alpha,ncol=1) %*% matrix(ones,nrow=1) +  matrix(beta,ncol=1) %*% matrix(kappa,nrow=1)
    old.log.like <- sum(Omat*logmu - Emat*exp(logmu))
    
    iter    <- 0L
    delta   <- 1 
    while ((iter < maxiter) && (delta > tol)) {
        iter <- iter + 1L

        # Update alpha
        Dhat     <- Emat*exp(matrix(alpha,ncol=1) %*% matrix(ones,nrow=1) + matrix(beta,ncol=1) %*% matrix(kappa,nrow=1))
        alpha    <- alpha + rowSums(Omat-Dhat)/rowSums(Dhat)

        # Update kappa
        Dhat     <- Emat*exp(matrix(alpha,ncol=1) %*% matrix(ones,nrow=1) + matrix(beta,ncol=1) %*% matrix(kappa,nrow=1))
        kappa    <- kappa + (matrix(beta,nrow=1) %*% (Omat - Dhat)) / (matrix(beta^2,nrow=1) %*% Dhat)
        meankap  <- mean(kappa)
        kappa    <- kappa - meankap
        alpha    <- alpha + meankap * beta
        
        # Update beta
        Dhat     <- Emat*exp(matrix(alpha,ncol=1) %*% matrix(ones,nrow=1) + matrix(beta,ncol=1) %*% matrix(kappa,nrow=1))
        beta     <- beta + ((Omat - Dhat) %*% matrix(kappa, ncol=1)) / ( Dhat %*% matrix(kappa^2,ncol=1) )
        sumbeta  <- sum(beta)
        beta     <- beta/sumbeta
        kappa    <- kappa*sumbeta

        logmu    <- matrix(alpha,ncol=1) %*% matrix(ones,nrow=1) +  matrix(beta,ncol=1) %*% matrix(kappa,nrow=1)
        log.like <- sum(Omat*logmu - Emat*exp(logmu))
        delta    <- log.like - old.log.like
        old.log.like <- log.like
    }
    if(iter==maxiter) cat("\n   > Estimation did NOT converge in", maxiter, "iterations. \n   > Last change in log-likelihood =", delta)
    else cat("\n   > Estimation converged after", iter, "iterations.")
    flush.console()
    list(alpha=c(alpha),beta=c(beta),kappa=c(kappa))
}

lc.par2mu <- function(alpha, beta, kappa) {
    stopifnot(is.vector(alpha), is.vector(beta), is.vector(kappa))
    stopifnot(is.numeric(alpha), is.numeric(beta), is.numeric(kappa))
    stopifnot(length(alpha) == length(beta))
    exp(matrix(alpha,ncol=1) %*% matrix(1,ncol = length(kappa)) + matrix(beta,ncol=1) %*% matrix(kappa,nrow=1))
}