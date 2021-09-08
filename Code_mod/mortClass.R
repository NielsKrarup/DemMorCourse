
library("R6")

###
###  mortClass
###

mortClass <- R6Class("mortClass", private = list(.histData = NULL, .par = NULL, .m = NULL, .s = NULL, .method = NULL), public = list())

# -----------------------------------------------------------------------------
#           Accessor functions to private fields
# -----------------------------------------------------------------------------

mortClass$set("active", "par"      ,  function(value) { if (missing(value)) { private$.par } else { stop("par is read-only.\n" , call.=FALSE) } } )
mortClass$set("active", "histData" ,  function(value) { if (missing(value)) { private$.histData } else { stop("histData is read-only.\n" , call.=FALSE) } } )

# -----------------------------------------------------------------------------
#           Constructor
# -----------------------------------------------------------------------------

# Input 'OEdata'
mortClass$set("public", "initialize", function(OEdata, method = "Lee-Carter", ...) {
    if (class(OEdata) != "OEdata") stop("First argument should be of type 'OEdata'.")
    private$.histData <- OEdata

    if (method != "Lee-Carter") stop("Only Lee-Carter is implemented")
    private$.method <- method
    
    cat("----------------- Estimating Lee-Carter model -----------------\n")
    cat("\n -  Ages:", OEdata$agelim[1],"-",OEdata$agelim[2])
    cat("\n - Years:", OEdata$timelim[1],"-",OEdata$timelim[2])
    cat("\n")
    for (gidx in seq_len(length(OEdata$groups))) {
        cat("\n - Estimating parameters for group", paste0("'",OEdata$groups[gidx],"'"), ".... ")
        tmp <- lc.estimate(OEdata$O[gidx,,], OEdata$E[gidx,,], ...)
        private$.par[[OEdata$groups[gidx]]] <- list("alpha" = tmp$alpha, "beta" = tmp$beta, "kappa" = tmp$kappa)
        private$.m[[OEdata$groups[gidx]]]   <- mean(diff(tmp$kappa))
        private$.s[[OEdata$groups[gidx]]]   <- sd(diff(tmp$kappa))
        cat("\n")
    }
    
    cat("\n--------------------- Estimation complete ---------------------\n")
    invisible(self)
})

# -----------------------------------------------------------------------------
#           Method: Print method
# -----------------------------------------------------------------------------

mortClass$set("public", "print", function(...) {
    cat("Mortality object based on the", private$.method,"method.\n")
    cat("The model is estimated for", length(private$.histData$groups), "group(s):", paste(private$.histData$groups, collapse = ", "),'\n')
    cat("The model is estimated on data for ages ",private$.histData$agelim[1],"-",private$.histData$agelim[2]," and years ",private$.histData$timelim[1],"-",private$.histData$timelim[2],".\n",sep='') 
    invisible(self)
})

# -----------------------------------------------------------------------------
#           Methods for argument checking
# -----------------------------------------------------------------------------

mortClass$set("private", ".verifyGroupName", function(groups) {
    if (is.null(groups)) {groups <- private$.histData$groups}
    else {
        if (!is.character(groups)) stop("Argument 'groups' must be of type character.", call. = F)
        groups <- properCase(groups)
        if (!all(groups %in% private$.histData$groups)) stop("Argument 'groups' must be one (or more) of the following:\n    ",paste(private$.histData$groups,collapse=", "), call. = F)
    }
    private$.histData$groups[private$.histData$groups %in% groups] # Flip groups so they always appear in the order they were loaded
})

mortClass$set("private", ".verifyAgelim", function(agelim) {
    if (is.null(agelim)) {agelim <- private$.histData$agelim}
    else { 
        if (!is.numeric(agelim)) stop("Argument 'agelim' must be of type numeric.", call. = F)
        if (length(agelim) > 2) stop("Argument 'agelim' must be of length 1 or 2.", call. = F)
        agelim  <- lim2lim(agelim)
        if(agelim[1] < private$.histData$agelim[1] | agelim[2] > private$.histData$agelim[2]) stop(paste(c("Ages must be in the interval ", private$.histData$agelim[1]," to ", private$.histData$agelim[2])), call. = F)
    }
    agelim
})

mortClass$set("private", ".verifyTimelim", function(timelim) {
    if (is.null(timelim)) {timelim <- private$.histData$timelim}
    else {
        if (!is.numeric(timelim)) stop("Argument 'timelim' must be of type numeric.", call. = F)
        if (length(timelim) > 2) stop("Argument 'timelim' must be of length 1 or 2.", call. = F)
        timelim <- lim2lim(timelim)
        if(timelim[1] < private$.histData$timelim[1]) stop(paste(c("Years must be beyond ", private$.histData$timelim[1])), call. = F)
    }
    timelim
})

# -----------------------------------------------------------------------------
#           Method: getMu
# -----------------------------------------------------------------------------

# If historic = T -> Return O/E-rates for years where data exists, rather than Lee-Carter estimate
mortClass$set("public", "getMu", function(groups = NULL, agelim = NULL, timelim = NULL, historic = F) {
    stopifnot(is.logical(historic))
    groups    <- private$.verifyGroupName(groups)
    agelim    <- private$.verifyAgelim(agelim)
    timelim   <- private$.verifyTimelim(timelim)
    
    gidx <- private$.histData$groups %in% groups
    mu   <- array(NA, dim  = c(length(groups), lim2length(agelim), lim2length(timelim)), dimnames = list("Group"=groups, "Age"=lim2vec(agelim), "Year"=lim2vec(timelim)))
    aidx <- lim2vec(private$.histData$agelim) %in% lim2vec(agelim)
    tidx <- lim2vec(private$.histData$timelim) %in% lim2vec(timelim)
    
    if (historic) {
        t0 <- 1
        t1 <- lim2length(timelim)
        if (timelim[1L]<=private$.histData$timelim[2L]) {
            mu.hist <- private$.histData$O/private$.histData$E
            mu.hist[private$.histData$E == 0] <- NaN
            mu[,,t0:sum(tidx)] <- mu.hist[gidx,aidx,tidx]
            t0 <- t0 + sum(tidx)
        }
        if (timelim[2L]>private$.histData$timelim[2L]) {
            tvec <- lim2vec(pmax(1,timelim-private$.histData$timelim[2]))
            for (group in groups) {
                par   <- private$.par[[group]]
                alpha <- par$alpha[aidx]
                beta  <- par$beta[aidx]
                kappa <- tail(par$kappa,1)+tvec*private$.m[[group]]
                mu[group,,t0:t1] <- lc.par2mu(alpha,beta,kappa)
            }
        }
    }
    else {
        tvec <- if (timelim[2]>private$.histData$timelim[2]) lim2vec(pmax(1,timelim-private$.histData$timelim[2])) else numeric(0)
        for (group in groups) {
            par   <- private$.par[[group]]
            alpha <- par$alpha[aidx]
            beta  <- par$beta[aidx]
            kappa <- c(par$kappa[tidx], tail(par$kappa,1)+tvec*private$.m[[group]])
            mu[group,,] <- lc.par2mu(alpha,beta,kappa)
        }
    }
    structure(list(groups = groups, agelim = agelim, timelim = timelim, mu = mu), class = "surf")
})

# -----------------------------------------------------------------------------
#           Method: plotdat
# -----------------------------------------------------------------------------

# Returns various statistics in a data frame with the format: Age Year Stat1 ... StatN
mortClass$set("public", "plotdat", function(groups=NULL, agelim=NULL, timelim=NULL) {
    groups    <- private$.verifyGroupName(groups)
    agelim    <- private$.verifyAgelim(agelim)
    timelim   <- private$.verifyTimelim(timelim)
    
    agevec  <- lim2vec(agelim)
    yearvec <- lim2vec(timelim)
    df.est  <- data.frame( Group  = rep(groups, times = length(agevec)*length(yearvec))
                         , Age    = rep(agevec, each = length(groups), times = length(yearvec))
                         , Year   = rep(yearvec, each = length(groups)*length(agevec))
                         , mu.est = as.vector(self$getMu(groups,agelim,timelim)$mu))
    
    ageHist  <- agevec[agevec %in% lim2vec(private$.histData$agelim)]
    yearHist <- yearvec[yearvec %in% lim2vec( private$.histData$timelim)]
    
    mu.hist <- private$.histData$O/private$.histData$E
    mu.hist[private$.histData$E == 0] <- NaN
    mu.hist <- mu.hist[private$.histData$groups %in% groups,lim2vec(private$.histData$agelim) %in% agevec, lim2vec(private$.histData$timelim) %in% yearvec]
    
    df.hist  <- data.frame(Group   = rep(groups, times = length(ageHist)*length(yearHist))
                         , Age     = rep(ageHist, each = length(groups), times = length(yearHist))
                         , Year    = rep(yearHist, each = length(groups)*length(ageHist))
                         , mu.hist = as.vector(mu.hist))
    
    merge(df.hist, df.est, by = c("Group","Age", "Year"), all = T)
})

mortClass$set("public", "plotdatEstimates", function(type = c("alpha", "beta", "kappa")) {
    type <- match.arg(type)
    
    ageOrYr    <- ifelse(type == "kappa", "Year", "Age")
    ageOrYrVec <- if (type == "kappa") lim2vec(private$.histData$timelim) else lim2vec(private$.histData$agelim)
    
    ret <- data.frame()
    for (group in private$.histData$groups) {
        par <- private$.par[[group]][[type]]
        string <- c("data.frame('Group' = group,", ageOrYr, "= ageOrYrVec,", properCase(type), "=par)")
        tmp <- eval(parse(text=paste(string, collapse="")))
        ret <- rbind(ret,tmp)
    }
    ret
})

