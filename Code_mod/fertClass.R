
library("R6")

###
###  fertClass
###

fertClass <- R6Class("fertClass", private = list(  # private variables for internal logic
                                                .estimated      = FALSE,  # is the fertility object estimated (and ready to project)?
                                                .method         = NULL,   # String identifying which method is used for projection
                                                .par            = NULL,   # List of method specific parameters used for projection  
                                                .lastage        = NULL),   # last age group on fertility surface; "integer">0 verified on creation
                                 public = list( # public variables for historic and projected data 
                                                histF           = NULL,   # Matrix containing historic fertility data (used for estimation)
                                                histAgelim      = NULL,   # Vector of the form (minage,maxage) corresponding to data in histF
                                                histTimelim     = NULL,   # Vector of the form (minyear,maxyear) corresponding to data in histF
                                                histSkewNormPar = NULL,   # Data.frame with estimated parameters for each year (obtained by curve fitting of skewed normal); seems somewhat unstable
                                                histNormalPar   = NULL,   # Data.frame with estimated parameters for each year (obtained by statistical fitting)
                                                histMoments     = NULL))  # Data.frame with historic moments (obtained by empirical moments)

# -----------------------------------------------------------------------------
#           Accessor functions to private fields
# -----------------------------------------------------------------------------

fertClass$set("active", "par" , function(value) { if (missing(value)) { private$.par } else { stop("par is read-only.\n" , call.=FALSE) } } )
fertClass$set("active", "estimated" , function(value) { if (missing(value)) { private$.estimated } else { stop("estimated is read-only.\n" , call.=FALSE) } } )

# -----------------------------------------------------------------------------
#           Constructor
# -----------------------------------------------------------------------------

# Overrides default behaviour of $new()
fertClass$set("public", "initialize", function(fertDat, method = "Cubic") {
    if (class(fertDat) != "fertdat") stop("First argument, 'fertDat', must be of class 'fertdat'.")

    self$histF       <- fertDat$fer/1000
    self$histAgelim  <- fertDat$agelim
    self$histTimelim <- fertDat$timelim
    private$.lastage <- fertDat$agelim[2L]
    if (method != "Cubic") stop("Only cubic projection of moments is implemented")
    private$.method   <- method
    
    cat('Imported historic data\n')
    cat('  - Age range  :', fertDat$agelim[1L],'-',fertDat$agelim[2L],"\n")
    cat('  - Year range :', fertDat$timelim[1L],'-',fertDat$timelim[2L],"\n")
    
    private$.normPar()
    
    invisible(self)
})

# -----------------------------------------------------------------------------
#         Method definitions
# -----------------------------------------------------------------------------

fertClass$set("public", "print", function(...) {
    cat("Fertility object with lastage of ",private$.lastage,".\n",sep='')
    cat("Historic data loaded for ages ",self$histAgelim[1],"-",self$histAgelim[2]," and years ",self$histTimelim[1],"-",self$histTimelim[2],".\n",sep='') 
    if (private$.estimated) { cat("Projection model estimated.\n") } else { cat("No projection model estimated.\n") }
    cat("Forecasting methodology:",private$.method,"\n")
    invisible(self)
})

# -----------------------------------------------------------------------------
#         Method: normPar
# -----------------------------------------------------------------------------

# Calculate parameters for smoothing and fitting distribution (scaled, skew normal resp. scaled, normal)

fertClass$set("private", ".normPar", function() {
    if (is.null(self$histNormalPar)) {
        yearvec    <- lim2vec(self$histTimelim)
        hagevec    <- lim2vec(self$histAgelim)
        self$histSkewNormPar <- data.frame(Year=yearvec, d=NA, xi=NA, omega=NA, alpha=NA)
        self$histNormalPar   <- data.frame(Year=yearvec,  MeanAgeMother=NA, SDAgeMother=NA, BruttoFer=NA)
        self$histMoments     <- data.frame(Year=yearvec,  MeanAgeMother=NA, SDAgeMother=NA, BruttoFer=NA)
        for (yidx in 1:length(yearvec)) {
            # Empirical moments
            fer       <- self$histF[,yidx]
            bruttofer <- sum(fer)
            fernorm   <- fer/bruttofer
            mage      <- sum(hagevec*fernorm)   # it could be argued to add 1/2 year to get the right interpretation
            sdage     <- sqrt(sum(fernorm*hagevec^2) - mage^2)
            self$histMoments[yidx,2:4] <- c(mage, sdage, bruttofer)
          
            # Curve fitting where we treat dskewnormal as a "curve" with no probabilistic interpretation
            obj.func  <- function(parvec) { sum((fer - dskewnormal(hagevec,parvec[1],parvec[2],parvec[3],parvec[4]))^2)}    
            par.opt   <- optim(par=c(bruttofer, mage, sdage, 0), obj.func)$par
            self$histSkewNormPar[yidx,2:5] <- par.opt
    
            #  Statistical fitting where we think of "age of mother" as a random variable following a distribution      
            obj.func  <- function(parvec) { -sum(fernorm*log(dnorm(hagevec,mean=parvec[1],sd=parvec[2])))}
            par.opt   <- optim(par=c(mage, sdage), obj.func)$par
            self$histNormalPar[yidx,2:4] <- c(par.opt, bruttofer)
        }
    }
    invisible(self)
})

# -----------------------------------------------------------------------------
#         Method: historicSmoothF
# -----------------------------------------------------------------------------

# Returns a smoothed, historic fertility surface
fertClass$set("public", "historicSmoothF", function() {
    agelim     <- c(0,private$.lastage)
    agevec     <- lim2vec(agelim)
    yearvec    <- lim2vec(self$histTimelim)    
    empty.fer  <- matrix(0,nrow=length(agevec),ncol=length(yearvec))
    dimnames(empty.fer) <- list("Age" = agevec, "Year" = yearvec)
    empty.surf <- structure(list(agelim=agelim, timelim=self$histTimelim, fer=empty.fer),class='surf.fer')
    ret        <- empty.surf
    for (yidx in 1:length(yearvec)) {
        normpar        <- self$histNormalPar[yidx,]
        ret$fer[,yidx] <- dnorm(agevec, normpar$MeanAgeMother, normpar$SDAgeMother)*normpar$BruttoFer
    }
    ret
})

# -----------------------------------------------------------------------------
#         Methods for plotdat extraction
# -----------------------------------------------------------------------------

fertClass$set("public", "plotdatEstimates", function(type = "Moments", horizon = 50) {
    stopifnot(horizon >= 1)
    yearhist <- lim2vec(self$histTimelim)
    yearvec  <- self$histTimelim[1]:(self$histTimelim[2] + horizon)
    if (type == "Moments") {
        NAprojVec <- rep(NA, horizon)
        ret <- data.frame(Year=yearvec,
                          HistMeanAgeMother = c(self$histNormalPar$MeanAgeMother, NAprojVec),
                          HistSDAgeMother   = c(self$histNormalPar$SDAgeMother  , NAprojVec),
                          HistBruttoFer     = c(self$histNormalPar$BruttoFer    , NAprojVec),
                          ProjMeanAgeMother = NA, ProjSDAgeMother = NA, ProjBruttoFer = NA)
        ret$ProjMeanAgeMother <- private$.par[["MeanAgeMother"]](yearvec)
        ret$ProjSDAgeMother   <- private$.par[["SDAgeMother"]](yearvec)
        ret$ProjBruttoFer     <- private$.par[["BruttoFer"]](yearvec)       
    }
    else stop('Unknown type')
    ret
})

# Returns various statistics in a data frame with the format: Age Year Stat1 ... StatN
fertClass$set("public", "plotdatFer", function(agelim=NULL, timelim=NULL) {
    if (is.null(agelim)) {agelim <- self$histAgelim}
    else { 
        if (!is.numeric(agelim)) stop("Argument 'agelim' must be of type numeric.")
        if (length(agelim) > 2) stop("Argument 'agelim' must be of length 1 or 2.")
        agelim  <- lim2lim(agelim)
        if(agelim[1] < 0 | agelim[2] > self$histAgelim[2]) stop(paste(c("Ages must be in the interval ", 0," to ", self$histAgelim[2])))
    }
    if (is.null(timelim)) {timelim <- self$histTimelim}
    else {
        if (!is.numeric(timelim)) stop("Argument 'timelim' must be of type numeric.")
        if (length(timelim) > 2) stop("Argument 'timelim' must be of length 1 or 2.")
        timelim <- lim2lim(timelim)
        if(timelim[1] < self$histTimelim[1]) stop(paste(c("Years must be beyond ", self$histTimelim[1])))
    }
    
    agevec   <- lim2vec(agelim)
    yearvec  <- lim2vec(timelim)
    ageHist  <- agevec[agevec %in% lim2vec(self$histAgelim)]
    yearHist <- yearvec[yearvec %in% lim2vec(self$histTimelim)]
    
    F.hist    <- self$histF[lim2vec(self$histAgelim) %in% agevec, lim2vec(self$histTimelim) %in% yearvec]
    F.smooth  <- self$getFer(agelim, timelim)$fer
    
    df.hist   <- data.frame(Age = rep(ageHist, times = length(yearHist)), Year = rep(yearHist, each = length(ageHist)), HistF = as.vector(F.hist))
    df.smooth <- data.frame(Age = rep(agevec, times = length(yearvec)), Year = rep(yearvec, each = length(agevec)), SmoothF = as.vector(F.smooth))
    
    merge(df.hist, df.smooth, by = c("Age", "Year"), all = T)
})

fertClass$set("public", "getFer", function(agelim = NULL, timelim = NULL){
    if (is.null(agelim)) {agelim <- c(0,self$histAgelim[2])}
    else { 
        if (!is.numeric(agelim)) stop("Argument 'agelim' must be of type numeric.")
        if (length(agelim) > 2) stop("Argument 'agelim' must be of length 1 or 2.")
        agelim  <- lim2lim(agelim)
        if(agelim[1] < 0 | agelim[2] > self$histAgelim[2]) stop(paste(c("Ages must be in the interval ", 0," to ", self$histAgelim[2])))
    }
    if (is.null(timelim)) {timelim <- self$histTimelim}
    else {
        if (!is.numeric(timelim)) stop("Argument 'timelim' must be of type numeric.")
        if (length(timelim) > 2) stop("Argument 'timelim' must be of length 1 or 2.")
        timelim <- lim2lim(timelim)
        if(timelim[1] < self$histTimelim[1]) stop(paste(c("Years must be beyond ", self$histTimelim[1])))
    }

    agevec  <- 0L:private$.lastage
    yearvec <- if (timelim[2]>self$histTimelim[2]) (self$histTimelim[2]+1):timelim[2] else numeric(0)

    ret <- array(NA, dim = c(length(agevec), lim2length(range(yearvec, self$histTimelim))), dimnames = list("Age"=agevec, "Year"=lim2vec(range(yearvec, self$histTimelim))))
    ret[,1:lim2length(self$histTimelim)] <- self$historicSmoothF()$fer
    for (year in yearvec) {
        meanAge <- private$.par[["MeanAgeMother"]](year)
        sdAge   <- private$.par[["SDAgeMother"]](year)
        bruttoF <- private$.par[["BruttoFer"]](year)
        ret[,year - self$histTimelim[1] + 1] <- dnorm(agevec, meanAge, sdAge)*bruttoF
    }
    
    aidx <- agevec %in% lim2vec(agelim) 
    tidx <- lim2vec(range(yearvec, self$histTimelim)) %in% lim2vec(timelim)
    ret  <- ret[aidx,tidx,drop=F]
    
    structure(list("agelim" = agelim, "timelim" = timelim, "fer" = ret), class = "surf.fer")
})

# -----------------------------------------------------------------------------
#         Methods for estimation
# -----------------------------------------------------------------------------

# Fits a cubic spline with the specified value and slope in jumpoff year (=last year of historic data)
# and the specified value and slope=0 at the horizon (which must be strictly larger than the jumpoff year).
# The arguments must be lists with the named arguments: MeanAgeMother, SDAgeMother, BruttoFer
#
# When arguments are missing, the following defaults are used:
#   If jumpoffValue is missing, the empiric value at the jumpoff year is used instead
#   If jumpoffSlope is missing, a slope of 0 is assumed
#   If horizonValue is missing, it is set equal to jumpoffValue
#   If jumpoffSlope and horizonValue are both missing, horizon can also missing, otherwise it must be specified and larger than jumpoff year
#
# A flat curve can be achieved by specifying only jumpoffValue (with jumpoffSlope, horizonValue and horizon all missing)
fertClass$set("public", "estimate", function(jumpoffValue=list(), jumpoffSlope=list(), horizonValue=list(), horizon=list()) {

  validNames <- c('MeanAgeMother','SDAgeMother','BruttoFer')
  errorMsg   <- paste0(' must be unique and from the list: ',paste(validNames,collapse=', '))
  validArgName <- function(arg) { (length(arg)==0) || ((length(names(arg))==length(arg)) && all(names(arg) %in% validNames) && (anyDuplicated(names(arg))==0)) }

  if (!all(is.list(jumpoffValue), is.list(jumpoffSlope), is.list(horizonValue), is.list(horizon))) stop('All arguments should be lists')
  if (!validArgName(jumpoffValue)) stop(paste0('Names of jumpoffValue',errorMsg))
  if (!validArgName(jumpoffSlope)) stop(paste0('Names of jumpoffSlope',errorMsg))
  if (!validArgName(horizonValue)) stop(paste0('Names of horizonValue',errorMsg))
  if (!validArgName(horizon))      stop(paste0('Names of horizon',errorMsg))

  getVal <- function(name,arg,default) { if (name %in% names(arg)) { return(arg[[name]]) } else { return(default) } }

  prettySpace <- list(MeanAgeMother = '', SDAgeMother='  ', BruttoFer = '    ')    
  private$.par <- list()  
  for (name in validNames) {
    x1     <- self$histTimelim[2]  # Note, x1 is the same for all names. Hence, no need to use 'force(x1)'
    y1     <- getVal(name, jumpoffValue, tail(self$histNormalPar[[name]],1))
    slope1 <- getVal(name, jumpoffSlope, 0)
    y2     <- getVal(name, horizonValue, y1)
    flat   <- FALSE
    if (!(name %in% names(horizon))) {
      if ((name %in% names(jumpoffSlope)) || (name %in% names(horizonValue))) { 
        stop(paste0('jumpoffSlope and/or horizonValue are specified for ',name,' but horizon is missing.')) }
      flat <- TRUE  # flat extrapolation, since all of jumpoffSlope, horizonValue and horizon are missing
    } else {
      x2 <- horizon[[name]]
      if (x2 <= x1) stop(paste0('Horizon of ',name,' must exceed ',x1))
    } 
    prefix <- paste0(name,prettySpace[[name]],':')
    if (flat) {
      constFunc <- function(val) { 
        force(val)  # This ensures that the current value of y2 is stored (rather than the final value)
        function(yearvec) { ret <- rep(NA,length(yearvec)); ret[yearvec>x1] <- val; ret }
      }
      private$.par[[name]] <- constFunc(y2) 
      cat(prefix,'Flat extrapolation at',y1,'\n')
    } else {
      cubicConstFunc <- function(leftval,leftslope,right,rightval) {
        force(leftval)
        force(leftslope)
        force(right)
        force(rightval)
        splinefunc <- fitCubicSpline(x1,leftval,leftslope,right,rightval,0) 
        function(yearvec) {
          ret <- rep(NA,length(yearvec))  # Returns NA where function is not defined
          ret[yearvec > x1] <- splinefunc(yearvec[yearvec > x1])
          ret[yearvec > right] <- rightval
          ret
        }
      }
      private$.par[[name]] <- cubicConstFunc(y1,slope1,x2,y2)
      cat(prefix,'Cubic spline starting at',y1,'with slope',slope1,'and ending in year',x2,'at',y2,'\n')      
    }
  }
  private$.estimated  <- TRUE
  invisible(self)
})