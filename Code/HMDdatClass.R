
library("R6")

###
###  HMDdatClass
###

HMDdatClass <- R6Class("HMDdatClass", private = list(.data = NULL), public = list())

# -----------------------------------------------------------------------------
#           Accessor functions to private fields
# -----------------------------------------------------------------------------

HMDdatClass$set("active", "data" ,  function(value) { if (missing(value)) { private$.data } else { stop("data is read-only.\n" , call.=FALSE) } } )

# -----------------------------------------------------------------------------
#           Constructor
# -----------------------------------------------------------------------------

HMDdatClass$set("public", "initialize", function(Otxtfile, Etxtfile, groups = c("Female", "Male")) {
    groups <- match.arg(groups, choices = c("Female", "Male", "Total"), several.ok = T)
    
    if (!file.exists(Otxtfile)) stop("The 'Otxtfile' ", Otxtfile, " does not exist.")
    if (!file.exists(Etxtfile)) stop("The 'Etxtfile' ", Etxtfile, " does not exist.")
    if (!grepl("\\.txt$", Otxtfile)) stop("The 'Otxtfile' argument should be a '.txt'-file.")
    if (!grepl("\\.txt$", Etxtfile)) stop("The 'Etxtfile' argument should be a '.txt'-file.")
    #Changed: added stringAsFactors = F argument
    Ofirstline <- paste(read.table(file = Otxtfile,header = F,nrows = 1, stringsAsFactors = F), collapse = " ")
    Efirstline <- paste(read.table(file = Etxtfile,header = F,nrows = 1), collapse = " ")

    if (!grepl("Deaths"  , Ofirstline, fixed = TRUE)) stop("The 'Otxtfile' argument should be a '.txt'-file from HMD containing deaths.")
    if (!grepl("Exposure", Efirstline, fixed = TRUE)) stop("The 'Etxtfile' argument should be a '.txt'-file from HMD containing exposures.")
    
    if (gsub(",.*$", "", Ofirstline) == gsub(",.*$", "", Efirstline)) private$.data$country <- gsub(",.*$", "", Ofirstline)
    else stop("\n Countries in the two files do not match.\n Country in 'Otxtfile': ", gsub(",.*$", "", Ofirstline), ".\n Country in 'Etxtfile': ", gsub(",.*$", "", Ofirstline))
    age2int <- function(x) as.integer(gsub("[-+]", "", unlist(lapply(strsplit(as.character(x), split = "-"), "[[", 1))))
    DFO     <- read.table(file = Otxtfile, header = TRUE, skip = 2, na.strings = ".", as.is = TRUE); DFO$Age <- age2int(DFO$Age)
    DFE     <- read.table(file = Etxtfile, header = TRUE, skip = 2, na.strings = ".", as.is = TRUE); DFE$Age <- age2int(DFE$Age)
    
    if (!all(range(DFO$Age) == range(DFE$Age))) stop("Age-range in first and second argument do not match. First argument has ", min(DFO$Age),"-",max(DFO$Age), " while second argument has ",  min(DFE$Age),"-",max(DFE$Age),".")
    if (!all(range(DFO$Year) == range(DFE$Year))) stop("Year-range in first and second argument do not match. First argument has ", min(DFO$Year),"-",max(DFO$Year), " while second argument has ",  min(DFE$Year),"-",max(DFE$Year),".")
    if (!(nrow(DFO)==nrow(DFE))) stop("Number of rows extracted from datasets do not match. First argument has ", nrow(DFO), " rows, while second argument has ", nrow(DFE), " rows.")
    
    agelim  <- range(DFO$Age)
    timelim <- range(DFO$Year)
    
    private$.data$groups  <- groups
    private$.data$agelim  <- agelim
    private$.data$timelim <- timelim
    private$.data$O       <- array(NaN, dim=c(length(groups),lim2length(agelim),lim2length(timelim)), dimnames=list("Group"=groups, "Age"=lim2vec(agelim), "Year"=lim2vec(timelim)))
    private$.data$E       <- private$.data$O
    
    gidx <- 1L
    for (group in groups) {
        yidx <- 1L
        for (year in lim2vec(timelim)) {
            subO <- DFO[[group]][DFO$Year == year]
            subE <- DFE[[group]][DFE$Year == year]
            if (length(subO) != lim2length(agelim)) stop('Wrong number of ages for group "',group,'" in year ',year,' for occurence data.')
            if (length(subE) != lim2length(agelim)) stop('Wrong number of ages for group "',group,'" in year ',year,' for exposure data.')
            private$.data$O[gidx,,yidx] <- subO
            private$.data$E[gidx,,yidx] <- subE
            yidx <- yidx + 1L
        }
        gidx <- gidx + 1L
    }
    
    cat('Imported historic data for',private$.data$country,'\n')
    cat('  - Groups     :', length(groups),'\n')
    cat('  - Age range  :', agelim[1L],'-',agelim[2L],"\n")
    cat('  - Year range :', timelim[1L],'-',timelim[2L],"\n")
    invisible(self)
})

# -----------------------------------------------------------------------------
#           Method: Print method
# -----------------------------------------------------------------------------

HMDdatClass$set("public", "print", function(...) {
    cat("HMD-data object for ", private$.data$country, " containing data for ", length(private$.data$groups)," group(s).\n",sep='')
    cat("Groups:", paste(private$.data$groups, collapse = ", "),'\n')
    cat("Historic data loaded for ages ",private$.data$agelim[1],"-",private$.data$agelim[2]," and years ",private$.data$timelim[1],"-",private$.data$timelim[2],".\n",sep='') 
    invisible(self)
})

# -----------------------------------------------------------------------------
#           Methods for argument checking
# -----------------------------------------------------------------------------

HMDdatClass$set("private", ".verifyGroupName", function(groups) {
    if (is.null(groups)) {groups <- private$.data$groups}
    else {
        if (!is.character(groups)) stop("Argument 'groups' must be of type character.")
        groups <- properCase(groups)
        if (!all(groups %in% private$.data$groups)) stop("Argument 'groups' must be one (or more) of the following:\n    ",paste(private$.data$groups,collapse=", "))
    }
    private$.data$groups[private$.data$groups %in% groups] # Flip groups so they always appear in the order they were loaded
})

HMDdatClass$set("private", ".verifyAgelim", function(agelim) {
    if (is.null(agelim)) {agelim <- private$.data$agelim}
    else { 
        if (!is.numeric(agelim)) stop("Argument 'agelim' must be of type numeric.")
        if (length(agelim) > 2) stop("Argument 'agelim' must be of length 1 or 2.")
        agelim  <- lim2lim(agelim)
        if(agelim[1] < private$.data$agelim[1] | agelim[2] > private$.data$agelim[2])
            stop(paste(c("Ages must be in the interval ", private$.data$agelim[1]," to ", private$.data$agelim[2])))
    }
    agelim
})

HMDdatClass$set("private", ".verifyTimelim", function(timelim) {
    if (is.null(timelim)) {timelim <- private$.data$timelim}
    else {
        if (!is.numeric(timelim)) stop("Argument 'timelim' must be of type numeric.")
        if (length(timelim) > 2) stop("Argument 'timelim' must be of length 1 or 2.")
        timelim <- lim2lim(timelim)
        if(timelim[1] < private$.data$timelim[1] | timelim[2] > private$.data$timelim[2]) 
            stop(paste(c("Years must be in the interval ", private$.data$timelim[1]," to ", private$.data$timelim[2])))
    }
    timelim
})

# -----------------------------------------------------------------------------
#           Method: Extracts OE data
# -----------------------------------------------------------------------------

HMDdatClass$set("public", "getOEdata", function(groups = NULL, agelim = NULL, timelim = NULL) {
    groups    <- private$.verifyGroupName(groups)
    agelim    <- private$.verifyAgelim(agelim)
    timelim   <- private$.verifyTimelim(timelim)
    ret       <- list("groups" = groups, "agelim" = agelim, "timelim" = timelim, "O" = NULL, "E" = NULL)
    gidx      <- private$.data$groups %in% groups
    aidx      <- lim2vec(private$.data$agelim) %in% lim2vec(agelim)
    tidx      <- lim2vec(private$.data$timelim) %in% lim2vec(timelim)
    ret$O     <- private$.data$O[gidx,aidx,tidx,drop=F]
    ret$E     <- private$.data$E[gidx,aidx,tidx,drop=F]
    structure(ret, class="OEdata")
})

# -----------------------------------------------------------------------------
#           Method: 'smoothMu' returns a smoothed, historic mu surface
# -----------------------------------------------------------------------------

HMDdatClass$set("public", "getSmoothMu", function(groups = NULL, agelim = NULL, timelim = NULL) {
    groups  <- private$.verifyGroupName(groups)
    agelim  <- private$.verifyAgelim(agelim)
    timelim <- private$.verifyTimelim(timelim)
    if ((private$.data$agelim[1] > 0) || (private$.data$agelim[2] < 99)) stop('Historic data should at least span the ages 0-99.')
    
    agevec     <- 0L:private$.data$agelim[2L] # The calculation is performed on the entire age-span regardless of input
    yearvec    <- lim2vec(timelim)
    old.idx    <- (agevec >= 80) & (agevec < 100)
    old.weight <- seq(0,0.95,by=0.05)  # 20 weights
    ext.ages   <- 100L:private$.data$agelim[2L] # lastage is at least 99 (guaranteed at initialization)
    empty.mu   <- array(NA, dim = c(length(groups), length(agevec), length(yearvec)), dimnames = list("Group"=groups, "Age"=agevec, "Year"=yearvec))
    ret        <- list("groups" = groups, "agelim" = agelim, "timelim" = timelim, "mu" = empty.mu)

    OEdat <- self$getOEdata(groups, range(agevec), timelim)

    ret$mu[,1L,] <- (OEdat$O/OEdat$E)[,1L,] # use empirical rate at age 0
    for (gidx in seq_len(length(groups))) {
        for (yidx in seq_len(length(yearvec))) {
            O     <- OEdat$O[gidx,,yidx]
            E     <- OEdat$E[gidx,,yidx]
            logmu <- log(O/E)
            logmu[1L] <- NA     # exclude age 0
            logmu[O == 0] <- NA # exclude ages with zero death count
            if (is.na(logmu[2L])) logmu[2L] <- log(1/E[2L]) # Make sure age 1 is not missing (since age 0 is always missing, this would cause problems in loess below)
            
            # Ages 1-79
            res <- loess('logmu ~ age',data.frame(logmu=logmu, age=agevec, E=E),weights=E,span=0.25,na.action=na.omit)
            ret$mu[gidx,2L:80L,yidx] <- exp(predict(res,1:79))  # use smoothed values from age 1 to age 79
          
            # Ages 80-99
            logit.idx <- old.idx[(E[old.idx]>0) & (O[old.idx]>0)]  
            q         <- 1-exp(-O[logit.idx]/E[logit.idx])  # death probability
            logit.dat <- data.frame(logitq = logit(q), age = agevec[logit.idx], E = E[logit.idx])
            res.old <- lm('logitq ~ age', logit.dat, weights=E)
            q.pred  <- expit(predict(res.old,data.frame(age=80:99)))
            ret$mu[gidx,81L:100L,yidx] <- exp(predict(res,80:99))*(1-old.weight) - log(1-q.pred)*old.weight # convert q to mu
            
            # Ages 99+
            if (private$.data$agelim[2L] > 99L) {
                ret$mu[gidx,ext.ages+1,yidx] <- -log(1-expit(predict(res.old,data.frame(age=ext.ages))))
            }
            
            # Sometimes there are missing values around age 100 (due to missing values in logmu)
            if (any(is.na(ret$mu[gidx,,yidx]))) {
                na.idx   <- (1:length(agevec))[is.na(ret$mu[gidx,,yidx])]
                prev.idx <- (min(na.idx)-4):(min(na.idx)-1)
                post.idx <- (max(na.idx)+1):(max(na.idx)+4)
                post.idx <- post.idx[post.idx <= length(agevec)]
                fill.dat <- data.frame(logmu = log(ret$mu[gidx,c(prev.idx,post.idx),yidx]), idx = c(prev.idx,post.idx))
                ret.na   <- lm('logmu ~ idx',fill.dat)
                ret$mu[gidx,na.idx,yidx] <- exp(predict(ret.na, data.frame(idx=na.idx)))
            }
        }
    }
    
    ret$mu <- ret$mu[,lim2vec(private$.data$agelim) %in% lim2vec(agelim),,drop=F] # Adjust output to chosen agelim
    structure(ret, class = "surf")
})

