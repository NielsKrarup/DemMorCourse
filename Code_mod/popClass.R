
library("R6")

###
###  popClass
###

popClass <- R6Class("popClass", private = list(  # private variables for internal logic
                                                .groups      = NULL,
                                                .country     = NULL,
                                                .mort        = NULL,
                                                .fert        = NULL,
                                                .groupapply  = NULL,   # Set by the projetion method
                                                .groupprop   = NULL),  # Set by the projetion method
                                public = list( # public variables for data 
                                                histGroups   = NULL,
                                                histAgelim   = NULL,   # Vector of the form (minage,maxage) corresponding to historic data 
                                                histTimelim  = NULL,   # Vector of the form (minyear,maxyear) corresponding to historic data 
                                                histN        = NULL))  # list of population surfaces, one for each group


# -----------------------------------------------------------------------------
#           Accessor functions to private fields
# -----------------------------------------------------------------------------

# Accessor functions to private fields
popClass$set("active", "country" ,  function(value) { if (missing(value)) { private$.country} else { stop("country is read-only.\n" , call.=FALSE) } } )
popClass$set("active", "groups" ,  function(value) { if (missing(value)) { private$.groups } else { stop("groups is read-only.\n" , call.=FALSE) } } )
popClass$set("active", "mort" ,  function(value) { if (missing(value)) { private$.mort } else { stop("groups is read-only.\n" , call.=FALSE) } } )


# -----------------------------------------------------------------------------
#           Constructor
# -----------------------------------------------------------------------------

# Overrides default behaviour of $new()
popClass$set("public", "initialize", function(popTxtFile, groups = c("Female", "Male"), maxage = 100, timelim = NULL) {
    groups <- match.arg(groups, choices = c("Female", "Male", "Total"), several.ok = T)
    private$.groups <- groups
    
    if (!file.exists(popTxtFile)) stop("First argument should be a filename.")
    if (!grepl("\\.txt$", popTxtFile)) stop("First argument should be a '.txt'-file.")
    
    firstline <- paste(read.table(file = popTxtFile,header = F,nrows = 1), collapse = " ")
    
    
    if (!grepl("Population", firstline, fixed = TRUE)) stop("First argument should be a '.txt'-file from HMD containing population counts.")
    private$.country <- gsub(",.*$", "", firstline)
    
    age2int <- function(x) as.integer(gsub("[-+]", "", unlist(lapply(strsplit(as.character(x), split = "-"), "[[", 1))))
    DF      <- read.table(file = popTxtFile, header = TRUE, skip = 2, na.strings = ".", as.is = TRUE); DF$Age  <- age2int(DF$Age)
    DF      <- DF[DF$Age <= maxage,]
    agelim  <- range(DF$Age)
    if ((agelim[1] > 0) || (agelim[2] < maxage)) stop('Historic data should at least span the ages 0-',maxage,'.\n')
    
    if (!is.null(timelim)) {
        timelim <- lim2lim(timelim)
        if (!all(lim2vec(timelim) %in% unique(DF$Year))) stop("Historic data spans the years ", min(DF$Year), "-", max(DF$Year),".\n")
        DF <- DF[DF$Year %in% lim2vec(timelim),]
    }
    else {timelim <- range(DF$Year)}
    
    N    <- array(NaN, dim=c(length(groups),lim2length(agelim),lim2length(timelim)), dimnames=list(groups, lim2vec(agelim), lim2vec(timelim)))
    gidx <- 1L
    for (group in groups) {
        yidx <- 1L
        for (year in lim2vec(timelim)) {
            subdat <- DF[[group]][DF$Year == year]
            if (length(subdat) != lim2length(agelim)) stop('Wrong number of ages for group "',group,'" in year ',year,'.')
            N[gidx,,yidx] <- subdat
            yidx <- yidx + 1L
        }
        histN <- N[gidx,,]
        if (any(is.na(histN)))  cat("There are",sum(is.na(histN)),"NA's in population data for group",groups[gidx],"\n")
        if (any(is.nan(histN))) cat("There are",sum(is.nan(histN)),"NaN's in population data for group",groups[gidx],"\n")    
        gidx <- gidx + 1L
    }
    
    self$histGroups  <- groups
    self$histAgelim  <- agelim
    self$histTimelim <- timelim
    self$histN       <- N

    cat('Imported historic data for',private$.country,'\n')
    cat('  - Groups     :', length(groups),'\n')
    cat('  - Age range  :', agelim[1L],'-',agelim[2L],"\n")
    cat('  - Year range :', timelim[1L],'-',timelim[2L],"\n")
    invisible(self)
})

popClass$set("public", "setMortAndFert", function(surf, surf.fer) {
    if (class(surf)[1] != "surf") stop('First argument must be of class surf')
    if (!all(lim2vec(self$histAgelim) %in% lim2vec(surf$agelim))) stop("surf object must span the ages ", self$histAgelim[1],"-", self$histAgelim[2],".\n")
    if (!all(names(private$.groups) %in% surf$groups)) stop('Mortality object must contain the group(s): ',paste(private$.groups,collapse=', '))
    if (class(surf.fer)[1] != "surf.fer") stop('Second argument must be of class surf.fer')
    timelim <- c(min(surf.fer$timelim[1], surf$timelim[1]),min(surf.fer$timelim[2], surf$timelim[2]))
    mu  <- surf2mu(surf, groups = self$histGroups, agelim = self$histAgelim, timelim = timelim)
    fer <- surf2fer(surf.fer, timelim = timelim)
    private$.mort <- mu2surf(self$histGroups, self$histAgelim, timelim, mu)
    private$.fert <- fer2surf(surf.fer$agelim, timelim, fer)
    invisible(self)
})


# -----------------------------------------------------------------------------
#         Method definitions
# -----------------------------------------------------------------------------

popClass$set("public", "print", function(...) {
    cat("Population object for ", private$.country, " containing data for ", length(private$.groups)," group(s).\n",sep='')
    cat("Groups:", paste(private$.groups, collapse = ", "),'\n')
    cat("Historic data loaded for ages ",self$histAgelim[1],"-",self$histAgelim[2]," and years ",self$histTimelim[1],"-",self$histTimelim[2],".\n",sep='') 
    if (!is.null(private$.mort)) {
        cat("\nMortality data loaded for ages", private$.mort$agelim[1],"-",private$.mort$agelim[2], "and years", private$.mort$timelim[1],"-",private$.mort$timelim[2])
        cat("\n Fertiliy data loaded for ages", private$.fert$agelim[1],"-",private$.fert$agelim[2], "and years", private$.fert$timelim[1],"-",private$.fert$timelim[2])
    }
    else{
        cat("Mortality and fertility data not loaded.")
    }
    invisible(self)
})

# Returns various statistics in a data frame with the format: Group Age Year Stat1 ... StatN
popClass$set("public", "plotdat", function(groups = NULL, agelim = NULL, timelim = NULL, groupapply = "Female", groupproportion = 0.487, stableyear=NULL) {
    if (is.null(agelim)) agelim  <- self$histAgelim
    else {
        agelim  <- lim2lim(agelim)
        if(agelim[1] < self$histAgelim[1] | agelim[2] > self$histAgelim[2]) stop(paste(c("Ages must be in the interval ", self$histAgelim[1]," to ", self$histAgelim[2])))
    }
    pop <- self$getPop(groups, self$histAgelim, timelim, groupapply, groupproportion, stableyear) # argument checks performed here

    groups  <- pop$groups
    agevec  <- lim2vec(agelim)
    yearvec <- lim2vec(pop$timelim)
    aidx    <- lim2vec(self$histAgelim) %in% agevec
    
    composition <- apply(pop$N,c(1,3), function(x) x/sum(x))
    composition <- aperm(composition, c(2,1,3))

    data.frame( Group  = rep(groups, times = length(agevec)*length(yearvec))
              , Age    = rep(agevec, each = length(groups), times = length(yearvec))
              , Year   = rep(yearvec, each = length(groups)*length(agevec))
              , D      = as.vector(pop$D[,aidx,])
              , B      = as.vector(pop$B[,aidx,])
              , N      = as.vector(pop$N[,aidx,])
              , c      = as.vector(composition[,aidx,])
              , stable = as.vector(pop$stable[,aidx,])
    )
})
 
popClass$set("public", "getPop", function(groups = NULL, agelim = NULL, timelim = NULL, groupapply = "Female", groupproportion = 0.487, stableyear=NULL) {
    if (is.null(groups)) groups <- self$histGroups
    else {
        if (!all(groups %in% self$histGroups)) stop("The following groups are not contained within 'surf': ", paste(groups[!(groups %in% self$histGroups)], sep = ", "))
        groups <- self$histGroups[self$histGroups %in% groups]
    }
    groups <- self$histGroups[self$histGroups %in% groups] # Flip groups so they always appear in the order they were loaded
    if (is.null(agelim)) agelim  <- self$histAgelim
    else {
        agelim  <- lim2lim(agelim)
        if(agelim[1] < self$histAgelim[1] | agelim[2] > self$histAgelim[2]) stop(paste(c("Ages must be in the interval ", self$histAgelim[1]," to ", self$histAgelim[2])))
    }
    if (is.null(timelim)) timelim <- self$histTimelim
    else {
        timelim <- lim2lim(timelim)
        targetLim <- self$histTimelim
        if (timelim[2] > targetLim[2] & (is.null(private$.mort) | is.null(private$.fert))) stop("You are trying to project the population beyond the last year of historical data but mort and fert objects not initialized. Execute method $setMortAndFert.\n")
        if (!is.null(private$.mort) & !is.null(private$.fert)) targetLim[2L] <- private$.mort$timelim[2L]+1L
        if(timelim[1] < targetLim[1] | timelim[2] > targetLim[2]) stop(paste(c("Years must be in the interval ", targetLim[1]," to ", targetLim[2])))
    }
    
    empty.array <- array(NaN, dim=c(length(groups),lim2length(agelim),lim2length(timelim)), dimnames=list("Group"=groups, "Age"=lim2vec(agelim), "Year"=lim2vec(timelim)))
    empty.mat   <- array(NaN, dim=c(length(groups),lim2length(timelim)), dimnames=list("Group"=groups, "Year"=lim2vec(timelim)))
    ret <- structure(list(groups = groups, agelim=agelim, timelim=timelim, D = empty.array, B = empty.array, N = empty.array, G = empty.mat, stable = empty.array), class='popdat')
    
    gidx <- self$histGroups %in% groups
    aidx <- lim2vec(self$histAgelim) %in% lim2vec(agelim)
    
    t0 <- 1L
    # Do we need historic data?
    if (timelim[1L]<=self$histTimelim[2L]) {
        tidx <- lim2vec(self$histTimelim) %in% lim2vec(timelim)
        ret$N[,,1:sum(tidx)] <- self$histN[gidx,aidx,tidx]
        t0 <- sum(tidx) # .. Do not +1, overwrite below
    }
    
    # Do we need to project?
    horizon <- timelim[2L] - self$histTimelim[2L] + 1L
    if (horizon > 0) {
        pop.for <- self$project(horizon, groupapply, groupproportion, stableyear)
        tidx <- lim2vec(pop.for$timelim) %in% lim2vec(timelim)
        ret$D[,,t0:lim2length(timelim)] <- pop.for$D[gidx,aidx,tidx]
        ret$B[,,t0:lim2length(timelim)] <- pop.for$B[gidx,aidx,tidx]
        ret$N[,,t0:lim2length(timelim)] <- pop.for$N[gidx,aidx,tidx]
        ret$G[,t0:lim2length(timelim)]  <- pop.for$G[gidx,tidx]
        ret$stable[,,t0:lim2length(timelim)] <- pop.for$stable[gidx,aidx,tidx]
    }

    ret
})


# Median projection based on given mortality and fertility object.
# The groups of the mortality object must contain  the ones in the population.
#
# There is only one fertility surface, i.e. no groups, but it can be interpreted in two different ways:
# 1) The fertility object is used on each group individually with the children being born belonging
#    to that group only. This is achieved by setting groupapply = "Groupwise" (and groupproportion=NA).
#    The idea is that the population is a total population of males and females and that the groups,
#    if more than one, correpond to a division along another dimension than sex, e.g. social groups
# 2) The fertility object is applied to one of two groups and the births are distributed between the two groups.
#    Groupportion determines the fraction of births of the 'reproducing' group, the rest belongs to the other group.
#    This is the functionality when:
#    - there are two groups in the population
#    - groupapply is the name of one of the groups
#    - groupproportion is set to a number between 0 and 1
#
# It is the responsibility of the user to scale the fertility object appropriately, i.e. brutto fertility
# in case 2) should be roughly twice as large as brutto fertility in case 1), since in case 1)
# the 'fertility' surface applies to both males and females, while in case 2) it only applies only to females.
#
# If stableyear is set, mortality and fertility rates for that year is used in the projection.
# To illustrate the convergence to different stable populations.
#
# Warning: The code performs no tests on fertility to judge whether the chosen functionality is the one intended.

popClass$set("public", "project", function(horizon, groupapply = NULL, groupproportion = NA, stableyear=NULL) {
    if (is.null(groupapply)) {
        usagestr <- "Groupapply must be set to either 'Groupwise' or the name of the group to which it applies, e.g. 'Female'.\n"
        usagestr <- paste0(usagestr,"In the latter case, there must be two groups in the population and the births are distributed between them.\n")
        usagestr <- paste0(usagestr,"Groupproportion determines the fraction of births of the 'reproducing' gender,\n")
        usagestr <- paste0(usagestr,"the remaining births are assumed to belong to the other group, e.g. 'Male'.\n")
        usagestr <- paste0(usagestr,"This functionality is included to support the usual definition of fertility.\n")
        cat(usagestr)
        return(NULL)
    }
    stopifnot(horizon>0)
    
    # Check of groupapply and groupproprotion are specified correctly
    legalnames <- c(private$.groups,'Groupwise')  # 'Groupwise' is reserved and cannot be used as a group name
    if (!any(groupapply %in% legalnames)) stop('Groupapply must be one of the following: ',paste(legalnames,collapse=' '))
    if (groupapply == 'Groupwise') {
    if (!is.na(groupproportion)) stop("Groupproportion must not be specified when groupapply='Groupwise'")
        cat('Groupwise mode: Fertility will be applied to each group separately (unisex fertility).\n')
        fertgrpidx    <- 0  # indicates that fert applies to all groups
        nonfertgrpidx <- 0  # indicates that fert applies to all groups
    } else {
        if (length(private$.groups) != 2) stop("The population must have two groups in two-group mode.")
        if (is.na(groupproportion)) stop("Groupproportion must be specified in two-group mode.")
        if ((groupproportion < 0) || (groupproportion>1)) stop("Groupproportion must be between 0 and 1.")
        fertgrpidx    <- which(groupapply == private$.groups)
        nonfertgrpidx <- 3-fertgrpidx
        cat('Two-group mode: Fertility will be applied to group',groupapply,'only (one-sex fertility).\n')
        cat('                Births are distributed with',groupproportion,'to', private$.groups[fertgrpidx],'and',1-groupproportion,'to', private$.groups[nonfertgrpidx],'\n')
    }
    if (is.null(private$.mort) | is.null(private$.fert)) stop("Mort and fert objects not initialized. Execute method $setMortAndFert.\n")
    
    # fertvec contains age-specific fertility rates from age 0 to length(fertvec)-1
    # survvec contains age-specific survival probabilities for the same ages (or longer)
    # The function returns the stable population growth, i.e. the solution r
    # to the eqaution sum_{age=0} exp(-ra)fertility(age)survivial(age) = 1 
    # Confer eq. (7.6) p.144 in Demography by Preston et al.
    calc.stable.growth.rate <- function(fertvec, survvec) {
        agevec <- 1:length(fertvec)
        wvec   <- fertvec*survvec[agevec]
        func  <- function(r) { 1-sum(exp(-r*agevec)*wvec) }
        uniroot(func,c(-0.01,0.01),extendInt="yes", tol = 1e-8)$root
    }
    
    groups     <- private$.groups
    agelim     <- self$histAgelim
    timelim    <- self$histTimelim[2L] + c(0,horizon-1)
    agelimFer  <- c(0, private$.fert$agelim[2])
    
    # check if mort and fert contain the required years
    if (is.null(stableyear)) {
        if (timelim[2] > (private$.mort$timelim[2]+1) | timelim[2] > (private$.fert$timelim[2]+1)) stop("Horizon exceeds available mortality and fertility data. End year for forecast is ", timelim[2], " but mort and fert only have data until ", private$.mort$timelim[2], " and ", private$.fert$timelim[2])
    }
    else {
        if(private$.mort$timelim[1] > stableyear | stableyear > private$.mort$timelim[2]) stop("Stable year exceeds available mortality data. Stable year is ", stableyear, " but there is only mortality data for ", private$.mort$timelim[1], " to ", private$.mort$timelim[2])
        if(private$.fert$timelim[1] > stableyear | stableyear > private$.fert$timelim[2]) stop("Stable year exceeds available fertility data. Stable year is ", stableyear, " but there is only fertility data for ", private$.fert$timelim[1], " to ", private$.fert$timelim[2])
    }
    
    m.surf <- if (is.null(stableyear)) surf2mu(private$.mort, groups, agelim, timelim) else surf2mu(private$.mort, groups, agelim, stableyear)
    f.surf <- if (is.null(stableyear)) surf2fer(private$.fert, agelimFer, timelim) else surf2fer(private$.fert, agelimFer, stableyear)
    
    empty.array <- array(NaN, dim=c(length(groups),lim2length(agelim),lim2length(timelim)), dimnames=list("Group"=groups, "Age"=lim2vec(agelim), "Year"=lim2vec(timelim)))
    N <- deaths <- births <- stable <- empty.array
    N[,,1] <- self$histN[,,lim2length(self$histTimelim)]
    G <- array(NaN, dim=c(length(groups), lim2length(timelim)), dimnames=list("Group"=groups, "Year"=lim2vec(timelim)))
    
    agevec      <- lim2vec(agelim)
    initpopaidx <- seq_len(lim2length(agelim)) # this is used to keep track of the initial population over time
    fertageidx  <- agevec %in% lim2vec(agelimFer) # Not every age is fertile
 
    # Order of group updates. For groupwise fertility, any order will do.
    # In the two-group mode, it is best to start with the fertility group and calculate all births in the first step.
    # Since then we can update N for both the fertility and non-fertility group in the same loop (instead of two loops).
    if (groupapply == 'Groupwise') { grpidxvec <- 1:length(private$.groups) } else { grpidxvec <- c(fertgrpidx,3-fertgrpidx) }
    for (tidx in seq_len(lim2length(timelim))) {
        tidx2 <- ifelse(is.null(stableyear), tidx, 1)
        for (gidx in grpidxvec) {
            primoN <- N[gidx,,tidx]
            
            mortvec <- c(m.surf[gidx,-lim2length(agelim),tidx2], Inf) # Add infinite mu for the last age
            fertvec <- f.surf[,tidx2]
            
            deaths[gidx,,tidx] <- primoN*(1-exp(-mortvec))
            births[gidx,,tidx] <- 0
            if ((groupapply == 'Groupwise') || (gidx == fertgrpidx)) {
                b <- rep(0,lim2length(agelim))
                b[fertageidx] <- primoN[fertageidx]*fertvec
                births[gidx,,tidx] <- b
            }
            
            grpprop <- 1.0
            # In two-group mode, only a fraction of births belong to the fertile group.
            if (gidx == fertgrpidx)    { grpprop <- groupproportion }
            if (gidx == nonfertgrpidx) { grpprop <- 1.0-groupproportion }
            
            # Calculate births attributable to the initial population
            G[gidx,tidx] <- grpprop * sum(births[gidx,initpopaidx,tidx])
            
            # Calculate the new primo N
            if (tidx < lim2length(timelim)) {
                survivors <- (primoN*exp(-mortvec))[-length(primoN)]  # the last age group dies completely
                totbirths <- sum(b)
                N[gidx,,tidx+1] <- c(grpprop*totbirths, survivors)
            }

            # Calculate stable equivalent population
            survprob <- c(1,exp(-cumsum(mortvec[-length(primoN)])))
            fertgrpprop <- 1.0
            
            # In two-group mode, groupportion should be used for both the fertile and non-fertile groups (not 1-groupprotion as one might initially suspect)
            if (groupapply != 'Groupwise') { fertgrpprop <- groupproportion }
            r <- calc.stable.growth.rate(fertgrpprop*fertvec, survprob)
            if (!is.null(stableyear) && (tidx == 1)) {
                fertagevec <- agevec[fertageidx]+1
                Astar      <- sum(fertagevec*exp(-r*fertagevec)*fertgrpprop*fertvec*survprob[1:length(fertvec)])
                wlifeexp   <- sum(exp(-r*(1:length(survprob)-1))*survprob)
                cat('Intrinsic growth rate for',groups[gidx],'     :',r,'\n')
                cat('Weighted mean age at birth for',groups[gidx],':',Astar,'\n')
                cat('Weighted life expectancy for',groups[gidx],'  :',wlifeexp,'\n')
                cat('Net reproduction rate for',groups[gidx],'     :',grpprop*sum(fertvec*survprob[1:length(fertvec)]),'\n')
            }
            wvec     <- exp(-r*agevec)*survprob
            stable[gidx,,tidx] <- wvec/sum(wvec)
        }
        initpopaidx <- initpopaidx[-1]
    }
    private$.groupapply <- groupapply
    private$.groupprop  <- groupproportion
    structure(list(groups = groups, agelim=agelim, timelim=timelim, D = deaths, B = births, N = N, G = G, stable = stable), class='popdat')
})
