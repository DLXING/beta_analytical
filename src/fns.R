############################################################################
## Functions for implementing and testing the analytical models of beta   ##
##  diversity of Xing & He (2019) "Analytical models for β-diversity and  ##
##  the power-law scaling of β-deviation"                                 ##
## Written by Dingliang Xing (XingDingliang@gmail.com)                    ##
## August 1st, 2019                                                       ##
############################################################################

############################################################################
## source the functions written in C++ for a fast implementation of Kraft ##
## et al's (2011) randomized null model.                                  ##
############################################################################
library(Rcpp)
Rcpp::sourceCpp('src/fns.cpp')

############################################################################
## function for finding the mete-derived logseries parameter (lambda)     ##
## based on state variables S and N (sensu Harte 2011)                    ##
############################################################################
meteSN=function(S0, N0, version='precise'){
    if(!length(S0) == length(N0)) stop("S and N must have the same length")
    if(!all(S0[] > 1,na.rm=T)) stop("S must be greater than 1")
    if(!all(N0[] > 0,na.rm=T)) stop("N must be greater than 0")
    if(!all(S0[]/N0[] < 1,na.rm=T)) stop("N must be greater than S")
    if(!version %in% c('precise', 'approx')) stop("version must be either 'precise' or 'approx'")
    
    res = rep(NA,length(N0))
    for(i in 1:length(S0)){
        S=S0[i]
        N=N0[i]
        SNid=paste(S,N)
        NSratio=N/S
        if(!is.na(S/N)){
            if(version == 'precise'){
                require(VGAM) # for the lerch() function 
                # eqn 7.29 of Harte (2011)
                y = function(x) S/N * x*(1-x^N)/(1-x) -
                    (-log(1-x)-x^(N+1)*lerch(x,1,N+1,iter = 1000))
                p = try(uniroot(y, lower=1e-15, upper=1-1e-7,tol=1e-15)$`root`)
                if(class(p)!='try-error')
                    res[i]= -log(p)
            }
            else if(version == 'approx'){
                require(pracma) # for the fsolve function
                # eqn 7.30 of Harte (2011)
                y = function(x) x * log(1/x) - S / N
                res[i]=fsolve(y,1/N)
            }
        }
    }
    return (res)
}

############################################################################
## function for calculating observed beta from a community composition    ##
## matrix                                                                 ##
############################################################################
beta_obs=function(comp){
    # @comp - S*M community matrix containing either occurence or abundance data
    comp[comp>0]=1 # convert it to occurence if comp is an abundance matrix
    if (any(rowSums(comp)==0)|any(colSums(comp)==0))
        stop('check comp: 0 row/col sums detected.')
    S=nrow(comp)
    M=ncol(comp)
    alphaMean=mean(colSums(comp))
    res=1-alphaMean/S
    return(res)
}

############################################################################
## function for calculating analytical null beta (eqn 2a in Xing & He)    ##
############################################################################
beta_null=function(M,S,N,version='precise'){
    # @M - num of sampling units (plots or quadrats)
    # @S - num of spp
    # @N - num of tot individuals
    # @version - the method used to find the logseries parameter
    lambda=meteSN(S,N,version=version)
    p=exp(-lambda)
    log(1-p+p/M)/log(1-p)
}
############################################################################
## function for calculating standard deviaiton of null beta               ##
## (eqn 2b in Xing & He)                                                  ##
############################################################################
beta_null_sd=function(M,S,N,version='precise'){
    # @M - num of sampling units (plots or quadrats)
    # @S - num of spp
    # @N - num of tot individuals
    # @version - the method used to find the logseries parameter
    lambda=meteSN(S,N,version=version)
    p=exp(-lambda)
    
    res=1/S/M/log(1-p)*((M-1)*log(1-p*(1-2/M)) +
                            log(1-p*(1-1/M)) -
                            M*log(1-p*(1-1/M)^2))
    sqrt(res)
}
############################################################################
## function for calculating analytical beta deviaiton                     ##
############################################################################
beta_dev=function(comp,N=NULL, version='precise'){
    # @comp - S*M community matrix containing either occurence or abundance data
    # @N - metacommunity size (total abundance), 
    #  is either specified (in cases comp is occurence data)
    #  or calculated from the community matrix (when comp is abundance data).
    # @version - the method used to find the logseries parameter
    if(is.null(N)){
        if(all(comp<=1))
            stop('N not available.')
        else N=sum(comp)
    }
    M=ncol(comp)
    S=nrow(comp)
    (beta_obs(comp) - beta_null(M,S,N,version)) / beta_null_sd(M,S,N,version)
}
############################################################################
## function for calculating standard deviaiton of beta deviation          ##
## (eqn 3 in Xing & He)                                                   ##
############################################################################
beta_dev_sd=function(M,S,N,version='precise'){
    # @M - num of sampling units (plots or quadrats)
    # @S - num of spp
    # @N - num of tot individuals
    # @version - the method used to find the logseries parameter
    lambda=meteSN(S,N,version=version)
    p=exp(-lambda)
    v_Phi=1/S/log(1-p)* (log(1-p*(1-1/M)^2) -
                             1/log(1-p) * log(1-p*(1-1/M))^2 )
    v_Pi=1/S/M/log(1-p)* ((M-1)*log(1-p*(1-2/M)) +
                              log(1-p*(1-1/M)) -
                              M*log(1-p*(1-1/M)^2))
    sqrt(v_Phi/v_Pi)
}

############################################################################
## function for calculating analytical NBD beta (eqn 4a in Xing & He)     ##
############################################################################
beta_NBD=function(M,S,N,version='precise',kappa=1){
    # @M - num of sampling units (plots or quadrats)
    # @S - num of spp
    # @N - num of tot individuals
    # @version - the method used to find the logseries parameter
    # @kappa - the aggregation parameter
    lambda=meteSN(S,N,version=version)
    p=exp(-lambda)
    n=1:max(N,1e5)
    -1/log(1-p)*sum(p^n/n*(1+n/M/kappa)^-kappa)
}
############################################################################
## function for calculating standard deviaiton of NBD beta                ##
## (eqn 4b in Xing & He)                                                  ##
############################################################################
beta_NBD_sd=function(M,S,N,version='precise',kappa=1,level=2){
    # @M - num of sampling units (plots or quadrats)
    # @S - num of spp
    # @N - num of tot individuals
    # @version - the method used to find the logseries parameter
    # @kappa - the aggregation parameter
    # @level - indicates whether only the variation due to spatial pattern 
    #  (level = 1) or both the variation due to spatial pattern and SAD 
    #  (level = 2, the default) will be calculated.
    lambda=meteSN(S,N,version=version)
    p=exp(-lambda)
    
    if(level==1)
        res=1/S/M/log(1-p)*(log(1-p*exp(-1/M)) -
                                log(1-p*exp(-2/M)) +
                                p*exp(-2/M)/M/(1-p*exp(-2/M)))
    else if(level==2){
        n=1:max(N,1e5)
        v_Phi=-1/log(1-p)/S* (sum(p^n/n*(1+n/M/kappa)^(-2*kappa)) +
                                  1/log(1-p)*sum(p^n/n*(1+n/M/kappa)^-kappa)^2)
        res=1/S/M/log(1-p)*(log(1-p*exp(-1/M)) -
                                log(1-p*exp(-2/M)) +
                                p*exp(-2/M)/M/(1-p*exp(-2/M))) + v_Phi
    }
    else stop('level can only take value of 1 or 2!')
    
    sqrt(res)
}

############################################################################
## function for estimating the aggregation parameter kappa (k)            ##
############################################################################
kappa_est=function(comp){
    # @comp - S*M community matrix containing abundance data
    M=ncol(comp)
    n=rowSums(comp)
    comp1=comp
    comp1[comp1>0]=1
    o=rowSums(comp1)
    negll=function(kappa,o,n)#negative log likelihood
        -sum(o*log(1-(1+n/M/kappa)^(-kappa))-
                 (M-o)*kappa*log(1+n/M/kappa))
    optimize(negll,c(0,1e5),o=o,n=n)[[1]]
}

############################################################################
## KS (Kolmogorov-Smirnov) test of empirical SADs vs. logseries           ##
############################################################################
library(dgof)
library(extraDistr)
sad_ks_test=function(comp,version='precise'){
    # @comp - S*M community matrix containing abundance data
    # @version - the method used to find the logseries parameter
    # the function returns a p-value of the KS test
    abund=rowSums(comp)
    N=sum(abund)
    S=length(abund)
    lambda=meteSN(S,N,version=version)
    if(!is.na(lambda)){
        p=exp(-lambda)
        x=1:max(abund)
        cdf_ls=stepfun(x,c(0,plgser(x,p)))
        dgof::ks.test(abund,cdf_ls)$p.value
    }
    else
        NA
}
############################################################################
## function for converting community matrix to individual-level data      ##
############################################################################
mat_to_ind=function(comp.dat){
    ## @comp.dat- list of community matrices with S*M abundance data
    for(i in 1:length(comp.dat)){
        mat=comp.dat[[i]]
        spnames=rownames(mat);abund=rowSums(mat)
        sp=rep(spnames,abund)
        plotnames=colnames(mat)
        subplot=NULL
        for(j in 1:nrow(mat))
            subplot=c(subplot,rep(plotnames,mat[j,]))
        comp.dat[[i]]=data.frame(sp=sp,subplot=subplot)
    }
    comp.dat
}
############################################################################
## function for R squred with respect to the 1:1 line                     ##
############################################################################
R2=function(obs,pred){
    na.id=c(which(is.na(obs)),which(is.na(pred)))
    if(length(na.id)>0){
        obs=obs[-na.id]
        pred=pred[-na.id]
    }
    1-sum((pred-obs)^2)/sum((obs-mean(obs))^2)
}
