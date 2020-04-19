##############################################################################
##   Simulation: robustness of the analytical results under lognormal SAD   ##
##############################################################################
library(mobsim)
library(abind)
source('src/fns.R')
# spatial extent: 500x1000
cellsize=25 # cellsize
nIND=20000 # total abundance
spp_levels=c(100)
sad_levels=c(.5, 1.5, 3)
agg_levels=c(100,20)
smp_levels=c(25,50,100,200,400,800)
N_sim=100 

# function to divide census data into #spp x #quadrats abundance matrix
divide=function(comm, cellsize=25){
    plotdata=comm$census
    plotdata$x[plotdata$x == max(plotdata$x)] = 
        plotdata$x[plotdata$x == max(plotdata$x)] - 1e-05
    plotdata$y[plotdata$y == max(plotdata$y)] = 
        plotdata$y[plotdata$y == max(plotdata$y)] - 1e-05
    codex = plotdata$x%/%cellsize + 1
    codey = plotdata$y%/%cellsize + 1
    code = codex + codey * 1000
    sp = plotdata$sp
    # abundances of each species in each quadrat
    comp = tapply(rep(1, length(code)), list(sp,code), sum)
    comp[is.na(comp)] = 0
    
    return(comp)
}
betaDev_smp=array(NA,c(length(sad_levels),
    length(agg_levels),
    length(smp_levels),
    N_sim),
    dimnames = list(cv_sad=sad_levels,
        sigma_agg=agg_levels,
        N_smp=smp_levels,
        simID=1:N_sim)
)
beta.sim= numeric(1000)


sample_comm=function(comp,n=10){
    #@comp: species * sites matrix/dataframe
    res=comp[,sample(ncol(comp),n)]
    res=res[rowSums(res)>0,]
    return(res)
}


set.seed(123)
for(i_sp in 1:length(spp_levels)){# only 1 sp level here
    cat('sp_level',i_sp,':\n')
    for(i_ab in 1:length(sad_levels)){# 3 sad levels
        cat('ab_level',i_ab,':\n')
        for(i_ag in 1:length(agg_levels)){# 2 agg levels
            cat('ag_level',i_ag,':\n')
            dat=sim_thomas_community(s_pool = spp_levels[i_sp],
                                     n_sim = nIND,
                                     sad_coef = list(cv_abund = sad_levels[i_ab]),
                                     fix_s_sim = T,
                                     xrange = c(0,500),
                                     yrange = c(0,1000),
                                     sigma = agg_levels[i_ag])
            comp=divide(dat,cellsize)
            
            for(i in 1:length(smp_levels)){# 6 smp levels
                cat('smp_level',i,':\n')
                if(smp_levels[i]==ncol(comp)){
                    ind=data.frame(sp=rep(rownames(comp),rowSums(comp)), 
                                       subplot=unlist(apply(comp,1,function(x)rep(colnames(comp),x))))
                    b=Beta(comp)
                    for(k in 1:1000){#this 1k is for Kraft's Null
                        comp.sim=commsimu(ind)
                        beta.sim[k]=Beta(comp.sim)
                    }
                    dev=(b-mean(beta.sim))/sd(beta.sim)
                    for(j in 1:N_sim)
                        betaDev_smp[i_ab, i_ag, i, j]=dev
                }
                else{
                    for(j in 1:N_sim){
                        comp.smp=sample_comm(comp,smp_levels[i]);
                        ind=data.frame(sp=rep(rownames(comp.smp),rowSums(comp.smp)), 
                                           subplot=unlist(apply(comp.smp,1,function(x)rep(colnames(comp.smp),x))))
                        b=Beta(comp.smp)
                        for(k in 1:1000){
                            comp.sim=commsimu(ind)
                            beta.sim[k]=Beta(comp.sim)
                        }
                        dev=(b-mean(beta.sim))/sd(beta.sim)
                        betaDev_smp[i_ab, i_ag, i, j] = dev
                        
                        
                        if(j%%20==1)
                            cat(j,'\t')
                    }
                }
                
                cat('\n')
            }
        }
        #save(betaDev_smp,file='res/logNorm_sim.RData')
    }
}


############ ploting ------------
# load('res/logNorm_sim.RData')
# smp_levels=c(25,50,100,200,400,800)
# for plotting purpose (on a log-log scale), ignore the few non-positive beta_dev (10 out of 3600).
betaDev_smp[betaDev_smp<=0]=NA
win.metafile(paste0('figs/smpl_logNorm_spp100x3sadx2agg.wmf'),width=7,height=5)
par(mfcol=c(2,3),mar=c(1.1,1.1,1.1,1.1),oma=c(2.5,2.5,1.5,1.5),mgp=c(2.5,.5,0))

for(i in 1:3){
    for(j in 1:2){
        boxplot(t(betaDev_smp[i,j,,]),at=smp_levels,log='xy',name=smp_levels,boxwex=.2)
        curve(sqrt(x/smp_levels[length(smp_levels)])*betaDev_smp[i,j,length(smp_levels),1],
              add=T,col='grey',lwd=2)
        
        op=par(xpd=NA)
        if(j==1&i==3)mtext(expression(paste('Low spatial aggregation, sigma = 100')),4,line=.8,cex=.7)
        if(j==2&i==3)mtext(expression(paste('High spatial aggregation, sigma = 20')),4,line=.8,cex=.7)
        
        if(j==1&i==1)mtext('cv_abund = 0.5, Even',3,cex=.7,line=.8)
        if(j==1&i==2)mtext('cv_abund = 1.5',3,cex=.7,line=.8)
        if(j==1&i==3)mtext('cv_abund = 3.0, Uneven',3,cex=.7,line=.8)
        
        par(op)
    }
}
mtext('Sampling effort (Number of quadrats)',side = 1,line=.8,outer = T)
mtext(expression(paste(beta,'-deviation')), side = 2,line=.8, outer=T)
dev.off()

