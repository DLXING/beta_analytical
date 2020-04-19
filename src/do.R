source('src/fns.R')
## beta and beta_dev for Gentry's data--------
load('data/gentry198.RData')
## gentry198.RData stores two objects:
## comp.dat: list of the 198 community matrices (each S*M, and M=10 for all).
## sites: data.frame with coordinates (and more) of all communities. 
ind.dat=mat_to_ind(comp.dat)
res_simu_beta=do_simu(ind.dat,nsim=999)
res_Gentry=data.frame(S=sapply(comp.dat,nrow),
                      N=sapply(comp.dat,sum),
                      M=sapply(comp.dat,ncol),
                      lat=abs(sites$lat),
                      beta_obs=sapply(comp.dat,beta_obs),
                      beta_rand_null=apply(res_simu_beta,1,mean),
                      beta_rand_null_sd=apply(res_simu_beta,1,sd))

## beta dev calculated from randomization (following Kraft et al. 2011):
res_Gentry$beta_rand_dev=with(res_Gentry,
                              (beta_obs-beta_rand_null)/beta_rand_null_sd)

## beta dev calculated using the analytical model:
## note that the SAD of the 17th plot, which has only 5 spp and 134 ind, cannot be fitted using the mete logseries. 
res_Gentry$beta_null=with(res_Gentry, mapply(beta_null,M=M,S=S,N=N))
res_Gentry$beta_null_sd=with(res_Gentry, mapply(beta_null_sd,M=M,S=S,N=N))
res_Gentry$beta_dev=with(res_Gentry, (beta_obs-beta_null)/beta_null_sd)

res_Gentry$beta_dev_sd=with(res_Gentry, mapply(beta_dev_sd,M=M,S=S,N=N))

res_Gentry$kappa=sapply(comp.dat,kappa_est)
res_Gentry$beta_NBD=with(res_Gentry, mapply(beta_NBD,M=M,S=S,N=N,kappa=kappa))
res_Gentry$beta_NBD_sd1=with(res_Gentry, mapply(beta_NBD_sd,M=M,S=S,N=N,kappa=kappa,level=1))
res_Gentry$beta_NBD_sd2=with(res_Gentry, mapply(beta_NBD_sd,M=M,S=S,N=N,kappa=kappa,level=2))

## KS test of the empirical SADs vs. logseries
res_Gentry$ks_pvalue=sapply(comp.dat,sad_ks_test)

## ploting the results -------------------------------
## fig 1 beta_null and SD_null --------
win.metafile('figs/beta_null.wmf',8,4)
par(mfrow=c(1,2),mar=c(4,3.5,1,1),mgp=c(2,.5,0))
with(res_Gentry,plot(beta_null~beta_rand_null,xlim= lim<- range(na.omit(cbind(beta_null,beta_rand_null))),ylim=lim,xlab=expression(paste('Randomized ',beta['null'])),ylab=expression(paste('Analytical ',beta['null'])),typ='n'));abline(0,1,lty='dashed',lwd=1)

with(res_Gentry,points(beta_rand_null,beta_null,pch=21,bg='white',lwd=.1,cex=(log(S)-min(log(S)))/(max(log(S))-min(log(S)))*2+.1))
r2=format(with(res_Gentry,R2(beta_rand_null,beta_null)),digits=3,nsmall=3)
rho2=format(with(res_Gentry,cor(beta_rand_null,beta_null,use='comp')^2),digits=3,nsmall=3)
text(par('usr')[2], par('usr')[3], substitute(paste(italic(R[1:1]^2),' = ',r2,';  ',italic(R[x:y]^2),' = ',rho2),list(r2=r2,rho2=rho2)),adj=c(1.05,-.5))
text(par('usr')[1], par('usr')[4], '(a)',adj=c(-.5,1.5),font=2)

with(res_Gentry,plot(beta_null_sd~beta_rand_null_sd,xlim= lim<- range(na.omit(cbind(beta_null_sd,beta_rand_null_sd))),ylim=lim,xlab=expression(paste('Randomized ',SD['null'])),ylab=expression(paste('Analytical ',SD['null'])),type='n'));abline(0,1,lty='dashed',lwd=1)
with(res_Gentry,points(beta_rand_null_sd,beta_null_sd,pch=21,bg='white',lwd=.1,cex=(log(S)-min(log(S)))/(max(log(S))-min(log(S)))*2+.1))
r2=format(with(res_Gentry,R2(beta_rand_null_sd,beta_null_sd)),digits=3,nsmall=3)
rho2=format(with(res_Gentry,cor(beta_rand_null_sd,beta_null_sd,use='comp')^2),digits=3,nsmall=3)
text(par('usr')[2], par('usr')[3], substitute(paste(italic(R[1:1]^2),' = ',r2,';  ',italic(R[x:y]^2),' = ',rho2),list(r2=r2,rho2=rho2)),adj=c(1.05,-.5))
text(par('usr')[1], par('usr')[4], '(b)',adj=c(-.5,1.5),font=2)
dev.off()

## fig 2 beta_dev --------
win.metafile('figs/beta_dev.wmf',2.5,6)
par(mfcol=c(3,1),mar=c(3.2,3.2,.5,.5),mgp=c(1.8,.5,0))
with(res_Gentry,plot(beta_dev~beta_rand_dev,xlim= lim<- range(na.omit(cbind(beta_dev,beta_rand_dev))),ylim=c(-10,30),xlab=expression(paste('Randomized ',italic(beta)['dev'])),ylab=expression(paste('Analytical ', italic(beta)['dev'])),type='n'));abline(0,1,lty='dashed',lwd=1)
with(res_Gentry,arrows(beta_rand_dev,beta_dev+1.96*beta_dev_sd, beta_rand_dev, beta_dev-1.96*beta_dev_sd,code=3,angle = 90,length=.02,col=grey(.7)),lwd=.3)
sig_id=with(res_Gentry,(abs(beta_dev-beta_rand_dev)>1.96*beta_dev_sd))
with(res_Gentry,points(beta_rand_dev,beta_dev,pch=21,bg=c('white','black')[sig_id+1]),lwd=.3)
r2=format(with(res_Gentry,R2(beta_rand_dev,beta_dev)),digits=3,nsmall=3)
# for !sig_id, r2 increases to -.0245,rho^2 increases to .541
rho2=format(with(res_Gentry,cor(beta_rand_dev,beta_dev,use='comp')^2),digits=3,nsmall=3)
text(par('usr')[2], par('usr')[3], substitute(paste(italic(R[1:1]^2),' = ',r2,';  ',italic(R[x:y]^2),' = ',rho2),list(r2=r2,rho2=rho2)),adj=c(1.05,-.5),cex=.8)
text(par('usr')[1], par('usr')[4], '(a)',adj=c(-.5,1.5),font=2)
## analytical beta_dev ~ lat
with(res_Gentry,plot(beta_dev~lat,type='n',xlab='Absolute latitude (°)',ylab=expression(paste('Analytical ',beta['dev']))))
with(res_Gentry,points(beta_dev~lat,pch=21,bg=c('white','black')[sig_id+1]))

rho2=format(with(res_Gentry,cor(beta_dev,lat,use='comp')^2), digits=3,nsmall=3)
p=round(with(res_Gentry,cor.test(beta_dev,lat,use='comp')$p.v),3)
p0=p;p[p0<.001]=' < 0.001';p[p0>=.001]=paste(' = ',format(p0[p0>=.001],nsmall=3),sep='');rm(p0)

text(par('usr')[2], par('usr')[4], substitute(paste(italic(R[x:y]^2),' = ', rho2,', ',italic(P),p),list(rho2=rho2,p=p)),adj=c(1.05,1.5),cex=.8)
text(par('usr')[1], par('usr')[4], '(b)',adj=c(-.5,1.5),font=2)

with(res_Gentry,clip(min(lat),max(lat),par('usr')[3],par('usr')[4]))
with(res_Gentry,abline(lm(beta_dev~lat),col='black',lwd=2))
## randimized beta_dev ~ lat
with(res_Gentry,plot(beta_rand_dev~lat,type='n',xlab='Absolute latitude (°)',ylab=expression(paste('Randomized ', beta['dev']))))
with(res_Gentry,points(beta_rand_dev~lat,pch=21,bg=c('white','black')[sig_id+1]))

rho2=format(with(res_Gentry,round(cor(beta_rand_dev,lat,use='comp')^2,3)), nsmall=3)
p=round(with(res_Gentry,cor.test(beta_rand_dev,lat,use='comp')$p.v),3)
p0=p;p[p0<.001]=' < 0.001';p[p0>=.001]=paste(' = ',format(p0[p0>=.001],nsmall=3),sep='');rm(p0)

text(par('usr')[2], par('usr')[4], substitute(paste(italic(R[x:y]^2),' = ', rho2,', ',italic(P),p),list(rho2=rho2,p=p)),adj=c(1.05,1.5),cex=.8)
text(par('usr')[1], par('usr')[4], '(c)',adj=c(-.5,1.5),font=2)

with(res_Gentry,clip(min(lat),max(lat),par('usr')[3],par('usr')[4]))
with(res_Gentry,abline(lm(beta_rand_dev~lat),col='black',lwd=2))
dev.off()

## fig 3 beta_dev vs p-value of the ks test-------
win.metafile('figs/beta_dev_KStest.wmf',3.75,6)
par(mfcol=c(2,1),mar=c(3.2,3.2,.5,.5),mgp=c(1.8,.5,0))
hist(res_Gentry$ks_pvalue,breaks = 20,xlab='',main = '')
text(par('usr')[1], par('usr')[4], '(a)',adj=c(-1.5,1.5),font=2)

with(res_Gentry,plot(ks_pvalue,(beta_dev-beta_rand_dev),xlab='P-value of Kolmogorov-Smirnov (KS) test',ylab=expression(paste('Analytical ', italic(beta)['dev'], ' – Randomized ',italic(beta)['dev'])),type='n'));abline(v=.05,lty='dashed',lwd=1)
sig_id=with(res_Gentry,(abs(beta_dev-beta_rand_dev)>1.96*beta_dev_sd))
with(res_Gentry,points(ks_pvalue,(beta_dev-beta_rand_dev),pch=21,bg=c('white','black')[sig_id+1]),lwd=.3)
text(par('usr')[1], par('usr')[4], '(b)',adj=c(-1.5,1.5),font=2)
dev.off()
## fig 4 beta_NBD---------
win.metafile('figs/beta_NBD.wmf',5,5)
par(mfrow=c(1,1),mar=c(3.2,3.2,.8,.8),mgp=c(1.8,.5,0))
with(res_Gentry,plot(beta_NBD~beta_obs,typ='n',xlim= lim<- range(na.omit(cbind(beta_NBD,beta_obs))),ylim=lim,xlab=expression(paste(beta['obs'])),ylab=expression(paste(beta['NBD']))));abline(0,1,lty='dashed',lwd=1)

with(res_Gentry,arrows(beta_obs,beta_NBD-1.96*beta_NBD_sd2, beta_obs,beta_NBD+1.96*beta_NBD_sd2,code=3,angle = 90,length=.02,lwd=1,col=grey(.7),lty=2))
with(res_Gentry,arrows(beta_obs,beta_NBD-1.96*beta_NBD_sd1, beta_obs,beta_NBD+1.96*beta_NBD_sd1,code=3,angle = 90,length=.0,lwd=1))
with(res_Gentry,points(beta_obs,beta_NBD,pch=21,bg='white',lwd=.3,cex=1))
sig_id=with(res_Gentry, which(abs(beta_NBD-beta_obs)>1.96*beta_NBD_sd2))
with(res_Gentry[sig_id,],points(beta_obs,beta_NBD,pch=21,bg='black',lwd=.3,cex=1))


r2=format(with(res_Gentry,R2(beta_obs,beta_NBD)),digits=3,nsmall=3)
rho2=format(with(res_Gentry,cor(beta_obs,beta_NBD,use='comp')^2),digits=3,nsmall=3)
text(par('usr')[2], par('usr')[3], substitute(paste(italic(R[1:1]^2),' = ',r2,';  ',italic(R[x:y]^2),' = ',rho2),list(r2=r2,rho2=rho2)),adj=c(1.05,-.5))

dev.off()



## fig 5 scaling of beta_dev ------
bennett=read.csv('data/Bennett&Gilbert.csv')
win.metafile('figs/beta_dev_scaling.wmf',5,5)
par(mfrow=c(1,1),mgp=c(1.8,.5,0),mar=c(3.2,3.2,.8,.8))
plot(beta_dev~nplots,bennett,log='xy',type='n',xlab='Sampling effort (Number of quadrats)', ylab=expression(paste(beta['dev'])), pch=as.numeric(dataset))

ab=matrix(NA,nrow = 4, ncol=2)#container of the coefs; 1st col a; 2nd col b

y=bennett$beta_dev[bennett$dataset=='meadow_fig1b']
x=bennett$nplots[bennett$dataset=='meadow_fig1b']
y_hat=sqrt(x)*(y[length(x)]/sqrt(x[length(x)]))
lines(x,y_hat,col='grey',lwd=2)
points(x,y,pch=1,col=1)

y=bennett$beta_dev[bennett$dataset=='ab_field_fig1f']
x=bennett$nplots[bennett$dataset=='ab_field_fig1f']
y_hat=sqrt(x)*(y[length(x)]/sqrt(x[length(x)]))
lines(x,y_hat,col='grey',lwd=2)
points(x,y,pch=2,col=2)

y=bennett$beta_dev[bennett$dataset=='forest_figS1b']
x=bennett$nplots[bennett$dataset=='forest_figS1b']
y_hat=sqrt(x)*(y[length(x)]/sqrt(x[length(x)]))
lines(x,y_hat,col='grey',lwd=2)
points(x,y,pch=3,col=3)

y=bennett$beta_dev[bennett$dataset=='diatom_figS1f']
x=bennett$nplots[bennett$dataset=='diatom_figS1f']
y_hat=sqrt(x)*(y[length(x)]/sqrt(x[length(x)]))
lines(x,y_hat,col='grey',lwd=2)
points(x,y,pch=4,col=4)

id=c('Meadow plants (Fig.1b)     ', 'Abandoned field plants (Fig.1f)','Forest understory (Fig.S1b)','Lake diatoms (Fig.S1f)        ')
legend('bottomright', legend = id ,bty = 'n', cex=1,pch=1:4, col=1:4,text.col = 1:4)
dev.off()
