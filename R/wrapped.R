
dwrappedg<-function(x, spec, K=K, ...)
{
        f<-function (x, ...) {do.call(paste("d",spec,sep=""),list(x, ...))}

	pdf<-rep(0,length(x))
        for (k in -K:K) pdf<-pdf+f(x+2*pi*k, ...)
        return(pdf)

}



pwrappedg<-function(x, spec, K=K, ...)
{
        F<-function (x, ...) {do.call(paste("p",spec,sep=""),list(x, ...))}

	cdf<-rep(0,length(x))
        for (k in -K:K) cdf<-cdf+F(x+2*pi*k, ...)-F(-pi+2*pi*k, ...)
        return(cdf)

}



rwrappedg<-function(n, spec, ...)
{
        rr<-function (n, ...) {do.call(paste("r",spec,sep=""),list(n, ...))}

	xx<-rr(n, ...)
        cc=cos(xx)
        ss=sin(xx)
        return(cbind(cc,ss))
}



qwrappedg<-function(p, spec, K=K, ...)
{
        fff<-function (x, ...) {pwrappedg(x, spec, K=K, ...)}
        
        qq=p
        for (i in 1:length(p))
        {
		ff=function (x) {fff(x)-p[i]}
        	qq[i]=uniroot(ff,lower=-pi,upper=pi)$root
        }
        return(qq)
}







mwrappedg<-function(g, data, starts, K=K, method="BFGS"){

if(g!="norm" & g!="gumbel" & g!="logis" & g!="t" & g!="cauchy" &
g!="sn" & g!="ald" & g!="nl" & g!="skewlap" & g!="glogis" & g!="4pl" & g!="sl"
& g!="normalp" & g!="SEP1" & g!="SEP2" & g!="SEP3" & g!="SEP4" & g!="NET" &
g!="SN2" & g!="sgt" & g!="skewhyp" & g!="ght" &  g!="logishp" & g!="kiener1" &
g!="laplacem" & g!="skewlaplace" & g!="asl" & g!="asla" & g!="al" & g!="PTL" & g!="gat" &
g!="vg" & g!="nig" & g!="sc" & g!="slash" & g!="exGAUS" & g!="lgamma")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="norm"){
den=function(par,x){mean=par[1]; sd=par[2]; dwrappedg(x,g,K=K,mean=mean,sd=sd)}
cum=function(par,x){mean=par[1]; sd=par[2]; pwrappedg(x,g,K=K,mean=mean,sd=sd)}
}

#gumbel - evd
if(g=="gumbel"){
den=function(par,x){mean=par[1]; sd=par[2]; dwrappedg(x,g,K=K,loc=mean,scale=sd)}
cum=function(par,x){mean=par[1]; sd=par[2]; pwrappedg(x,g,K=K,loc=mean,scale=sd)}
}


if(g=="logis"){
den=function(par,x){mean=par[1]; sd=par[2]; dwrappedg(x,g,K=K,location=mean,scale=sd)}
cum=function(par,x){mean=par[1]; sd=par[2]; pwrappedg(x,g,K=K,location=mean,scale=sd)}
}

if(g=="t"){
den=function(par,x){mean=par[1]; sd=par[2]; df=par[3]; (1/sd)*dwrappedg((x-mean)/sd,g,K=K,df=df)}
cum=function(par,x){mean=par[1]; sd=par[2]; df=par[3]; pwrappedg((x-mean)/sd,g,K=K,df=df)}
}



if(g=="cauchy"){
den=function(par,x){mean=par[1]; sd=par[2]; dwrappedg(x,g,K=K,location=mean,scale=sd)}
cum=function(par,x){mean=par[1]; sd=par[2]; pwrappedg(x,g,K=K,location=mean,scale=sd)}
}


#skew normal - sn
if(g=="sn"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha)}
}



#skew t - sn
if(g=="sn"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; df=par[4]; dwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha,nu=df)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; df=par[4]; pwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha,nu=df)}
}



#asymmetric laplace - ald
if(g=="ald"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,p=alpha)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,p=alpha)}
}




#normal laplace - NormalLaplace
if(g=="nl"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
}



#Skew Laplace - GeneralizedHyperbolic
if(g=="skewlap"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
}



#Generalized logistic - glogis
if(g=="glogis"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,location=mean,scale=sd,shape=alpha)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,location=mean,scale=sd,shape=alpha)}
}



#Four parameter logistic - irtProb
if(g=="4pl"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,a=mean,b=sd,c=alpha,d=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,a=mean,b=sd,c=alpha,d=beta)}
}

#Skew logistic - sld
if(g=="sl"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,parameters=c(mean,sd,alpha))}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,parameters=c(mean,sd,alpha))}
}



#exponential power - normalp
if(g=="normalp"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,mu=mean,sigmap=sd,p=alpha)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,mu=mean,sigmap=sd,p=alpha)}
}



#SEP1 - gamlss.dist
if(g=="SEP1"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SEP2 - gamlss.dist
if(g=="SEP2"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SEP3 - gamlss.dist
if(g=="SEP3"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SEP4 - gamlss.dist
if(g=="SEP4"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}



#NET - gamlss.dist
if(g=="NET"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SN2 - gamlss.dist
if(g=="SN2"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha)}
}


#sgt - sgt
if(g=="sgt"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; gamma=par[5]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,lambda=alpha,p=beta,q=gamma)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; gamma=par[5]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,lambda=alpha,p=beta,q=gamma)}
}


#skew hyperbolic - SkewHyperbolic
if(g=="skewhyp"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,delta=sd,beta=alpha,nu=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,delta=sd,beta=alpha,nu=beta)}
}



#generalized hyperbolic Student t - fBasics
if(g=="ght"){
den=function(par,x){beta=par[1]; delta=par[2]; mu=par[3]; nu=par[4]; dwrappedg(x,g,K=K,beta=beta,delta=delta,mu=mu,nu=nu)}
cum=function(par,x){beta=par[1]; delta=par[2]; mu=par[3]; nu=par[4]; pwrappedg(x,g,K=K,beta=beta,delta=delta,mu=mu,nu=nu)}
}


#logishp - FatTailsR
if(g=="logishp"){
den=function(par,x){k=k; dwrappedg(x,g,K=K,k=k)}
cum=function(par,x){k=k; pwrappedg(x,g,K=K,k=k)}
}



#Asymmetric Laplace Distribution - cubfits
if(g=="asl"){
den=function(par,x){theta=par[1]; mu=par[2]; sigma=par[3]; dwrappedg(x,g,K=K,theta=theta,mu=mu,sigma=sigma)}
cum=function(par,x){theta=par[1]; mu=par[2]; sigma=par[3]; pwrappedg(x,g,K=K,theta=theta,mu=mu,sigma=sigma)}
}


#Asymmetric Laplace Distribution - cubfits
if(g=="asla"){
den=function(par,x){theta=par[1]; kappa=par[2]; sigma=par[3]; dwrappedg(x,g,K=K,theta=theta,kappa=kappa,sigma=sigma)}
cum=function(par,x){theta=par[1]; kappa=par[2]; sigma=par[3]; pwrappedg(x,g,K=K,theta=theta,kappa=kappa,sigma=sigma)}
}


#dal - lqmm
if(g=="al"){
den=function(par,x){mu=par[1]; sigma=par[2]; tau=par[3]; dwrappedg(x,g,K=K,mu=mu,sigma=sigma,tau=tau)}
cum=function(par,x){mu=par[1]; sigma=par[2]; tau=par[3]; pwrappedg(x,g,K=K,mu=mu,sigma=sigma,tau=tau)}
}

#PTL - LCA
if(g=="PTL"){
den=function(par,x){alpha=par[1]; beta=par[2]; gamma=par[3]; dwrappedg(x,g,K=K,alpha=alpha,beta=beta,gamma=gamma)}
cum=function(par,x){alpha=par[1]; beta=par[2]; gamma=par[3]; dwrappedg(x,g,K=K,alpha=alpha,beta=beta,gamma=gamma)}
}


#gat - GEVStableGarch
if(g=="gat"){
den=function(par,x){mean=par[1]; sd=par[2]; nu=par[3]; d=par[4]; xi=par[5]; dwrappedg(x,g,K=K,mean=mean,sd=sd,nu=nu,d=d,xi=xi)}
cum=function(par,x){mean=par[1]; sd=par[2]; nu=par[3]; d=par[4]; xi=par[5]; pwrappedg(x,g,K=K,mean=mean,sd=sd,nu=nu,d=d,xi=xi)}
}


#VarianceGammaDistribution - VarianceGamma
if(g=="vg"){
den=function(par,x){vgC=par[1]; sigma=par[2]; theta=par[3]; nu=par[4]; dwrappedg(x,g,K=K,vgC=vgC,sigma=sigma,theta=theta,nu=nu)}
cum=function(par,x){vgC=par[1]; sigma=par[2]; theta=par[3]; nu=par[4]; pwrappedg(x,g,K=K,vgC=vgC,sigma=sigma,theta=theta,nu=nu)}
}



#NIG - GeneralizedHyperbolic
if(g=="nig"){
den=function(par,x){mu=par[1]; delta=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
cum=function(par,x){mu=par[1]; delta=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
}

#sc - sn
if(g=="sc"){
den=function(par,x){xi=par[1]; omega=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,xi=xi,omega=omega,alpha=alpha)}
cum=function(par,x){xi=par[1]; omega=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,xi=xi,omega=omega,alpha=alpha)}
}


#Slash - VGAM
if(g=="slash"){
den=function(par,x){mu=par[1]; sigma=par[2]; dwrappedg(x,g,K=K,mu=mu,sigma=sigma)}
cum=function(par,x){mu=par[1]; sigma=par[2]; pwrappedg(x,g,K=K,mu=mu,sigma=sigma)}
}


#exGAUS  - gamlss.dist
if(g=="exGAUS"){
den=function(par,x){mu=par[1]; sigma=par[2]; nu=par[3]; dwrappedg(x,g,K=K,mu=mu,sigma=sigma,nu=nu)}
cum=function(par,x){mu=par[1]; sigma=par[2]; nu=par[3]; pwrappedg(x,g,K=K,mu=mu,sigma=sigma,nu=nu)}
}



#lgamma - ordinal
if(g=="lgamma"){
den=function(par,x){lambda=par[1]; dwrappedg(x,g,K=K,lambda=lambda)}
cum=function(par,x){lambda=par[1]; pwrappedg(x,g,K=K,lambda=lambda)}
}



med=suppressWarnings(goodness.fit(pdf=den, cdf=cum, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else{"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}
