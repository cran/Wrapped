
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
		ff=function (x) {pwrappedg(x, spec, K=K, ...)-p[i]}
        	qq[i]=uniroot(ff,lower=-pi,upper=pi)$root
        }
        return(qq)
}







mwrappedg<-function(g, data, starts, K=K, method="BFGS"){

if(g!="norm" & g!="gumbel" & g!="logis" & g!="t.scaled" & g!="cauchy" &
g!="sn" & g!="ALD" & g!="nl" & g!="skewlap" & g!="glogis" & g!="sl"
& g!="normp" & g!="SEP1" & g!="SEP2" & g!="SEP3" & g!="SEP4" & g!="NET" &
g!="SN2" & g!="sgt" & g!="skewhyp" & g!="kiener1" & g!="st" & g!="laplace" &
g!="skewlaplace" & g!="asl" & g!="asla" & g!="al" & g!="PTL" & g!="gat" &
g!="vg" & g!="nig" & g!="sc" & g!="slash" & g!="exGAUS" & g!="lgamma" & g!="ST1"
& g!="ST3" & g!="ST4" & g!="ST5" & g!="SHASH" & g!="SHASHo" & g!="SHASH2" & g!="EGB2"
& g!="JSU" & g!="JUSo" & g!="betanorm" & g!="tikuv" & g!="hyperb")
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

#metRology
if(g=="t.scaled"){
den=function(par,x){mean=par[1]; sd=par[2]; df=par[3]; dwrappedg(x,g,K=K,mean=mean,sd=sd,df=df)}
cum=function(par,x){mean=par[1]; sd=par[2]; df=par[3]; pwrappedg(x,g,K=K,mean=mean,sd=sd,df=df)}
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
if(g=="st"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; df=par[4]; dwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha,nu=df)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; df=par[4]; pwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha,nu=df)}
}



#sc - sn
if(g=="sc"){
den=function(par,x){xi=par[1]; omega=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,xi=xi,omega=omega,alpha=alpha)}
cum=function(par,x){xi=par[1]; omega=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,xi=xi,omega=omega,alpha=alpha)}
}



#asymmetric laplace - ald
if(g=="ALD"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,p=alpha)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,p=alpha)}
}




#normal Laplace - NormalLaplace
if(g=="nl"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
}


#Generalized logistic - glogis
if(g=="glogis"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,location=mean,scale=sd,shape=alpha)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,location=mean,scale=sd,shape=alpha)}
}



#Skew logistic - sld
if(g=="sl"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; dwrappedg(x,g,K=K,parameters=c(mean,sd,alpha))}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; pwrappedg(x,g,K=K,parameters=c(mean,sd,alpha))}
}



#exponential power - normalp
if(g=="normp"){
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



#exGAUS  - gamlss.dist
if(g=="exGAUS"){
den=function(par,x){mu=par[1]; sigma=par[2]; nu=par[3]; dwrappedg(x,g,K=K,mu=mu,sigma=sigma,nu=nu)}
cum=function(par,x){mu=par[1]; sigma=par[2]; nu=par[3]; pwrappedg(x,g,K=K,mu=mu,sigma=sigma,nu=nu)}
}


#ST1 - gamlss.dist
if(g=="ST1"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#ST3 - gamlss.dist
if(g=="ST3"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#ST4 - gamlss.dist
if(g=="ST4"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#ST5 - gamlss.dist
if(g=="ST5"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}




#SHASH - gamlss.dist
if(g=="SHASH"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}



#SHASHo - gamlss.dist
if(g=="SHASHo"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SHASH2 - gamlss.dist
if(g=="SHASH2"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#EGB2 - gamlss.dist
if(g=="EGB2"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#JSU - gamlss.dist
if(g=="JSU"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#JSUo - gamlss.dist
if(g=="JSUo"){
den=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1]; sd=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
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


#Hyperbolic - GeneralizedHyperbolic
if(g=="hyperb"){
den=function(par,x){mu=par[1]; delta=par[2]; alpha=par[3]; beta=par[4]; dwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
cum=function(par,x){mu=par[1]; delta=par[2]; alpha=par[3]; beta=par[4]; pwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
}



#Skew Laplace - GeneralizedHyperbolic
if(g=="skewlap"){
den=function(par,x){mean=par[1]; alpha=par[2]; beta=par[3]; dwrappedg(x,g,K=K,mu=mean,alpha=alpha,beta=beta)}
cum=function(par,x){mean=par[1]; alpha=par[2]; beta=par[3]; pwrappedg(x,g,K=K,mu=mean,alpha=alpha,beta=beta)}
}



#Slash - VGAM
if(g=="slash"){
den=function(par,x){mu=par[1]; sigma=par[2]; dwrappedg(x,g,K=K,mu=mu,sigma=sigma)}
cum=function(par,x){mu=par[1]; sigma=par[2]; pwrappedg(x,g,K=K,mu=mu,sigma=sigma)}
}


#betanorm - VGAM
if(g=="betanorm"){
den=function(par,x){mean=par[1]; sd=par[2]; shape1=par[3]; shape2=par[4]; dwrappedg(x,g,K=K,mean=mean,sd=sd,shape1=shape1,shape2=shape2)}
cum=function(par,x){mean=par[1]; sd=par[2]; shape1=par[3]; shape2=par[4]; pwrappedg(x,g,K=K,mean=mean,sd=sd,shape1=shape1,shape2=shape2)}
}



#tikuv - VGAM
if(g=="tikuv"){
den=function(par,x){mean=par[1]; sigma=par[2]; d=par[3]; dwrappedg(x,g,K=K,mean=mean,sigma=sigma,d=d)}
cum=function(par,x){mean=par[1]; sigma=par[2]; d=par[3]; pwrappedg(x,g,K=K,mean=mean,sigma=sigma,d=d)}
}


#laplace - VGAM
if(g=="laplace"){
den=function(par,x){mean=par[1]; sigma=par[2]; dwrappedg(x,g,K=K,location=mean,scale=sigma)}
cum=function(par,x){mean=par[1]; sigma=par[2]; pwrappedg(x,g,K=K,location=mean,scale=sigma)}
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





plotfour<-function(g,K=K,para,plotit){

if(g=="norm"){
den=function(par,x){mean=par[1];sd=par[2];dwrappedg(x,g,K=K,mean=mean,sd=sd)}
cum=function(par,x){mean=par[1];sd=par[2];pwrappedg(x,g,K=K,mean=mean,sd=sd)}
qut=function(par,x){mean=par[1];sd=par[2];qwrappedg(x,g,K=K,mean=mean,sd=sd)}
ran=function(par,n){mean=par[1];sd=par[2];rwrappedg(n,g,mean=mean,sd=sd)}
}

#gumbel-evd
if(g=="gumbel"){
den=function(par,x){mean=par[1];sd=par[2];dwrappedg(x,g,K=K,loc=mean,scale=sd)}
cum=function(par,x){mean=par[1];sd=par[2];pwrappedg(x,g,K=K,loc=mean,scale=sd)}
qut=function(par,x){mean=par[1];sd=par[2];qwrappedg(x,g,K=K,loc=mean,scale=sd)}
ran=function(par,n){mean=par[1];sd=par[2];rwrappedg(n,g,loc=mean,scale=sd)}
}


if(g=="logis"){
den=function(par,x){mean=par[1];sd=par[2];dwrappedg(x,g,K=K,location=mean,scale=sd)}
cum=function(par,x){mean=par[1];sd=par[2];pwrappedg(x,g,K=K,location=mean,scale=sd)}
qut=function(par,x){mean=par[1];sd=par[2];qwrappedg(x,g,K=K,location=mean,scale=sd)}
ran=function(par,n){mean=par[1];sd=par[2];rwrappedg(n,g,location=mean,scale=sd)}
}

#metRology
if(g=="t.scaled"){
den=function(par,x){mean=par[1];sd=par[2];df=par[3];dwrappedg(x,g,K=K,mean=mean,sd=sd,df=df)}
cum=function(par,x){mean=par[1];sd=par[2];df=par[3];pwrappedg(x,g,K=K,mean=mean,sd=sd,df=df)}
qut=function(par,x){mean=par[1];sd=par[2];df=par[3];qwrappedg(x,g,K=K,mean=mean,sd=sd,df=df)}
ran=function(par,n){mean=par[1];sd=par[2];df=par[3];rwrappedg(n,g,mean=mean,sd=sd,df=df)}
}



if(g=="cauchy"){
den=function(par,x){mean=par[1];sd=par[2];dwrappedg(x,g,K=K,location=mean,scale=sd)}
cum=function(par,x){mean=par[1];sd=par[2];pwrappedg(x,g,K=K,location=mean,scale=sd)}
qut=function(par,x){mean=par[1];sd=par[2];qwrappedg(x,g,K=K,location=mean,scale=sd)}
ran=function(par,n){mean=par[1];sd=par[2];rwrappedg(n,g,location=mean,scale=sd)}
}


#skewnormal-sn
if(g=="sn"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];dwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];pwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];qwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];rwrappedg(n,g,xi=mean,omega=sd,alpha=alpha)}
}



#skewt-sn
if(g=="st"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];df=par[4];dwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha,nu=df)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];df=par[4];pwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha,nu=df)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];df=par[4];qwrappedg(x,g,K=K,xi=mean,omega=sd,alpha=alpha,nu=df)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];df=par[4];rwrappedg(n,g,xi=mean,omega=sd,alpha=alpha,nu=df)}
}



#sc-sn
if(g=="sc"){
den=function(par,x){xi=par[1];omega=par[2];alpha=par[3];dwrappedg(x,g,K=K,xi=xi,omega=omega,alpha=alpha)}
cum=function(par,x){xi=par[1];omega=par[2];alpha=par[3];pwrappedg(x,g,K=K,xi=xi,omega=omega,alpha=alpha)}
qut=function(par,x){xi=par[1];omega=par[2];alpha=par[3];qwrappedg(x,g,K=K,xi=xi,omega=omega,alpha=alpha)}
ran=function(par,n){xi=par[1];omega=par[2];alpha=par[3];rwrappedg(n,g,xi=xi,omega=omega,alpha=alpha)}
}



#asymmetriclaplace-ald
if(g=="ALD"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];dwrappedg(x,g,K=K,mu=mean,sigma=sd,p=alpha)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];pwrappedg(x,g,K=K,mu=mean,sigma=sd,p=alpha)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];qwrappedg(x,g,K=K,mu=mean,sigma=sd,p=alpha)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];rwrappedg(n,g,mu=mean,sigma=sd,p=alpha)}
}




#normalLaplace-NormalLaplace
if(g=="nl"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,alpha=alpha,beta=beta)}
}


#Generalizedlogistic-glogis
if(g=="glogis"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];dwrappedg(x,g,K=K,location=mean,scale=sd,shape=alpha)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];pwrappedg(x,g,K=K,location=mean,scale=sd,shape=alpha)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];qwrappedg(x,g,K=K,location=mean,scale=sd,shape=alpha)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];rwrappedg(n,g,location=mean,scale=sd,shape=alpha)}
}



#Skewlogistic-sld
if(g=="sl"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];dwrappedg(x,g,K=K,parameters=c(mean,sd,alpha))}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];pwrappedg(x,g,K=K,parameters=c(mean,sd,alpha))}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];qwrappedg(x,g,K=K,parameters=c(mean,sd,alpha))}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];rwrappedg(n,g,parameters=c(mean,sd,alpha))}
}



#exponentialpower-normalp
if(g=="normp"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];dwrappedg(x,g,K=K,mu=mean,sigmap=sd,p=alpha)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];pwrappedg(x,g,K=K,mu=mean,sigmap=sd,p=alpha)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];qwrappedg(x,g,K=K,mu=mean,sigmap=sd,p=alpha)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];rwrappedg(n,g,mu=mean,sigmap=sd,p=alpha)}
}



#SEP1-gamlss.dist
if(g=="SEP1"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SEP2-gamlss.dist
if(g=="SEP2"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SEP3-gamlss.dist
if(g=="SEP3"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SEP4-gamlss.dist
if(g=="SEP4"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}



#NET-gamlss.dist
if(g=="NET"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SN2-gamlss.dist
if(g=="SN2"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha)}
}



#exGAUS-gamlss.dist
if(g=="exGAUS"){
den=function(par,x){mu=par[1];sigma=par[2];nu=par[3];dwrappedg(x,g,K=K,mu=mu,sigma=sigma,nu=nu)}
cum=function(par,x){mu=par[1];sigma=par[2];nu=par[3];pwrappedg(x,g,K=K,mu=mu,sigma=sigma,nu=nu)}
qut=function(par,x){mu=par[1];sigma=par[2];nu=par[3];qwrappedg(x,g,K=K,mu=mu,sigma=sigma,nu=nu)}
ran=function(par,n){mu=par[1];sigma=par[2];nu=par[3];rwrappedg(n,g,mu=mu,sigma=sigma,nu=nu)}
}


#ST1-gamlss.dist
if(g=="ST1"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#ST3-gamlss.dist
if(g=="ST3"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#ST4-gamlss.dist
if(g=="ST4"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#ST5-gamlss.dist
if(g=="ST5"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}




#SHASH-gamlss.dist
if(g=="SHASH"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}



#SHASHo-gamlss.dist
if(g=="SHASHo"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#SHASH2-gamlss.dist
if(g=="SHASH2"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#EGB2-gamlss.dist
if(g=="EGB2"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#JSU-gamlss.dist
if(g=="JSU"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#JSUo-gamlss.dist
if(g=="JSUo"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,sigma=sd,nu=alpha,tau=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,sigma=sd,nu=alpha,tau=beta)}
}


#sgt-sgt
if(g=="sgt"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];gamma=par[5];dwrappedg(x,g,K=K,mu=mean,sigma=sd,lambda=alpha,p=beta,q=gamma)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];gamma=par[5];pwrappedg(x,g,K=K,mu=mean,sigma=sd,lambda=alpha,p=beta,q=gamma)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];gamma=par[5];qwrappedg(x,g,K=K,mu=mean,sigma=sd,lambda=alpha,p=beta,q=gamma)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];gamma=par[5];rwrappedg(n,g,mu=mean,sigma=sd,lambda=alpha,p=beta,q=gamma)}
}


#skewhyperbolic-SkewHyperbolic
if(g=="skewhyp"){
den=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mean,delta=sd,beta=alpha,nu=beta)}
cum=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mean,delta=sd,beta=alpha,nu=beta)}
qut=function(par,x){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mean,delta=sd,beta=alpha,nu=beta)}
ran=function(par,n){mean=par[1];sd=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mean,delta=sd,beta=alpha,nu=beta)}
}



#AsymmetricLaplaceDistribution-cubfits
if(g=="asl"){
den=function(par,x){theta=par[1];mu=par[2];sigma=par[3];dwrappedg(x,g,K=K,theta=theta,mu=mu,sigma=sigma)}
cum=function(par,x){theta=par[1];mu=par[2];sigma=par[3];pwrappedg(x,g,K=K,theta=theta,mu=mu,sigma=sigma)}
qut=function(par,x){theta=par[1];mu=par[2];sigma=par[3];qwrappedg(x,g,K=K,theta=theta,mu=mu,sigma=sigma)}
ran=function(par,n){theta=par[1];mu=par[2];sigma=par[3];rwrappedg(n,g,theta=theta,mu=mu,sigma=sigma)}
}


#AsymmetricLaplaceDistribution-cubfits
if(g=="asla"){
den=function(par,x){theta=par[1];kappa=par[2];sigma=par[3];dwrappedg(x,g,K=K,theta=theta,kappa=kappa,sigma=sigma)}
cum=function(par,x){theta=par[1];kappa=par[2];sigma=par[3];pwrappedg(x,g,K=K,theta=theta,kappa=kappa,sigma=sigma)}
qut=function(par,x){theta=par[1];kappa=par[2];sigma=par[3];qwrappedg(x,g,K=K,theta=theta,kappa=kappa,sigma=sigma)}
ran=function(par,n){theta=par[1];kappa=par[2];sigma=par[3];rwrappedg(n,g,theta=theta,kappa=kappa,sigma=sigma)}
}


#dal-lqmm
if(g=="al"){
den=function(par,x){mu=par[1];sigma=par[2];tau=par[3];dwrappedg(x,g,K=K,mu=mu,sigma=sigma,tau=tau)}
cum=function(par,x){mu=par[1];sigma=par[2];tau=par[3];pwrappedg(x,g,K=K,mu=mu,sigma=sigma,tau=tau)}
qut=function(par,x){mu=par[1];sigma=par[2];tau=par[3];qwrappedg(x,g,K=K,mu=mu,sigma=sigma,tau=tau)}
ran=function(par,n){mu=par[1];sigma=par[2];tau=par[3];rwrappedg(n,g,mu=mu,sigma=sigma,tau=tau)}
}

#PTL-LCA
if(g=="PTL"){
den=function(par,x){alpha=par[1];beta=par[2];gamma=par[3];dwrappedg(x,g,K=K,alpha=alpha,beta=beta,gamma=gamma)}
cum=function(par,x){alpha=par[1];beta=par[2];gamma=par[3];dwrappedg(x,g,K=K,alpha=alpha,beta=beta,gamma=gamma)}
qut=function(par,x){alpha=par[1];beta=par[2];gamma=par[3];qwrappedg(x,g,K=K,alpha=alpha,beta=beta,gamma=gamma)}
ran=function(par,n){alpha=par[1];beta=par[2];gamma=par[3];rwrappedg(n,g,alpha=alpha,beta=beta,gamma=gamma)}
}


#gat-GEVStableGarch
if(g=="gat"){
den=function(par,x){mean=par[1];sd=par[2];nu=par[3];d=par[4];xi=par[5];dwrappedg(x,g,K=K,mean=mean,sd=sd,nu=nu,d=d,xi=xi)}
cum=function(par,x){mean=par[1];sd=par[2];nu=par[3];d=par[4];xi=par[5];pwrappedg(x,g,K=K,mean=mean,sd=sd,nu=nu,d=d,xi=xi)}
qut=function(par,x){mean=par[1];sd=par[2];nu=par[3];d=par[4];xi=par[5];qwrappedg(x,g,K=K,mean=mean,sd=sd,nu=nu,d=d,xi=xi)}
ran=function(par,n){mean=par[1];sd=par[2];nu=par[3];d=par[4];xi=par[5];rwrappedg(n,g,mean=mean,sd=sd,nu=nu,d=d,xi=xi)}
}


#VarianceGammaDistribution-VarianceGamma
if(g=="vg"){
den=function(par,x){vgC=par[1];sigma=par[2];theta=par[3];nu=par[4];dwrappedg(x,g,K=K,vgC=vgC,sigma=sigma,theta=theta,nu=nu)}
cum=function(par,x){vgC=par[1];sigma=par[2];theta=par[3];nu=par[4];pwrappedg(x,g,K=K,vgC=vgC,sigma=sigma,theta=theta,nu=nu)}
qut=function(par,x){vgC=par[1];sigma=par[2];theta=par[3];nu=par[4];qwrappedg(x,g,K=K,vgC=vgC,sigma=sigma,theta=theta,nu=nu)}
ran=function(par,n){vgC=par[1];sigma=par[2];theta=par[3];nu=par[4];rwrappedg(n,g,vgC=vgC,sigma=sigma,theta=theta,nu=nu)}
}



#NIG-GeneralizedHyperbolic
if(g=="nig"){
den=function(par,x){mu=par[1];delta=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
cum=function(par,x){mu=par[1];delta=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
qut=function(par,x){mu=par[1];delta=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
ran=function(par,n){mu=par[1];delta=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mu,delta=delta,alpha=alpha,beta=beta)}
}


#Hyperbolic-GeneralizedHyperbolic
if(g=="hyperb"){
den=function(par,x){mu=par[1];delta=par[2];alpha=par[3];beta=par[4];dwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
cum=function(par,x){mu=par[1];delta=par[2];alpha=par[3];beta=par[4];pwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
qut=function(par,x){mu=par[1];delta=par[2];alpha=par[3];beta=par[4];qwrappedg(x,g,K=K,mu=mu,delta=delta,alpha=alpha,beta=beta)}
ran=function(par,n){mu=par[1];delta=par[2];alpha=par[3];beta=par[4];rwrappedg(n,g,mu=mu,delta=delta,alpha=alpha,beta=beta)}
}



#SkewLaplace-GeneralizedHyperbolic
if(g=="skewlap"){
den=function(par,x){mean=par[1];alpha=par[2];beta=par[3];dwrappedg(x,g,K=K,mu=mean,alpha=alpha,beta=beta)}
cum=function(par,x){mean=par[1];alpha=par[2];beta=par[3];pwrappedg(x,g,K=K,mu=mean,alpha=alpha,beta=beta)}
qut=function(par,x){mean=par[1];alpha=par[2];beta=par[3];qwrappedg(x,g,K=K,mu=mean,alpha=alpha,beta=beta)}
ran=function(par,n){mean=par[1];alpha=par[2];beta=par[3];rwrappedg(n,g,mu=mean,alpha=alpha,beta=beta)}
}



#Slash-VGAM
if(g=="slash"){
den=function(par,x){mu=par[1];sigma=par[2];dwrappedg(x,g,K=K,mu=mu,sigma=sigma)}
cum=function(par,x){mu=par[1];sigma=par[2];pwrappedg(x,g,K=K,mu=mu,sigma=sigma)}
qut=function(par,x){mu=par[1];sigma=par[2];qwrappedg(x,g,K=K,mu=mu,sigma=sigma)}
ran=function(par,n){mu=par[1];sigma=par[2];rwrappedg(n,g,mu=mu,sigma=sigma)}
}


#betanorm-VGAM
if(g=="betanorm"){
den=function(par,x){mean=par[1];sd=par[2];shape1=par[3];shape2=par[4];dwrappedg(x,g,K=K,mean=mean,sd=sd,shape1=shape1,shape2=shape2)}
cum=function(par,x){mean=par[1];sd=par[2];shape1=par[3];shape2=par[4];pwrappedg(x,g,K=K,mean=mean,sd=sd,shape1=shape1,shape2=shape2)}
qut=function(par,x){mean=par[1];sd=par[2];shape1=par[3];shape2=par[4];qwrappedg(x,g,K=K,mean=mean,sd=sd,shape1=shape1,shape2=shape2)}
ran=function(par,n){mean=par[1];sd=par[2];shape1=par[3];shape2=par[4];rwrappedg(n,g,mean=mean,sd=sd,shape1=shape1,shape2=shape2)}
}



#tikuv-VGAM
if(g=="tikuv"){
den=function(par,x){mean=par[1];sigma=par[2];d=par[3];dwrappedg(x,g,K=K,mean=mean,sigma=sigma,d=d)}
cum=function(par,x){mean=par[1];sigma=par[2];d=par[3];pwrappedg(x,g,K=K,mean=mean,sigma=sigma,d=d)}
qut=function(par,x){mean=par[1];sigma=par[2];d=par[3];qwrappedg(x,g,K=K,mean=mean,sigma=sigma,d=d)}
ran=function(par,n){mean=par[1];sigma=par[2];d=par[3];rwrappedg(n,g,mean=mean,sigma=sigma,d=d)}
}


#laplace-VGAM
if(g=="laplace"){
den=function(par,x){mean=par[1];sigma=par[2];dwrappedg(x,g,K=K,location=mean,scale=sigma)}
cum=function(par,x){mean=par[1];sigma=par[2];pwrappedg(x,g,K=K,location=mean,scale=sigma)}
qut=function(par,x){mean=par[1];sigma=par[2];qwrappedg(x,g,K=K,location=mean,scale=sigma)}
ran=function(par,n){mean=par[1];sigma=par[2];rwrappedg(n,g,location=mean,scale=sigma)}
}






#lgamma-ordinal
if(g=="lgamma"){
den=function(par,x){lambda=par[1];dwrappedg(x,g,K=K,lambda=lambda)}
cum=function(par,x){lambda=par[1];pwrappedg(x,g,K=K,lambda=lambda)}
qut=function(par,x){lambda=par[1];qwrappedg(x,g,K=K,lambda=lambda)}
ran=function(par,n){lambda=par[1];rwrappedg(n,g,K=K,lambda=lambda)}
}


par(mfrow=c(2,2),oma=rep(0,4))

if(plotit=="pdf")
{x=-pi+2*pi*seq(0.01,0.99,0.01)
y1=den(para[[1]],x)
y2=den(para[[2]],x)
y3=den(para[[3]],x)
y4=den(para[[4]],x)
plot(x,y1,xlab="x",ylab="PDF",type="l",main=toString(para[[1]]))
plot(x,y2,xlab="x",ylab="PDF",type="l",main=toString(para[[2]]))
plot(x,y3,xlab="x",ylab="PDF",type="l",main=toString(para[[3]]))
plot(x,y4,xlab="x",ylab="PDF",type="l",main=toString(para[[4]]))}



if(plotit=="cdf")
{x=-pi+2*pi*seq(0.01,0.99,0.01)
y1=cum(par=para[[1]],x)
y2=cum(par=para[[2]],x)
y3=cum(par=para[[3]],x)
y4=cum(par=para[[4]],x)
plot(x,y1,xlab="x",ylab="CDF",type="l",main=toString(para[[1]]))
plot(x,y2,xlab="x",ylab="CDF",type="l",main=toString(para[[2]]))
plot(x,y3,xlab="x",ylab="CDF",type="l",main=toString(para[[3]]))
plot(x,y4,xlab="x",ylab="CDF",type="l",main=toString(para[[4]]))}


if(plotit=="quantile")
{x=seq(0.01,0.99,0.01)
y1=qut(par=para[[1]],x)
y2=qut(par=para[[2]],x)
y3=qut(par=para[[3]],x)
y4=qut(par=para[[4]],x)
plot(x,y1,xlab="p",ylab="Quantile",type="l",main=toString(para[[1]]))
plot(x,y2,xlab="p",ylab="Quantile",type="l",main=toString(para[[2]]))
plot(x,y3,xlab="p",ylab="Quantile",type="l",main=toString(para[[3]]))
plot(x,y4,xlab="p",ylab="Quantile",type="l",main=toString(para[[4]]))}


if(plotit=="random")
{n=100
y1=rep(0,n)
y2=y1
y3=y1
y4=y1
yy=ran(par=para[[1]],n)
for (i in 1:n) {y1[i]=atan(yy[i,2]/yy[i,1])}
yy=ran(par=para[[2]],n)
for (i in 1:n) {y2[i]=atan(yy[i,2]/yy[i,1])}
yy=ran(par=para[[3]],n)
for (i in 1:n) {y3[i]=atan(yy[i,2]/yy[i,1])}
yy=ran(par=para[[4]],n)
for (i in 1:n) {y4[i]=atan(yy[i,2]/yy[i,1])}
hist(y1,xlab="Radius",ylab="Frequency",main=toString(para[[1]]),nclass=20)
hist(y2,xlab="Radius",ylab="Frequency",main=toString(para[[2]]),nclass=20)
hist(y3,xlab="Radius",ylab="Frequency",main=toString(para[[3]]),nclass=20)
hist(y4,xlab="Radius",ylab="Frequency",main=toString(para[[4]]),nclass=20)}


}
