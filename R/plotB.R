"plotB" <-
function(y,dose,parm,sigma2,
	count=rep(1,length(y)),dgrid=sort(unique(c(seq(0,max(dose),length=50), dose))),
	predict= TRUE,plotDif=FALSE,plotMed=FALSE,
	plotResid=FALSE,clev=0.8,
	binary=c('no','logit','probit','BinRes'),BinResLev,
	BinResDir=c('>','<'),	activeControl=FALSE,ac,yac,dac,countac=rep(1,length(yac)),
	labac='Act Comp',shapeac=8,colac='red',
	symbol,symbolLabel='Group',symbolShape=8,symbolColor='red',symbolSize=4,
	xlim,ylim,xlab="Dose",ylab=ifelse(plotDif,"Diff with Comparator","Mean"),
	modelFun=emaxfun,makePlot=TRUE, ...) {

	ny<-length(y)
	
	if(missing(symbol)){
		nolegend<-TRUE
		symbol<-rep(1,ny)
	}else nolegend<-FALSE

	nsym<-length(symbol)
	if(nsym!=ny)stop('Symbol variable has invalid length')
	symbol<-as.character(symbol)  ### make sure symbol is character until output
	symlev<-sort(unique(symbol))
	nsymlev<-length(symlev)
	nshape<-length(symbolShape)
	if(nshape!=nsymlev){
		if(nshape==1){
			symbolShape<-rep(symbolShape,nsymlev)
		}else if(nshape<nsymlev)stop('symbol and symbolShape incompatible')
	}
	ncolor<-length(symbolColor)
	if(ncolor!=nsymlev){
		if(ncolor==1){
			symbolColor<-rep(symbolColor,nsymlev)
		}else if(ncolor<nsymlev)stop('symbol and symbolColor are incompatible')
	}

	binary<-match.arg(binary)
	BinResDir<-match.arg(BinResDir)

	if(binary=='BinRes' && missing(BinResLev))stop('BinResLev must be specified with the binary=BinRes option')
	if(any(count>1) && binary=='BinRes')stop('Binary responder computation cannot be specified with count>1')


		if(binary=='logit')transp<-function(x)plogis(x)
		if(binary=='probit')transp<-function(x)pnorm(x)

		### create binary y if cutlevel specified
		if(binary=='BinRes'){
		if(BinResDir=='>')y<- 1*(y>BinResLev)
		if(BinResDir=='<')y<- 1*(y<BinResLev)
		if(activeControl){
			if(BinResDir=='>')yac<- 1*(yac>BinResLev)
			if(BinResDir=='<')yac<- 1*(yac<BinResLev)
		}
	}


	doselev <- sort(unique(dose))
	ndose<-length(doselev)
	combG<-c(doselev,dgrid)
	nG<-length(combG)
	ndes<-tapply(count[!is.na(y)],dose, sum, na.rm = TRUE)
	Ngrid<-round(rep(mean(ndes),nG))
	Ngrid[1:ndose]<-ndes
	nsim<-length(sigma2)
	if(doselev[1]==0){npbo<-ndes[1]
	}else npbo<-round(mean(ndes))
	if(activeControl)nac<-sum(countac)

	cadj <- abs(qnorm((1 - clev)/2))
	low.clev <- (1 - clev)/2
	up.clev <- clev + (1 - clev)/2

	### placebo  and active control means
	### this pbo is used only if there is no pbo group
	### in the study.  ac only if there is active control
	if(binary=='logit' | binary=='probit' ){
		pbo<-modelFun(0,parm)
		pbo<-transp(pbo)
		pbopred<-rbinom(nsim,npbo,pbo)/npbo
		if(activeControl)acpred<-rbinom(nsim,nac,ac)/nac
	}else{
		pbo<-modelFun(0,parm)
		pbopred<-rnorm(nsim,pbo,sqrt(sigma2/npbo))
		if(activeControl)acpred<-rnorm(nsim,ac,sqrt(sigma2/nac))
	}
	if(!missing(BinResLev)){
		if(BinResDir=='>')pbo<-1-pnorm((BinResLev-pbo)/sqrt(sigma2))
		else pbo<-pnorm((BinResLev-pbo)/sqrt(sigma2))
		pbopred<-rbinom(nsim,npbo,pbo)/npbo
		if(activeControl){
			if(BinResDir=='>')ac<-1-pnorm((BinResLev-ac)/sqrt(sigma2))
			else ac<-pnorm((BinResLev-ac)/sqrt(sigma2))
			acpred<-rbinom(nsim,nac,ac)/nac
		}
	}

	Bmean <- NULL
	Bmed <- NULL
	lowL <- NULL
	upL <- NULL
	lowLpred <- NULL
	upLpred <- NULL

	Bdifmean <- NULL
	Bdifmed  <- NULL
	lowLdif <- NULL
	upLdif <- NULL
	lowLdifpred <- NULL
	upLdifpred <- NULL

	i<-0
	for (dlev in combG) {
		i<-i+1
		tmp <- modelFun(dlev, parm)
		if(binary=='logit' | binary=='probit')tmp<-transp(tmp)
		if(!missing(BinResLev)){
			if(BinResDir=='>')tmp<-1-pnorm((BinResLev-tmp)/sqrt(sigma2))
			else tmp<-pnorm((BinResLev-tmp)/sqrt(sigma2))
		}
		Bmean <- c(Bmean, mean(tmp))
		Bmed<-c(Bmed,median(tmp))
		lowL <- c(lowL, quantile(tmp, low.clev))
		upL <- c(upL, quantile(tmp, up.clev))
		if(binary=='logit' | binary=='probit' | !missing(BinResLev)){
			pred<-rbinom(nsim,Ngrid[i],tmp)/Ngrid[i]
		}else{
			pred<-rnorm(nsim,tmp,sqrt(sigma2/Ngrid[i]))
		}
		lowLpred <- c(lowLpred, quantile(pred, low.clev))
		upLpred <- c(upLpred, quantile(pred, up.clev))

		### assign pbo if part of design
		if(i==1 && doselev[1]==0){
			pbo<-tmp
			pbopred<-pred
		}

		if(activeControl){
			dif<-tmp-ac
			difpred<-pred-acpred
		}else{ 
			dif<- tmp-pbo
			difpred<-pred-pbopred
		}

		Bdifmean <- c(Bdifmean, mean(dif))
		Bdifmed <- c(Bdifmed, median(dif))
		lowLdif <- c(lowLdif, quantile(dif, low.clev))
		upLdif <- c(upLdif, quantile(dif, up.clev))
		lowLdifpred <- c(lowLdifpred, quantile(difpred, low.clev))
		upLdifpred <- c(upLdifpred, quantile(difpred, up.clev))
		
	}

	### split study doses from grid
	BmeanG<-Bmean[(ndose+1):nG];  Bmean<-Bmean[1:ndose]
	BmedG<-Bmed[(ndose+1):nG];  Bmed<-Bmed[1:ndose]
	lowLG<-lowL[(ndose+1):nG];  lowL<-lowL[1:ndose]
	upLG<-upL[(ndose+1):nG];  upL<-upL[1:ndose]
	lowLGpred<-lowLpred[(ndose+1):nG];  lowLpred<-lowLpred[1:ndose]
	upLGpred<-upLpred[(ndose+1):nG];  upLpred<-upLpred[1:ndose]

	BdifmeanG<-Bdifmean[(ndose+1):nG];  Bdifmean<-Bdifmean[1:ndose]
	BdifmedG<-Bdifmed[(ndose+1):nG];  Bdifmed<-Bdifmed[1:ndose]
	lowLdifG<-lowLdif[(ndose+1):nG];  lowLdif<-lowLdif[1:ndose]
	upLdifG<-upLdif[(ndose+1):nG];  upLdif<-upLdif[1:ndose]
	lowLdifGpred<-lowLdifpred[(ndose+1):nG];  lowLdifpred<-lowLdifpred[1:ndose]
	upLdifGpred<-upLdifpred[(ndose+1):nG];  upLdifpred<-upLdifpred[1:ndose]

	ym<-NULL
	symsub<-NULL
	dosesub<-NULL
	for(i in 1:ndose){
		for(j in 1:nsymlev){
			ind<-((dose==doselev[i]) & (symbol==symlev[j]))
			ng0<-sum(ind)
			if(ng0>0){
				ym<-c(ym,weighted.mean(y[ind],
						count[ind],na.rm=TRUE))
				symsub<-c(symsub,symlev[j])
				dosesub<-c(dosesub,doselev[i])
			}
		}
	}

	if(activeControl){
		ymac<-weighted.mean(yac,countac)
		ymdif<-ym-ymac
		acmean <- mean(ac)
		acmed <- median(ac)
		aclowL <- quantile(ac,low.clev) 
		acupL <-  quantile(ac,up.clev)
		aclowLpred <- quantile(acpred,low.clev)
		acupLpred<- quantile(acpred,up.clev)
		AC<-c(acmean,acmed,aclowL,acupL,aclowLpred,
			acupLpred,ymac,dac)
		names(AC)<-c('acmean','acmed','aclowL','acupL','aclowLpred','acupLpred','yac','dac')
	}else{
		AC<-NULL
		ymdif<-ym - ym[1]
	}

	pairwise<-cbind(ym, ymdif)
	rownames(pairwise)<-dosesub
	modelABS<-cbind(Bmean, Bmed,lowL,upL, lowLpred, upLpred)
	rownames(modelABS)<-doselev
	modelDIF<-cbind(Bdifmean, Bdifmed,lowLdif, 
					upLdif, lowLdifpred, upLdifpred)
	rownames(modelDIF)<-doselev
	modelABSG<-cbind(BmeanG, BmedG, lowLG, upLG, lowLGpred, upLGpred)
	rownames(modelABSG)<-dgrid
	modelDIFG<-cbind(BdifmeanG, BdifmedG, lowLdifG, upLdifG,
					lowLdifGpred,upLdifGpred)
	rownames(modelDIFG)<-dgrid

    plotBobj<-structure(list(pairwise=pairwise,dose=dosesub,symbol=factor(symsub),nolegend=nolegend, 
				modelABS=modelABS,
				modelDIF=modelDIF, 
				modelABSG=modelABSG,
				modelDIFG=modelDIFG,
				activeControl=activeControl,AC=AC,
				dgrid=dgrid,clev=clev 
				),class='plotB')

	if(makePlot){
		print(plot(plotBobj,
		plotDif=plotDif,plotMed=plotMed,plotResid=plotResid, 
		predict=predict,
		xlim=xlim, ylim=ylim, xlab = xlab, ylab = ylab, 
		labac=labac,shapeac=shapeac,colac=colac,
		symbol=symbol,symbolLabel=symbolLabel,symbolShape=symbolShape,
		symbolColor=symbolColor,symbolSize=symbolSize))
	}

	return(plotBobj)
}


