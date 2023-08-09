"print.emaxPrior"<- 
	function(x, doc=FALSE, diffuse=NULL, file="", 
					 modType=c('4','3'), docType=c('sap','protocol'), ... )
	{
		
		if(!doc){
			class(x)<-NULL
			print(x)
			return(invisible(x))
		}
		
		if(!is.logical(diffuse))stop('diffuse should be TRUE/FALSE')
		
	
		docType<-match.arg(docType,c('sap','protocol'))
		prot<- (docType=='protocol')
		modType<-match.arg(modType,c('4','3'))
		binary<-x$binary
		effDF<-x$effDF
		parmDF<-x$parmDF
		epsca<-x$epsca
		difTargetsca<-x$difTargetsca
		p50<-x$p50
		sigmalow<-x$sigmalow
		sigmaup<-x$sigmaup
		loged50sca<-x$loged50sca
		loglamsca<-x$loglamsca
		parmCor<-x$parmCor
		
		if(!x$default)stop(paste0('No prior description printed if non-',
															'default values for the meta-analysis parameters are specified'))
		
		
		if(length(effDF)>1){
			e0df<-effDF[1]
			diftdf<-effDF[2]
		} else{e0df<-effDF; diftdf<-effDF}
		
		if(diffuse && binary){
			if(epsca<2.0 || difTargetsca<2.0)warning(paste0('epsca and ',
																											'difTargetSca not consistent with diffuse prior specification'))
		}
		
		baseline<-(!is.null(x$basemu))
		if(baseline){
			basemu<-x$basemu
			basevar<-x$basevar
			lenb<-length(basemu)
			if(lenb>1)multiv<-'multivariate ' else multiv<-''
			bchar<-paste0(round(basemu,3),collapse=',')
			if(lenb==1){
				btxt1<-'the baseline covariate parameter'
			}else btxt1<-paste0('the ',lenb,' baseline covariate parameters')
			if(lenb>1){
				btxt2<-paste0('(',bchar,
											'), and the prior variance-covariance matrix is ',
											'(user must specify, but usually diffuse and independent)'  
				)
			}else{
				btxt2<-paste0(bchar,
											', and the prior variance is ',basevar)
			}	
		}
		
		
		textout<- paste0('The specification of the Emax model and the ', 
										 'prior distribution for some of its parameters are ',
										 'based on three meta-analyses of clinical dose response ', 
										 'that include more than 200 compounds (Thomas, 2014; ', 
										 'Thomas & Roy, 2017; Wu, 2017).  ' 
		)
		if(modType==4){
			textout<-paste0(textout,
											'The model parameters are the placebo response (E0), ',
											'the ED50, the drug effect at the targeted dose ',
											'denoted by difTarget, and the Hill (slope) parameter, ',
											'lambda.  '
			)
		}else{
			textout<-paste0(textout,
											'The model parameters are the placebo response (E0), ',
											'the ED50, and the drug effect at the targeted dose ',
											'denoted by difTarget. '
			)		
		}
		
		if(!prot){
			if(binary){
				if(diffuse){
					textout<-paste0(textout,
													'\n\nThe logit of the placebo response will have a t-distribution ',
													'with ',e0df,' degrees of freedom (df), centered at logit(',
													round(plogis(x$epmu),3),').  The scale parameter for the t-distribution will be ',
													epsca,', which yields a diffuse prior distribution.  A ',
													't-distribution with ', diftdf,' df is also used for the effect ',
													'parameter, difTarget, which is also specified on the ',
													'logit scale.  It is centered at ',
													x$difTargetmu, 
													', with a scale parameter of ',x$difTargetsca,', yielding ',
													'a diffuse prior distribution.  ' 
					)	
				}else{
					textout<-paste0(textout,
													'\n\nThe logit of the placebo response will have a t-distribution ',
													'with ',e0df,' degrees of freedom (df), centered at logit(',
													round(plogis(x$epmu),3),').  The scale parameter for the t-distribution will be ',
													epsca,', which yields an informative prior distribution', 
													' supported by (insert text here).  A ',
													't-distribution with ', diftdf,' df is also used for the effect ',
													'parameter, difTarget, which is also specified on the ',
													'logit scale.  It is centered at ',
													x$difTargetmu, 
													', with a scale parameter of ',x$difTargetsca,', yielding ',
													'an informative prior distribution supported by (insert here).  ' 
					)					
				}
				
			}else{
				if(diffuse){
					textout<-paste0(textout,
													'\n\nThe placebo response will have a t-distribution ',
													'with ',e0df,' degrees of freedom (df), centered at ',
													round(plogis(x$epmu),3),'.  The scale parameter ',
													'for the t-distribution will be ',
													epsca,', which yields a diffuse prior distribution.  A ',
													't-distribution with ', diftdf,' df is also used for the effect ',
													'parameter, difTarget, which is ',
													'centered at ', x$difTargetmu, 
													', with a scale parameter of ',x$difTargetsca,', yielding ',
													'a diffuse prior distribution. The specified scale parameters ',
													'produce diffuse prior distributions relative to the ',
													'variability in the treatment group responses because ',
													'(give justification here, typically scale exceeds ',
													'individual patient variability).  '
					)	
				}else{
					textout<-paste0(textout,
													'\n\nThe placebo response will have a t-distribution ',
													'with ',e0df,' degrees of freedom (df), centered at ',
													round(plogis(x$epmu),3),'.  The scale parameter for the t-distribution will be ',
													epsca,', which yields an informative prior distribution', 
													' supported by (insert text here).  A ',
													't-distribution with ', diftdf,' df is also used for the effect ',
													'parameter, difTarget, which is centered at ',
													x$difTargetmu, 
													'), with a scale parameter of ',x$difTargetsca,', yielding ',
													'an informative prior distribution supported by (insert here).  ' 
					)					
				}
				textout<-paste0(textout,
												'The residual standard deviation has a uniform ',
												'distribution with bounds of (',
												sigmalow,',',sigmaup,').  '
				)
			}
			if(baseline){
				textout<-paste0(textout,
												'The prior distribution for ',btxt1,
												' is a ',multiv, 'normal distribution.  The prior mean is ',
												btxt2,
												'.  The covariate regression parameters are a priori independent ',
												'of all other parameters.'
				)			
			}	
			if(modType=='4'){
				textout<-paste0(textout,
												'\n\nThe prior distribution for the (log(ED50),log(lambda)) is ',
												'a bivariate t-distribution with ',parmDF,' df, which was ', 
												'derived from the meta-analysis of past dose response studies.  ',  
												'The current projected value for the ED50 is p50=',p50,
												' mg, which is based on (insert the basis for the projection ',
												'here).  The prior distribution for ',
												'log(ED50) is centered at the log(P50), with substantial ',
												'uncertainty specified by the scale parameter of ',
												loged50sca,', which was derived from the meta-analysis.  The 80% prior prediction interval for the ',
												'ED50 is approximately (P50/10,10P50).   The log(lambda) prior ',
												'distribution is centered at 0, with variation determined ',
												'by the scale parameter of ',loglamsca,
												', which was derived from the meta-analysis.  The 80% prior ',
												'prediction interval for lambda is approximately (0.5,2).  ',
												'The prior correlation between the parameters is ',
												parmCor,', which is also determined from the meta-analysis.'
				)
			}else{
				textout<-paste0(textout,
												'\n\nThe prior distribution for the log(ED50) is ',
												'a t-distribution with ',parmDF,' df, which was ', 
												'derived from the meta-analysis of past dose response studies.  ',  
												'The current projected value for the ED50 is p50=',p50,
												' mg, which is based on (insert the basis for the projection ',
												'here).  The prior distribution for ',
												'log(ED50) is centered at the log(P50), with substantial ',
												'uncertainty specified by a scale parameter of ',
												loged50sca,', which was derived from the meta-analysis. ', 
												'The 80% prior prediction interval for the ',
												'ED50 is roughly (P50/10,10P50).  '
				)
			}	
		}else{
			if(binary){
				if(diffuse){
					textout<-paste0(textout,
													'\n\nThe logit of the placebo response will have a t-distribution ',
													'with ',e0df,' degrees of freedom (df) centered at logit(',
													round(plogis(x$epmu),3),') with a scale parameter specified to ',
													'yield a diffuse prior distribution.  A ',
													't-distribution with ', diftdf,' df is also used for the effect ',
													'parameter, difTarget, which is also specified on the ',
													'logit scale with a scale parameter specified to yield a ',
													'diffuse prior distribution.  ' 
					)	
				}else{
					textout<-paste0(textout,
													'\n\nThe logit of the placebo response will have a t-distribution ',
													'with ',e0df,' degrees of freedom (df) centered at logit(',
													round(plogis(x$epmu),3),').  The scale parameter for the ',
													't-distribution is ',epsca,', which yields an informative prior ', 
													'distribution supported by (insert text here).  A ',
													't-distribution with ', diftdf,' df is also used for the effect ',
													'parameter, difTarget, which is specified on the ',
													'logit scale.  It is centered at ',
													x$difTargetmu, 
													', with a scale parameter of ',x$difTargetsca,', yielding ',
													'an informative prior distribution supported by (insert here).  ' 
					)					
				}
				
			}else{
				if(diffuse){
					textout<-paste0(textout,
													'\n\nThe placebo response will have a t-distribution ',
													'with ',e0df,' degrees of freedom (df), centered at ',
													round(plogis(x$epmu),3),'.  The scale parameter ',
													'will be chosen to yield a diffuse prior distribution.  A ',
													't-distribution with ', diftdf,' df is also used for the effect ',
													'parameter, difTarget, which is ',
													'centered at ', x$difTargetmu, 
													', with a scale parameter chosen to yield ',
													'a diffuse prior distribution.  ' 
					)	
				}else{
					textout<-paste0(textout,
													'\n\nThe placebo response will have a t-distribution ',
													'with ',e0df,' degrees of freedom (df) centered at ',
													round(plogis(x$epmu),3),'.  The scale parameter for the ',
													't-distribution is ',epsca,', which yields an informative ', 
													'prior distribution supported by (insert text here).  A ',
													't-distribution with ', diftdf,' df is also used for the effect ',
													'parameter, difTarget, which is centered at ',
													x$difTargetmu, 
													'), with a scale parameter of ',x$difTargetsca,', yielding ',
													'an informative prior distribution supported by (insert here).  ' 
					)					
				}
				textout<-paste0(textout,
												'The residual standard deviation has a uniform ',
												'distribution with bounds of (',
												sigmalow,',',sigmaup,').  '
				)
			}
			if(baseline){
				textout<-paste0(textout,
												'The prior distribution for ',btxt1,
												' is a ',multiv, 'normal distribution.  The prior mean is ',
												btxt2,
												'.  The covariate regression parameters are a priori independent ',
												'of all other parameters.'
				)			
			}	
			if(modType=='4'){
				textout<-paste0(textout,
												'\n\nThe prior distribution for the (log(ED50),log(lambda)) is ',
												'a bivariate t-distribution with ',parmDF,' df, which was ', 
												'derived from the meta-analysis of past dose response studies.  ',  
												'The current projected value for the ED50 is p50=',p50,
												' mg, which is based on (insert the basis for the projection ',
												'here).  The prior distribution for ',
												'log(ED50) is centered at the log(P50), with an approximate 80% ', 
												'prior prediction interval ',
												'for the ED50 of (P50/10,10P50).   The log(lambda) prior ',
												'distribution is centered at 0, with an approximate 80% prior ',
												'prediction interval for lambda of (0.5,2).  ',
												'The prior correlation between the parameters is ',
												parmCor,', which is also determined from the meta-analysis.'
				)
			}else{
				textout<-paste0(textout,
												'\n\nThe prior distribution for the log(ED50) is ',
												'a t-distribution with ',parmDF,' df, which was ', 
												'derived from the meta-analysis of past dose response studies.  ',  
												'The current projected value for the ED50 is p50=',p50,
												' mg, which is based on (insert the basis for the projection ',
												'here).  The prior distribution for ',
												'log(ED50) is centered at the log(P50), with an ',
												'an approximate 80% prior prediction interval for the ',
												'ED50 of (P50/10,10P50).  '
				)
			}	
		}
		
		
		
		if(prot){
			textout<-paste0(textout,	
											'\n\nThe prior distributions may be updated in the ',
											'statistical analysis plan if updated information ',
											'becomes available.  ',
											'For all of the planned Bayesian analyses, graphical ', 
											'displays and model checking diagnostics will be ', 
											'evaluated.  '
			)
		}else{
			textout<-paste0(textout,	
											'\n\n', 
											'The Bayesian model will be evaluated using ',
											'Markov Chain Monte ',
											'Carlo (MCMC) simulation.  The fitted model ',
											'and posterior prediction intervals will be ',
											'plotted along with the dose group ',
											'sample means/proportions.  ',
											'Convergence will be checked ',
											'using trace plots, ',
											'auto-correlation plots, and Gelman-Rubin ',
											'divergence statistics. '
			)	
		}
		
		textout<-paste0(textout,
										'\n\n\nThomas N, Roy D. Analysis of clinical dose-response in ',
										'small-molecule drug development: 2009-2014. Statistics in ', 
										'Biopharmaceutical Research 2017;9(2):137-46.\n\nThomas N, ',
										'Sweeney K, Somayaji V. Meta-analysis of clinical dose ',
										'response in a large drug development portfolio. Statistics ',
										'in Biopharmaceutical Research 2014;6(4):302-17.\n\n',
										'Wu J, Banerjee A, Jin B et al. Clinical dose-response for ',
										'a broad set of biological products: A model-based ',
										'meta-analysis. Statistical Methods in Medical Research ',
										'2018;27(9):2694-2721.\n\n'
		)
		cat(textout,file=file)
	}