"print.emaxPrior"<- 
  function(x, doc=FALSE, diffuse=c(pbo=TRUE,eff=TRUE), 
           file=paste('emaxPrior',format(Sys.time(), 
                                               "%b%d%H%M%S"),'.docx',sep=''),
           modType=c('4','3'), docType=c('sap','protocol'),
           font.family = "Courier New", rounding=3, ...) 
  {
    
    if(!doc){
      class(x)<-NULL
      print(x)
      return(invisible(x))
    }
    
    if(length(diffuse)!=2)stop('diffuse must be specified for pbo and eff')
    if(any(!is.logical(diffuse)))stop('diffuse should have TRUE/FALSE values')
    
    
    docType<-match.arg(docType,c('sap','protocol'))
    prot<- (docType=='protocol')
    modType<-match.arg(modType,c('4','3'))
    binary<-x$binary
    effDF<-x$effDF
    parmDF<-x$parmDF
    epmu <- x$epmu                     #update
    mixP <- x$mixP                     #update
    mu_ep <- x$mu_ep                   #update
    sd_ep <- x$sd_ep                   #update
    w_ep <- x$w_ep                     #update
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
    
    #Update
    if(diffuse[1] && mixP > 1)warning(paste0('Specifying a ', 
                                 'diffuse prior for E0 using a multi-component',
                                 ' mixture (mixP > 1) is unusual'))
    
    if(any(diffuse) && binary){
      options(warn=1)
      if(mixP > 0 && diffuse[2]){
        if(difTargetsca<2.0)warning(paste0('difTargetSca not consistent ',
        'with diffuse prior specification'))
      } else if(mixP==0){
        if(diffuse[1]){
          if(epsca<2.0 || difTargetsca<2.0)warning(paste0('epsca ',
                                    'is not consistent with diffuse',
                                    ' prior specification'))
        }
        if(diffuse[2]){
          if(difTargetsca<2.0)warning(paste0(
                                'difTargetSca is not consistent with diffuse',
                                ' prior specification'))                 
        }
      }
    }
    
    baseline<-(!is.null(x$basemu))
    if(baseline && diffuse[1])warning(paste0('Combining covariate ',
                                'adjustment and an informative placebo ',
                                'prior distribution based on historical ',
                                'data is not recommended'))
    if(baseline){
      basemu<-x$basemu
      basevar<-x$basevar
      lenb<-length(basemu)
      if(lenb>1)multiv<-'multivariate ' else multiv<-''
      bchar<-paste0(round(basemu,rounding),collapse=',')
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
        if(diffuse[1]){
          if(mixP >0){
            print("sap/binary/diffuse/mixp")
            
            textout <- paste0(textout,
                            '\n\nThe prior distribution for placebo response,',
                            ' will be a mixture of a normal distributions ',
                            'with means of (',paste(mu_ep,collapse=','),'), ', 
                            'and scale parameters (',
                            paste(sd_ep,collapse=','),'), ',
                            'which yields a diffuse distribution. ')
            if(diffuse[2]){
              textout <- paste0(textout,
                               'A t-distribution with ', diftdf,' df is used for the effect ',
                               'parameter, difTarget, which is also specified on the ',
                               'logit scale.  It is centered at ',
                               round(x$difTargetmu,rounding), 
                               ', with a scale parameter of ',
                               round(x$difTargetsca,rounding),', yielding ',
                               'a diffuse prior distribution. (! reference)  ' 
                          )
            }else{
              textout <- paste0(textout,
                              'A t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget, which is also specified on the ',
                              'logit scale.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',round(x$difTargetsca,rounding),
                              ', yielding ',
                              'an informative prior distribution supported by (insert here).  ')
           } 
            
          }else{
            print("sap/binary/diffuse/t")
            textout<-paste0(textout,
                          '\n\nThe logit of the placebo response will have a t-distribution ',
                          'with ',e0df,' degrees of freedom (df), centered at logit(',
                          round(plogis(epmu),rounding),').  The scale parameter for the t-distribution will be ',
                          epsca,', which yields a diffuse prior distribution (!reference).  '
                          )	
            if(diffuse[2]){
              textout <- paste0(textout,
                               'A t-distribution with ', diftdf,' df is used for the effect ',
                               'parameter, difTarget, which is also specified on the ',
                               'logit scale.  It is centered at ',
                               round(x$difTargetmu,rounding), 
                               ', with a scale parameter of ',
                               round(x$difTargetsca,rounding),', also yielding ',
                               'a diffuse prior distribution. (! reference)  ' 
                          )
            }else{
              textout <- paste0(textout,
                              'A t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget, which is also specified on the ',
                              'logit scale.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',round(x$difTargetsca,rounding),', yielding ',
                              'an informative prior distribution supported by (insert here).  ')
           } 
          }
        }else{
          if(mixP>0){
            print("sap/binary/informative/mixp")
            textout <- paste0(textout,
                              '\n\nThe logit of the placebo response will be approximated by a mixture of ',
                              mixP,' normal distributions.  ',
                              'It consists of a MAP prior derived from historical data. ',
                              'The mixture weights are (',paste0(w_ep, collapse = ", "),')',
                              ' and the mean and SD of each component are ',
                              '(',paste0(round(mu_ep,rounding),collapse=", "),')', ' and ',
                              '(',paste0( round(sd_ep,rounding),collapse=", "),').',
                              '  The historical information is in (insert text here). '
                              )					
            
            if(diffuse[2]){
              textout <- paste0(textout,
                               'At-distribution with ', diftdf,' df is used for the effect ',
                               'parameter, difTarget, which is also specified on the ',
                               'logit scale.  It is centered at ',
                               round(x$difTargetmu,rounding), 
                               ', with a scale parameter of ',
                               round(x$difTargetsca,rounding),', yielding ',
                               'a diffuse prior distribution. (! reference)  ' 
                          )
            }else{
              textout <- paste0(textout,
                              'A t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget, which is also specified on the ',
                              'logit scale.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',round(x$difTargetsca,rounding),
                              ', yielding ',
                              'an informative prior distribution supported by (insert here).  ')
           } 
          } else{
            print("sap/binary/informative/t")
            textout <- paste0(textout,
                          '\n\nThe logit of the placebo response will have a t-distribution ',
                          'with ',e0df,' degrees of freedom (df), centered at logit(',
                          round(plogis(epmu),rounding),').  The scale parameter for the t-distribution will be ',
                          round(epsca,rounding),', which yields an informative prior distribution', 
                          ' supported by (insert text here).  '
                          )					
              if(diffuse[2]){
                textout <- paste0(textout,
                               'At-distribution with ', diftdf,' df is used for the effect ',
                               'parameter, difTarget, which is also specified on the ',
                               'logit scale.  It is centered at ',
                               round(x$difTargetmu,rounding), 
                               ', with a scale parameter of ',
                               round(x$difTargetsca,rounding),', yielding ',
                               'a diffuse prior distribution. (! reference)  ' 
                          )
              }else{
                textout <- paste0(textout,
                              'A t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget, which is also specified on the ',
                              'logit scale.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'an informative prior distribution supported by (insert here).  ')
           } 
          }
        }
      }else{
        if(diffuse[1]){
          if(mixP>0){
            print("sap/continuous/diffuse/mixp")
            textout <- paste0(textout,
                            '\n\nThe prior distribution for placebo response ',
                            'will be a mixture of normal distributions ',
                            'with means of (',paste(mu_ep,collapse=','),'), ', 
                            'and scale parameters (',
                            paste(sd_ep,collapse=','),'), ',
                            'which yields a diffuse distribution ',
                            'with variability ',  
                            'large relative to the variability in ',
                            'response (!references, and marginal mixture SD).')
            if(diffuse[2]){
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'a diffuse prior distribution (typical:  because ',
                              'the prior scale parameter is large compared to ',
                              'the response variability.)'
                              )
            }else{
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,
                              ' df is used for the effect ',
                              'parameter, difTarget.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'an informative prior distribution ',
                              'supported by (insert here).  '
                          )
            }
          }else{
            print("sap/continuous/diffuse/t")
            textout<-paste0(textout,
                            '\n\nThe placebo response will have a t-distribution ',
                            'with ',e0df,' degrees of freedom (df) centered at ',
                            round(epmu,rounding),'.  The scale parameter ',
                            'for the t-distribution will be ',
                            round(epsca,rounding),', which yields a diffuse ',
                            'prior distribution  (typical:  because ',
                            'the prior scale parameter is large compared to ',
                            'the response variability.)'
                           )
             if(diffuse[2]){
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'a diffuse prior distribution (typical:  because ',
                              'the prior scale parameter is large compared to ',
                              'the response variability.)  '
                              )
             }else{
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,
                              ' df is used for the effect ',
                              'parameter, difTarget.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'an informative prior distribution ',
                              'supported by (insert here).  '
                          )
            }                        	
          }
          
        }else{
          if(mixP>0){
            print("sap/continuous/informative/mixp")
            textout <- paste0(textout,
                              '\n\nThe prior distribution of placebo response ',
                              'is a MAP (meta-analytic predictive) prior ',
                              'consisting of a mixture of ',
                              mixP,' normal distributions.  ',
                              'The mixture weights are (',
                              paste0(w_ep, collapse = ", "),'), ',
                              'and the means and SDs of the component are ',
                              '(',paste0(round(mu_ep,rounding), collapse = ", "),')',', and ',
                              '(',paste0(round(sd_ep,rounding), collapse = ", "),')',  '.  ',
                              'The MAP prior was derived from (text, references).  '
                        )
             if(diffuse[2]){
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'a diffuse prior distribution (typical:  because ',
                              'the prior scale parameter is large compared to ',
                              'the response variability.)  '
                              )
             }else{
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,
                              ' df is used for the effect ',
                              'parameter, difTarget.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'an informative prior distribution ',
                              'supported by (insert here).  '
                          )
            }                        	
          }else{
            print("sap/continuous/informative/t")
            textout<-paste0(textout,
                          '\n\nThe placebo response will have a t-distribution ',
                          'with ',e0df,' degrees of freedom (df), centered at ',
                          round(epmu,rounding),'.  The scale parameter for the ',
                          't-distribution will be ',
                          round(epsca,rounding),
                          ', which yields an informative prior distribution', 
                          ' supported by (insert text here).  '
                         )
             if(diffuse[2]){
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'a diffuse prior distribution (typical:  because ',
                              'the prior scale parameter is large compared to ',
                              'the response variability.)  '
                              )
             }else{
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,
                              ' df is used for the effect ',
                              'parameter, difTarget.  It is centered at ',
                              round(x$difTargetmu,rounding), 
                              ', with a scale parameter of ',
                              round(x$difTargetsca,rounding),', yielding ',
                              'an informative prior distribution ',
                              'supported by (insert here).  '
                          )
            }                        	                        
          }
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
                        loged50sca,', which was derived from the meta-analysis. ',  
                        'An approximate 80% prior prediction interval for the ',
                        'ED50 is approximately (P50/10,10P50).   The log(lambda) prior ',
                        'distribution is centered at 0, with variation determined ',
                        'by the scale parameter of ',loglamsca,
                        ', which was derived from the meta-analysis. ', 
                        'An approximate 80% prior ',
                        'prediction interval for lambda is approximately (0.5,2).  ',
                        'The prior correlation between the parameters is ',
                        parmCor,', which was also determined from the meta-analysis.'
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
                        'An approximate 80% prior prediction interval for the ',
                        'ED50 is roughly (P50/10,10P50).  '
        )
      }	
    } else {##start of protocol
      if(binary){
        if(diffuse[1]){
          if(mixP>0){
            print("protocol/binary/diffuse/mixp")
            textout <- paste0(textout,
                              '\n\nThe prior distribution for placebo response ',
                              'will be a mixture of a normal distributions ',
                              'with means and SDs selected to ',
                              'yield a diffuse distribution. ',
                              'The prior distribution is specified on the ',
                              'logit scale.  ')
            if(diffuse[2]){
              textout <- paste0(textout,
                                'A t-distribution with ', diftdf,' df is used for the effect ',
                                'parameter, difTarget, which is also specified on the ',
                                'logit scale.  The mean and scale parameters ',
                                'will be selected to yield ',
                                'a diffuse prior distribution.  ' 
              )
            }else{
              textout <- paste0(textout,
                                'A t-distribution with ', diftdf,' df is used for the effect ',
                                'parameter, difTarget, which is also specified on the ',
                                'logit scale.  An informative prior distribution ',
                                'is planned for this parameter.  The empirical basis for the ',
                                'informative prior will be described in the SAP.  '
                                )
            } 
          } else {
            print("protocol/binary/diffuse/t")
            textout<-paste0(textout,
                          '\n\nThe logit of the placebo response will have a t-distribution ',
                          'with ',e0df,' degrees of freedom (df) centered at logit(',
                          round(plogis(epmu),rounding),') with a scale parameter specified to ',
                          'yield a diffuse prior distribution.  ' 
                            )
             if(diffuse[2]){
              textout <- paste0(textout,
                                'A t-distribution with ', diftdf,' df is used for the effect ',
                                'parameter, difTarget, which is also specified on the ',
                                'logit scale.  The prior mean will be ',
                                round(x$difTargetmu,rounding),
                                ' and the scale parameter ',
                                'will be selected to yield ',
                                'a diffuse prior distribution.  ' 
              )
            }else{
              textout <- paste0(textout,
                                'A t-distribution with ', diftdf,' df is used for the effect ',
                                'parameter, difTarget, which is also specified on the ',
                                'logit scale.  An informative prior distribution ',
                                'is planned for this parameter.  The empirical basis for the ',
                                'informative prior will be described in the SAP.  '
                                )
            }
          }
        } else {
          if(mixP>0){
            print("protocol/binary/informative/mixp")
            textout <- paste0(textout,
                       '\n\nThe logit of the placebo response will ', 
                        'be approximated by a ',
                        mixP,'-component normal mixture.  ',
                          "The derivation of the informative prior ", 
                          "will be described in the SAP.  "
                          )
             if(diffuse[2]){
              textout <- paste0(textout,
                                'A t-distribution with ', diftdf,' df is used for the effect ',
                                'parameter, difTarget, which is also specified on the ',
                                'logit scale.  The prior mean will be ',
                                round(x$difTargetmu,rounding),
                                ' and the scale parameter ',
                                'will be selected to yield ',
                                'a diffuse prior distribution.  ' 
                          )
            }else{
              textout <- paste0(textout,
                                'A t-distribution with ', diftdf,' df is used for the effect ',
                                'parameter, difTarget, which is also specified on the ',
                                'logit scale.  An informative prior distribution ',
                                'is planned for difTarget.  The empirical basis for the ',
                                'informative prior will be described in the SAP.  '
                                )
            }
            
          } else {
            print("protocol/binary/informative/t")
            textout<-paste0(textout,
                          '\n\nThe logit of the placebo response will have a t-distribution ',
                          'with ',e0df,' degrees of freedom (df). ',
                          'An informative prior distribution ',
                          'is planned.  The empirical basis for the ',
                           'informative prior will be described in the SAP.  '
                          )
            if(diffuse[2]){
              textout <- paste0(textout,
                                'A t-distribution with ', diftdf,' df is used for the effect ',
                                'parameter, difTarget, which is also specified on the ',
                                'logit scale.  The prior mean will be ',
                                round(x$difTargetmu,rounding),
                                ' and the scale parameter ',
                                'will be selected to yield ',
                                'a diffuse prior distribution.  ' 
              )
            }else{
              textout <- paste0(textout,
                                'A t-distribution with ', diftdf,' df is used for the effect ',
                                'parameter, difTarget, which is also specified on the ',
                                'logit scale.  An informative prior distribution ',
                                'is planned for this parameter.  The empirical basis for the ',
                                'informative prior will be described in the SAP.  '
                                )
            }
          }
        }
        
      }else{
        if(diffuse[1]){
          if(mixP>0){
            print("protocol/continuous/diffuse/mixp")
              textout <- paste0(textout,
                            '\n\nThe prior distribution for placebo response ',
                            'will be a mixture of normal distributions ',
                            'which yields a diffuse distribution ',
                            'with variability ',  
                            'large relative to the variability in ',
                            'response. Details of this prior distribution ',
                            'will be in the SAP.  ')
            if(diffuse[2]){
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget.  Its scale parameter will ',
                              'be large enough to create a diffuse distribution.  ',
                              'Details of these prior distributions ',
                              'will be in the SAP.  '
                              )
            }else{
              textout <- paste0(textout,
                              '\n\nAn informative t-distribution with ', diftdf,
                              ' df is used for the effect ',
                              'parameter, difTarget.  Details of this ',
                              'distribution and its derivation will be ',
                              'in the SAP.  '
                          )
            } 
          } else {
            print("protocol/continuous/diffuse/t")
            textout<-paste0(textout,
                            '\n\nThe placebo response will have a t-distribution ',
                            'with ',e0df,' degrees of freedom (df).  ',
                            'A large scale parameter ',
                            'for the t-distribution will be ',
                            'specified that yields a diffuse ',
                            'prior distribution.  '
                           )
             if(diffuse[2]){
              textout <- paste0(textout,
                              'A t-distribution with ', diftdf,
                              ' df is used for the effect ',
                              'parameter, difTarget.  ', 
                              'A large scale parameter of ',
                              'will also be specified to yield ',
                              'a diffuse prior distribution.  ',
                              'Details of these ',
                              'prior distributions will be in the SAP. '
                              )
              }else{
              textout <- paste0(textout,
                              'An informative t-distribution with ', diftdf,
                              ' df is used for the effect ',
                              'parameter, difTarget. Details of these ',
                              'distributions and their derivation will be',
                              'in the SAP.  ',
                          )
              }                        	
          }
        }else{
          if(mixP>0){
            print("protocol/continuous/informative/mixp")
              textout <- paste0(textout,
                              '\n\nThe prior distribution of placebo response ',
                              'is a MAP (meta-analytic predictive) prior ',
                              'consisting of a mixture of ',
                              'normal distributions.  ',
                              'Details of the MAP prior and its derivation ',
                              'will be in the SAP'
                        )
             if(diffuse[2]){
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget.  The scale parameter of ',
                              'the t-distribution will be large compared to ',
                              'the response variability yielding a ',
                              'diffuse prior distribution.  '
                              )
             }else{
              textout <- paste0(textout,
                              '\n\nAn informative t-distribution with ', diftdf,
                              ' df is used for the effect ',
                              'parameter, difTarget.  Details of the ',
                              'prior distribution and the empirical ',
                              'support for its use will be in the SAP.  '
                          )
            }                        	 
          } else {
            print("protocol/continuous/informative/t")
            textout <- paste0(textout,
                              '\n\nThe placebo response will have an ',
                              'informative t-distribution ',
                              'with ',e0df,' degrees of freedom (df).  ',
                              'The details of this distribution and ',
                              'the empirical support for an informative ',
                              'prior distribution will be in the SAP. '
                       )
            if(diffuse[2]){
              textout <- paste0(textout,
                              '\n\nA t-distribution with ', diftdf,' df is used for the effect ',
                              'parameter, difTarget.  The scale parameter of ',
                              'the t-distribution will be large compared to ',
                              'the response variability yielding a ',
                              'diffuse prior distribution.  '
                              )
             }else{
              textout <- paste0(textout,
                              '\n\nAn informative t-distribution with ', diftdf,
                              ' df is also used for the effect ',
                              'parameter, difTarget.  Details of the ',
                              'prior distribution and the empirical ',
                              'support for its use will be in the SAP.  '
                          )
            }                        	                              
          }
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
                        parmCor,', which was also determined from the meta-analysis.'
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
    
    out.text <- capture.output(cat(textout))
    # Create a new Word document
    doc.out <- read_docx()
    
    # Add each line separately to preserve spacing and blank lines
    for (line in out.text) {
      if (line == "") {
        # Insert a blank paragraph for empty lines
        doc.out <- body_add_par(doc.out, value = " ", style = "Normal")
      } else {
        # Insert text with monospace font for alignment
        doc.out <- body_add_fpar(
          doc.out,
          fpar(
            ftext(line, prop = fp_text(font.family = font.family))
          )
        )
      }
    }
    
    # Save the document
    if(!is.null(file)){
      print(doc.out, target = file)
      message("Word file saved to: ", normalizePath(file))
    }
    return(invisible(textout)) 
  }

