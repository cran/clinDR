data {
    int<lower=0> N; // number of patients/groups 
    int<lower=0> nprot; // number of protocols-strata
		int<lower=0> protv[N];
		int<lower=0,upper=1> cont;  //indicator of data type
		int<lower=0,upper=1> sigmoid;  //indicator of model type
		int<lower=0,upper=1> gp;
		int<lower=0,upper=1> intercept;
		int<lower=0> nbase;
		int<lower=0> dimFit;
		
		//continuous data
    vector[N] yv; //  outcome variable
		vector<lower=0>[N] nv; // number patients per group; real for arithmetic
		//binary data
    int yvb[N]; //  outcome variable
		int<lower=0> nvb[N]; // number patients per group
		
    vector<lower=0>[N] dv;  
    matrix[nbase ? N : 0,nbase ? nbase : 0] xbase;
    
		real<lower=0> df2;   // saturated df divided by 2
		real<lower=0> ssy;
    real epmu;
		real<lower=0> epsca;
    real difTargetmu;
		real<lower=0> difTargetsca;
		real<lower=0> dTarget;
		real<lower=0> sigmalow;
		real<lower=0> sigmaup;
		real<lower=0> p50;
    real loged50mu;
    real<lower=0> loged50sca;
    real<lower=0> e0DF;
    real<lower=0> diftDF;
    real<lower=0> parmDF;
    real<lower=0> loglammu;
    real<lower=0> loglamsca;
		real <lower=-1,upper=1>parmCor;
		real lowled50;
		real highled50;
		real lowllam;
		real highllam;
		vector[nbase ? nbase : 0] basemu;
		cov_matrix[nbase ? nbase : 1] basevar;	
    cov_matrix[dimFit ? dimFit : 1] vcest;
}

transformed data{
		cholesky_factor_cov[nbase ? nbase : 1] cbvar;
		cholesky_factor_cov[dimFit ? dimFit : 1] cholvc;

		cbvar=cholesky_decompose(basevar);
		cholvc=cholesky_decompose(vcest);
}


parameters{
    vector[intercept ? nprot : 0] e0;
    real difTarget;
		real<lower=lowled50,upper=highled50> loged50;      //implicitly normed by p50
		real<lower=lowllam,upper=highllam> llam;
		vector<lower=sigmalow,upper=sigmaup>[cont ? 1 : 0] sigma;
		vector[nbase ? nbase : 0] bslope;
}

transformed parameters{
    real <lower=0> lambda;
		real loglambda;
    real<lower=0> ed50;
		real led50;
		real emax;
		// the bivariate t-prior distribution for log(ed50/p50) and log(lambda)
		// is represented as a univariate t for log(ed50/p50) and a 
		// conditional univariate t for log(lambda) given log(ed50/p50)
		// using the results in ding (tas, 2016).  this allows univariate
		// truncation of the distributions of log(ed50/p50) and log(lambda)
		real bpreg;
		real spreg;
		vector<lower=0>[cont ? 1 : 0] tau2;
		
		bpreg=parmCor*loglamsca/loged50sca;
		spreg=loglamsca*loglamsca*(1-pow(parmCor,2));

		ed50=p50*exp(loged50);
		led50=log(ed50);                 // includes p50 for output
   	if(sigmoid){
			lambda=exp(llam);
		}else{
			lambda=1;
		}
		loglambda=log(lambda);

		if(cont){
			tau2[1]=1/(2*sigma[1]*sigma[1]);
		}
		emax=(difTarget)*(ed50^lambda+dTarget^lambda)/dTarget^lambda;
}

model{
    vector[N] emx;
		vector[N] sex;
		real csd;
		real cmu;

		cmu=loglammu+bpreg*(loged50-loged50mu);
		csd=sqrt(spreg*(parmDF+pow(loged50-loged50mu,2)/pow(loged50sca,2))/(parmDF+1));
		
		if(intercept){
			for(i in 1:nprot){
				e0[i]~student_t(e0DF,epmu,epsca);
			}
		}
		
		loged50~student_t(parmDF,loged50mu,loged50sca);

    if(sigmoid){
			llam~student_t(parmDF+1,cmu,csd);	
		}

    if(cont){
			sigma[1]~uniform(sigmalow,sigmaup);
    }

    difTarget~student_t(diftDF,difTargetmu,difTargetsca);
    
    for(i in 1:N){
			if(intercept){
					emx[i] = e0[protv[i]] + (emax * pow(dv[i],lambda))/(pow(ed50,lambda)
										+pow(dv[i],lambda));
			}else{
					emx[i] = (emax * pow(dv[i],lambda))/(pow(ed50,lambda)
										+pow(dv[i],lambda));
			}
			if(cont) sex[i]=sigma[1]/sqrt(nv[i]);
    }

		if(nbase){
			bslope~multi_normal_cholesky(basemu,cbvar);
			emx = emx + xbase*bslope; //emx;
		}
		
		if(dimFit){
			yv ~ multi_normal_cholesky(emx,cholvc);
    }else if(cont){
	   	yv ~ normal(emx,sex);
			if(gp){
				ssy ~ gamma(df2,tau2[1]);
			}		
		}else{
			yvb ~ binomial(nvb,inv_logit(emx));
		}
}


