data {
    int<lower=0> N; // number of patients/groups 
    int<lower=0> nprot; // number of protocols-strata
		int<lower=0> protv[N];
		int<lower=0,upper=1> cont;  //indicator of data type
		int<lower=0,upper=1> sigmoid;  //indicator of model type
		int<lower=0,upper=1> gp;
		int<lower=0,upper=1> intercept;
		int<lower=0> nbase;
		
		//continuous data
    real yv[N]; //  outcome variable
		real<lower=0> nv[N]; // number patients per group; real for arithmetic
		//binary data
    int yvb[N]; //  outcome variable
		int<lower=0> nvb[N]; // number patients per group
		
    real<lower=0> dv[N];  
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
    real<lower=0> loged50mu;
    real<lower=0> loged50sca;
    real<lower=0> e0DF;
    real<lower=0> diftDF;
    real<lower=0> parmDF;
    real<lower=0> loglammu;
    real<lower=0> loglamsca;
		real parmCor;
		vector[nbase ? nbase : 0] basemu;
		cov_matrix[nbase ? nbase : 1] basevar;	
}
transformed data{
		cholesky_factor_cov[nbase ? nbase : 1] cbvar;
		vector[2] zvec;
		vector[2] muvec;
		cov_matrix[2] scavar;
		cholesky_factor_cov[2] cscavar;

		cbvar=cholesky_decompose(basevar);

		if(sigmoid){
				scavar[1,1]=square(loged50sca);
				scavar[1,2]=loged50sca*loglamsca*parmCor;
				scavar[2,1]=scavar[1,2];
				scavar[2,2]=square(loglamsca);
				cscavar=cholesky_decompose(scavar);
				zvec[1]=0.0;
				zvec[2]=0.0;
				muvec[1]=loged50mu;
				muvec[2]=loglammu;
		// assign defaults to avoid compiler error
		}else{  				
				scavar[1,1]=1.0;
				scavar[1,2]=0.0;
				scavar[2,1]=scavar[1,2];
				scavar[2,2]=1.0;
				cscavar=cholesky_decompose(scavar);
		}
}
parameters{
    vector[intercept ? nprot : 0] e0;
    real difTarget;
		vector[sigmoid ? 2 : 1]parmvec;
		vector<lower=sigmalow,upper=sigmaup>[cont ? 1 : 0] sigma;
		vector[nbase ? nbase : 0] bslope;
		vector<lower=0>[sigmoid ? 1 : 0] chi2var;
}
transformed parameters{
    real <lower=0> lambda;
		real loglambda;
    real<lower=0> ed50;
		real led50;
		real emax;
		real<lower=0> tau2[cont ? 1 : 0];
		vector[sigmoid ? 2 : 1]parmvect;
		
    if(sigmoid){
			parmvect=muvec+parmvec/sqrt(chi2var[1]/parmDF);
			ed50=p50*exp(parmvect[1]);
			loglambda=parmvect[2];
    	lambda=exp(parmvect[2]);
    }else{
			parmvect=loged50mu+parmvec;
			ed50=p50*exp(parmvect[1]);
			loglambda=0.0;
			lambda=1.0;
		}
		led50=log(ed50);
		if(cont){
			tau2[1]=1/(2*sigma[1]*sigma[1]);
		}
		emax=(difTarget)*(ed50^lambda+dTarget^lambda)/dTarget^lambda;
}

model{
    vector[N] emx;
		vector[N] sex;
		
		if(intercept){
			for(i in 1:nprot){
				e0[i]~student_t(e0DF,epmu,epsca);
			}
		}
    if(sigmoid){
			parmvec~multi_normal_cholesky(zvec,cscavar);
			chi2var[1]~chi_square(parmDF);
		}else parmvec[1]~student_t(parmDF,0.0,loged50sca);
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
		
    if(cont){
	   	yv ~ normal(emx,sex);
			if(gp){
				ssy ~ gamma(df2,tau2[1]);
			}		
		}else{
			yvb ~ binomial(nvb,inv_logit(emx));
		}
}


