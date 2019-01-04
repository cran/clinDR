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
		real<lower=0> epsd;
    real emaxmu;
		real<lower=0> emaxsd;
		real<lower=0> sigmalow;
		real<lower=0> sigmaup;
		real<lower=0> p50;
    real<lower=0> led50mu;
    real<lower=0> led50sca;
    real<lower=0> edDF;
    real<lower=0> lama;
    real<lower=0> lamb;
    real<lower=0> lamsca;
		vector[nbase ? nbase : 0] basemu;
		cov_matrix[nbase ? nbase : 1] basevar;	
}
transformed data{
		cholesky_factor_cov[nbase ? nbase : 1] cbvar;
		cbvar=cholesky_decompose(basevar);
}
parameters{
    vector[intercept ? nprot : 0] e0;
    real emax;
    vector<lower=0,upper=1>[sigmoid ? 1 : 0] lamt;
    real ed50t;
		vector<lower=sigmalow,upper=sigmaup>[cont ? 1 : 0] sigma;
		vector[nbase ? nbase : 0] bslope;
}
transformed parameters{
    real<lower=0,upper=lamsca> lambda;
    real<lower=0> ed50;
		real led50;
		real<lower=0> tau2[cont ? 1 : 0];
		
    if(sigmoid){
    	lambda=lamsca*lamt[1];
    }else lambda=1;
    ed50=p50*exp(ed50t);
		led50=log(ed50);
		if(cont){
			tau2[1]=1/(2*sigma[1]*sigma[1]);
		}
}

model{
    vector[N] emx;
		vector[N] sex;
		
		if(intercept){
			for(i in 1:nprot){
				e0[i]~normal(epmu,epsd);
			}
		}
    if(sigmoid)lamt[1]~beta(lama,lamb);
    if(cont){
			sigma[1]~uniform(sigmalow,sigmaup);
    }

    ed50t~student_t(edDF,led50mu,led50sca);
    emax~normal(emaxmu,emaxsd);
    
    for(i in 1:N){
			emx[i] = (emax * pow(dv[i],lambda))/(pow(ed50,lambda)
								+pow(dv[i],lambda));
			if(intercept)emx[i] = emx[i]+e0[protv[i]];
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


