//
// fit model to impute binary longitudinal missing data
//

data{
	int<lower=1>		N;				// Tot Num of measurement 
	int<lower=1>    nsubj;		// Num of subject
	int<lower=1>    nvisit;	// Num of visits 
	int<lower=1>		ntrt;			//Num of trt groups
	int<lower=1>		ntrt1;			// num treatments minus 1

	int<lower=1,upper=nsubj>		id[N];			// patient index (sequential 1,2,3...,nsubj)
	int<lower=0,upper=ntrt>			trt[N];		// trt(0=cntl, 1=trt)
	int<lower=1,upper=nvisit>		visit[N];	// visit num (sequential, 1,2,3...,nvisit)

	int<lower=0,upper=1>				y[N];		//1=resp,  0=no resp

														// priors on logit scale
	real						prmean0[nvisit];		// prior mean for cntl at each visit 
	real						prsd0[nvisit];			// prior sd for cntl at each visit 
	real						prmean[nvisit];		// prior mean for trt effect at each visit 
	real						prsd[nvisit];			// prior sd for trt effect at each visit 
	real						gparm1;           // sigma gamma prior shape parameter
	real						gparm2;           // sigma gamma prior rate parameter

}

parameters{
	
	matrix[nvisit,ntrt1]  trtdel;
	vector[nvisit] cntl;

	real						sigmasq;

	real						theta[nsubj];  // random patient level
}

transformed parameters{
	matrix[nvisit,ntrt]  beta;
	real sigma;

	sigma = sqrt(sigmasq);

	beta[,1] = cntl;
	for(i in 2:ntrt){
		beta[,i] = cntl + trtdel[,i-1];
	}
}

model{
	
	cntl ~ normal(prmean0,prsd0);
	for (i in 1:ntrt1){
			trtdel[,i]~normal(prmean,prsd);
	}


	sigmasq~inv_gamma(gparm1,gparm2);
	theta~normal(0,1.0);    // centered parameterization much better

	{
	vector[N] logitx;
	for (i in 1:N){
		logitx[i] = beta[visit[i],trt[i]] + sigma*theta[id[i]];
	}
	y~bernoulli_logit( logitx );
	}
}

