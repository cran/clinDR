//
// fit model to impute binary longitudinal missing data
//

data{
	int<lower=1>		N;				// Tot Num of measurement 
	int<lower=1>    nsubj;		// Num of subject
	int<lower=1>    nvisit;	// Num of visits 
	int<lower=1>		ntrt;			//Num of trt groups
	int<lower=1>		ntrt1;			// num treatments minus 1

	array[N] int<lower=1,upper=nsubj>	id; // patient index (sequential 1,2,3...,nsubj)
	array[N] int<lower=0,upper=ntrt>	trt;		// trt(0=cntl, 1=trt)
	array[N] int<lower=1,upper=nvisit> visit;	// visit num (sequential, 1,2,3...,nvisit)

	array[N] int<lower=0,upper=1>				y;		//1=resp,  0=no resp

														// priors on logit scale
	array[nvisit] real						prmean0;		// prior mean for cntl at each visit 
	array[nvisit] real						prsd0;			// prior sd for cntl at each visit 
	array[nvisit] real			prmean;		// prior mean for trt effect at each visit 
	array[nvisit] real			prsd;			// prior sd for trt effect at each visit 
	real						gparm1;           // sigma gamma prior shape parameter
	real						gparm2;           // sigma gamma prior rate parameter

}

parameters{
	
	matrix[nvisit,ntrt1]  trtdel;
	vector[nvisit] cntl;

	real						sigmasq;

	array[nsubj] real						theta;  // random patient level
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

