data {
    int<lower=0> N; // number of patients/groups 
    int<lower=0> nprot; // number of protocols-strata
		int<lower=0> protv[N];
    int yv[N]; //  outcome variable
		int<lower=0> nv[N]; // number patients per group; real for arithmetic
    real<lower=0> dv[N];  
    real epmu;
		real<lower=0> epsd;
    real emaxmu;
		real<lower=0> emaxsd;
		real<lower=0> p50;
    real<lower=0> led50mu;
    real<lower=0> led50sca;
    real<lower=0> edDF;
    real<lower=0> lama;
    real<lower=0> lamb;
    real<lower=0> lamsca;
}
parameters{
    real e0[nprot];
    real emax;
    real<lower=0,upper=1> lamt;
    real ed50t;
}
transformed parameters{
    real<lower=0,upper=lamsca> lambda;
    real<lower=0> ed50;
		real led50;

    lambda=lamsca*lamt;
    ed50=p50*exp(ed50t);
		led50=log(ed50);
}

model{
    real emx[N];

		for(i in 1:nprot){
			e0[i]~normal(epmu,epsd);
		}
    emax~normal(emaxmu,emaxsd);
    lamt~beta(lama,lamb);
    ed50t~student_t(edDF,led50mu,led50sca);

    for (i in 1:N) {
        emx[i] = inv_logit( e0[protv[i]] 
        + (emax * pow(dv[i],lambda))/(pow(ed50,lambda)+pow(dv[i],lambda)) );
    }
    yv ~ binomial(nv,emx);
}

