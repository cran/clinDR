data {
    int<lower=0> N; // number of patients/groups 
    int<lower=0> nprot; // number of protocols-strata
		int<lower=0> prot[N];
    real yv[N]; //  outcome variable
		real<lower=0> count[N]; // number patients per group; real for arithmetic
    real<lower=0> dv[N];  
		int<lower=0,upper=1> gp;
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
}
parameters{
    real e0[nprot];
    real emax;
    real ed50t;
		real<lower=sigmalow,upper=sigmaup> sigma;
}
transformed parameters{
    real<lower=0> ed50;
		real led50;
		real<lower=0> tau2;

    ed50=p50*exp(ed50t);
		led50=log(ed50);
		tau2=1/(2*sigma*sigma);
}

model{
    real emx[N];
		real sex[N];

		for(i in 1:nprot){
			e0[i]~normal(epmu,epsd);
		}

    emax~normal(emaxmu,emaxsd);
    ed50t~student_t(edDF,led50mu,led50sca);
		sigma~uniform(sigmalow,sigmaup);

    for (i in 1:N) {
        emx[i] = e0[prot[i]] 
					+ (emax*dv[i])/(ed50+dv[i]);
				sex[i]=sigma/sqrt(count[i]);
    }
    yv ~ normal(emx,sex);
		if(gp){
			ssy ~ gamma(df2,tau2);
		}
}



