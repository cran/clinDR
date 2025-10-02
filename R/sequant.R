"sequant" <-
function(p, sig2, n)
{
#
#  compute se for pth quantile
#
#  sig2 is distribution variance estimate
#  (preferably by .75*IQR)
#
#  n number in sample
#
	pq <- p * (1 - p)
	phip <- (qnorm(p))^2
	seqret <- sig2 * pq * 2 * pi * exp(phip)
	seqret <- sqrt(seqret/n)
	return(seqret)
}

