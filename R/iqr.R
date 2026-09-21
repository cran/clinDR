"iqr" <-
function(y){
	return(quantile(y,.75,na.rm= TRUE)-quantile(y,.25,na.rm= TRUE))
}

