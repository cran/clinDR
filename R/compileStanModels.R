compileStanModels <- function()
{
	
	fileloc=system.file(package="clinDR","models")
	
	if(file.access(fileloc,mode=0)<0)stop(paste('Could not locate package directory ',
																							fileloc))
	
	if( file.access(fileloc,mode=2)<0 )
		stop(paste('You need write privileges for directory ', fileloc))
	
	if(isTRUE(grep("64",Sys.getenv("R_ARCH"))>0)){
		filelocnew<-file.path(fileloc,'comp64')
	}else filelocnew<-file.path(fileloc,'comp32')
	
	dir.create(filelocnew,showWarnings=FALSE)

	cat('Model compilation requires several minutes\n')
	
	capture.output(
	estancont4<-stan_model(file=file.path(fileloc,'modelcont4.stan'),
												 save_dso=TRUE,auto_write=FALSE,model_name='cont4')
	,file=NULL)
	
	saveRDS(estancont4,file=file.path(filelocnew,"estancont4.rds"))
	
	cat('estancont4.rds created\n')
	flush.console()
	
	capture.output(	
	estancont3<-stan_model(file=file.path(fileloc,'modelcont3.stan'),
												 save_dso=TRUE,auto_write=FALSE,model_name='cont3')
	,file=NULL)
	
	saveRDS(estancont3,file=file.path(filelocnew,"estancont3.rds"))
	cat('estancont3.rds created\n')
	flush.console()
	
	capture.output(	
	estancont4NoPbo<-stan_model(file=file.path(fileloc,'modelcont4NoPbo.stan'),
														save_dso=TRUE,auto_write=FALSE,model_name='cont4NoPbo')
	,file=NULL)
	
	saveRDS(estancont4NoPbo,file=file.path(filelocnew,"estancont4NoPbo.rds"))
	cat('estancont4NoPbo.rds created\n')
	flush.console()
	
	capture.output(	
	estancont3NoPbo<-stan_model(file=file.path(fileloc,'modelcont3NoPbo.stan'),
															save_dso=TRUE,auto_write=FALSE,model_name='cont3NoPbo')
	,file=NULL)
	
	saveRDS(estancont3NoPbo,file=file.path(filelocnew,"estancont3NoPbo.rds"))
	cat('estancont3NoPbo.rds created\n')
	flush.console()
	
	capture.output(	
	estanbin4<-stan_model(file=file.path(fileloc,'modelbin4.stan'),
												save_dso=TRUE,auto_write=FALSE,model_name='bin4')
	,file=NULL)
	
	saveRDS(estanbin4,file=file.path(filelocnew,"estanbin4.rds"))
	cat('estanbin4.rds created\n')
	flush.console()	
	
	capture.output(	
	estanbin3<-stan_model(file=file.path(fileloc,'modelbin3.stan'),
												save_dso=TRUE,auto_write=FALSE,model_name='bin3')
	,file=NULL)
	
	saveRDS(estanbin3,file=file.path(filelocnew,"estanbin3.rds"))
	cat('estanbin3.rds created\n')
	flush.console()	

	return(invisible())
	
}
