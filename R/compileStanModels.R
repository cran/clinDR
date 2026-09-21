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

	cat('Model compilation may take several minutes\n')
	
	capture.output(
	basemodel<-stan_model(file=file.path(fileloc,'basemodel.stan'),
												 save_dso=TRUE,auto_write=FALSE,model_name='base')
	,file=NULL)
	
		capture.output(
	mrmodel<-stan_model(file=file.path(fileloc,'mrmodel.stan'),
												 save_dso=TRUE,auto_write=FALSE,model_name='mrmod')
	,file=NULL)
		
	saveRDS(basemodel,file=file.path(filelocnew,"basemodel.rds"))
	saveRDS(mrmodel,file=file.path(filelocnew,"mrmodel.rds"))
	
	cat('basemodel.rds and mrmodel.rds created\n')
	flush.console()

	return(invisible())
	
}
