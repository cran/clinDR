#' @export
runSimulations<- function() {
	appDir <- system.file("shiny", "", package = "clinDR")
	if (appDir == "") {
		stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
	}
	
	shiny::runApp(appDir, display.mode = "normal")
}