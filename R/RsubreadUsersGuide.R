RsubreadUsersGuide <- function()
{
	f <- system.file("doc","SubreadUsersGuide.pdf",package="Rsubread")
	system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	return(f)
}
