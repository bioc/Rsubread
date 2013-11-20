RsubreadUsersGuide <- function()
{
	f <- system.file("usersguide","SubreadUsersGuide.pdf",package="Rsubread")
	system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	return(f)
}
