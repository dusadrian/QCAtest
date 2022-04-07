`runGUI` <- function(
    x
) {
    if (missing(x)) {
        x <- system.file("gui", PACKAGE = "QCAtest")
    }
    Sys.setenv(userwd = getwd())
    runApp(x, launch.browser = TRUE)
}
