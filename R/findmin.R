`findmin` <- function(
    chart, ...
) {
    dots <- list(...)
    if (!methods::is(chart, "QCA_pic")) {
        if (!is.matrix(chart) | (!is.logical(chart) & length(setdiff(chart, 0:1)) > 0)) {
            admisc::stopError(
                "Use a logical, TRUE/FALSE matrix. See makeChart()'s output."
            )
        }
    }
    if (all(colSums(chart) > 0)) {
        result <- .Call("C_findmin", matrix(as.logical(chart), nrow = nrow(chart)), PACKAGE = "QCAtest")
    }
    else {
        result <- 0
    }
    class(result) <- c("numeric", "QCA_findmin")
    return(result)
}
