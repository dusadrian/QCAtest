`replaceCategories` <- function(x, categories = NULL) {
    mv <- any(grepl("\\[|\\{", x)) && all(grepl("\\]|\\}", x))
    target <- replacement <- c()
    nms <- names(categories)
    for (i in seq(length(categories))) {
        if (mv) {
            values <- seq(length(categories[[i]])) - 1
            target <- c(target, paste0(nms[i], "[", values, "]"))
        }
        else {
            target <- c(target, paste0("~", nms[i]), nms[i])
        }
        replacement <- c(replacement, categories[[i]])
    }
    for (i in seq(length(x))) {
        x[i] <- replaceText(x[i], target, replacement)
    }
    return(x)
}
