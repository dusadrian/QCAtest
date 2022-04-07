`fuzzyor` <- function(
    ..., na.rm = FALSE
) {
    funargs <- unlist(
        lapply(
            lapply(match.call(), deparse)[-1],
            function(x) gsub("\"|[[:space:]]", "", x)
        )
    )
    if (!is.na(rem <- match("na.rm", names(funargs)))) {
        funargs <- funargs[-rem]
    }
    dots <- vector(mode = "list", length = length(funargs))
    funargs <- gsub(rawToChar(as.raw(c(226, 128, 147))), "-", funargs)
    negated <- grepl("1-", funargs)
    funargs <- gsub("1-", "", funargs)
    tildenegated <- badnames <- cols <- logical(length(funargs))
    for (i in seq(length(funargs))) {
        badnames[i] <- grepl("\\(|:", funargs[i])
        cols[i] <- admisc::getName(admisc::notilde(funargs[i]))
        tildenegated[i] <- admisc::tilde1st(funargs[i])
        funargs[i] <- admisc::notilde(funargs[i])
    }
    if (sum(badnames) > 0) {
        if (sum(badnames) > length(LETTERS) | any(is.element(cols, LETTERS))) {
            cols[badnames] <- paste("X", seq(sum(badnames)), sep = "")
        }
        else {
            cols[badnames] <- LETTERS[seq(sum(badnames))]
        }
    }
    for (i in seq(length(funargs))) {
        tc <- tryCatch(
            eval.parent(
                parse(text = funargs[i])
            ),
            error = function(e) e,
            warning = function(w) w
        )
        if (is.function(tc) | inherits(tc, "error")) {
            admisc::stopError(sprintf(
                "Object '%s' not found.", funargs[i])
            )
        }
        else {
            dots[[i]] <- eval.parent(parse(text = funargs[i]), n = 1)
        }
    }
    if (is.element("name", names(attributes(dots[[1]])))) {
        dots[[1]] <- as.vector(dots[[1]])
    }
    if (is.vector(dots[[1]])) {
        if (
            any(
                !unlist(
                    lapply(
                        dots,
                        function(x) is.numeric(x) | is.logical(x)
                    )
                )
            )
        ) {
            admisc::stopError(
                "Input vectors should be numeric or logical."
            )
        }
        dots <- as.data.frame(dots)
    }
    else if (is.matrix(dots[[1]])) {
        dots <- dots[[1]]
        if (is.null(colnames(dots))) {
            if (ncol(dots) > length(LETTERS)) {
                cols <- paste("X", seq(ncol(dots)), sep = "")
            }
            else {
                cols <- LETTERS[seq(ncol(dots))]
            }
        }
        dots <- as.data.frame(dots)
        negated <- logical(ncol(dots))
        tildenegated <- logical(ncol(dots))
        if (
            !all(
                unlist(
                    lapply(
                        dots,
                        function(x) is.numeric(x) | is.logical(x)
                    )
                )
            )
        ) {
            admisc::stopError(
                "Input should be numeric or logical."
            )
        }
    }
    else if (is.data.frame(dots[[1]])) {
        dots <- dots[[1]]
        negated <- logical(ncol(dots))
        tildenegated <- logical(ncol(dots))
        cols <- colnames(dots)
        if (
            !all(
                unlist(
                    lapply(
                        dots,
                        function(x) is.numeric(x) | is.logical(x)
                    )
                )
            )
        ) {
            admisc::stopError(
                "Some columns are not numeric or logical."
            )
        }
    }
    else {
        admisc::stopError(
            "The input should be vectors, or a matrix or a dataframe."
        )
    }
    for (i in seq(length(cols))) {
        if (tildenegated[i]) {
            dots[[i]] <- 1 - dots[[i]]
        }
        if (negated[i]) {
            dots[[i]] <- 1 - dots[[i]]
        }
        if (negated[i] + tildenegated[i] == 1) {
            cols[i] <- paste("~", cols[i], sep = "")
        }
    }
    result <- apply(dots, 1, max, na.rm = na.rm)
    attr(result, "names") <- NULL
    attr(result, "name") <- paste(cols, collapse = " + ")
    class(result) <- c("numeric", "QCA_fuzzy")
    return(result)
}
