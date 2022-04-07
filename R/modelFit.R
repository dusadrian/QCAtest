`modelFit` <- function(
    model, theory = ""
) {
    if (!(methods::is(model, "QCA_min") | methods::is(model, "admisc_deMorgan"))) {
        admisc::stopError(
            "The model should be a minimization object or its negation."
        )
    }
    theory <- admisc::recreate(substitute(theory))
    if (is.character(theory)) {
        if (length(theory) != 1) {
            admisc::stopError(
                "Theory should be a single character expression."
            )
        }
    }
    else {
        admisc::stopError(
            "Theory should be a character expression or its negation."
        )
    }
    noflevels <- model$tt$noflevels
    snames <- model$tt$options$conditions
    if (model$tt$options$use.letters) {
        snames <- LETTERS[seq(length(snames))]
    }
    pims <- model$pims
    if (is.element("i.sol", names(model))) {
        pims <- lapply(model$i.sol, function(x) x$pims)
        names(pims) <- NULL
        pims <- do.call("cbind", pims)
        solutions <- lapply(model$i.sol, function(x) x$solution)
    }
    else {
        solutions <- list(model$solution)
    }
    models <- unlist(lapply(solutions, function(x) unlist(lapply(x, paste, collapse = " + "))))
    slengths <- unlist(lapply(solutions, length))
    if (is.null(names(solutions))) {
        names(models) <- "M"
        if (slengths > 1) {
            names(models) <- paste("M", seq(slengths), sep = "")
        }
    }
    else {
        mnum <- unlist(lapply(slengths, function(x) {
            mnum <- ""
            if (x > 1) {
                mnum <- seq(x)
            }
            paste("M", mnum, sep = "")
        }))
        names(models) <- paste(mnum, rep(names(solutions), slengths), sep = "-")
    }
    result <- intersections <- vector(mode = "list", length = length(models))
    arglist <- list(snames = snames, noflevels = noflevels)
    for (i in seq(length(models))) {
        expression <- models[i]
        cpims <- pims[, unlist(strsplit(expression, split = " \\+ ")), drop = FALSE]
        cpims$model <- admisc::compute(expression, data = model$tt$initial.data)
        cpims$theory <- admisc::compute(theory, data = model$tt$initial.data)
        intersections <- rep("", 4)
        intersections[1] <- do.call(
            admisc::intersection,
            c(
                list(theory, expression),
                arglist
            )
        )
        intersections[2] <- do.call(
            admisc::intersection,
            c(
                list(
                    negate(theory, snames = snames)[[1]][1],
                    expression
                ),
                arglist
            )
        )
        intersections[3] <- do.call(
            admisc::intersection,
            c(
                list(
                    theory,
                    negate(expression, snames = snames)[[1]][1]
                ),
                arglist
            )
        )
        intersections[4] <- do.call(
            admisc::intersection,
            c(
                list(
                    negate(theory, snames = snames)[[1]][1],
                    negate(expression, snames = snames)[[1]][1]
                ),
                arglist
            )
        )
        intnms <- c("model*theory", "model*~theory", "~model*theory", "~model*~theory")
        for (nm in seq(4)) {
            int <- intersections[nm]
            if (int == "") {
                cpims[[intnms[nm]]] <- rep(0, nrow(model$tt$initial.data))
            }
            else {
                cpims[[intnms[nm]]] <- admisc::compute(int, data = model$tt$initial.data)
            }
        }
        intersections[intersections == ""] <- "-"
        names(intersections) <- intnms
        neg.out <- admisc::hastilde(model$tt$options$outcome)
        pofobj <- pof(
            cpims,
            model$tt$initial.data[, admisc::notilde(model$tt$options$outcome)],
            relation = "sufficiency",
            neg.out = neg.out
        )
        pofobj$incl.cov <- pofobj$incl.cov[, 1:3]
        pofobj$incl.cov[is.na(pofobj$incl.cov[, 1]), 3] <- NA
        pofobj$modelfit <- list(
            model = expression,
            theory = theory,
            intersections = intersections
        )
        result[[i]] <- pofobj
    }
    if (length(result) == 1) {
        return(result[[1]])
    }
    return(structure(result, class = "QCA_modelFit"))
}
