`getSolution` <- function(
    expressions, mv, collapse, inputt, row.dom, initial, all.sol, indata, curly, categorical, ...
) {
    mtrx <- NULL
    sol.matrix <- NULL
    dots <- list(...)
    enter <- ifelse (is.element("enter", names(dots)), dots$enter, "\n")
    pi.cons <- if (is.element("pi.cons", names(dots))) dots$pi.cons else 0
    outcome <- if (is.element("outcome", names(dots))) dots$outcome else ""
    complex <- FALSE
    if (is.list(expressions)) {
        mtrx <- expressions[[2]]
        sol.matrix <- expressions[[3]]
        if (length(expressions) > 3) {
            complex <- expressions[[4]]
        }
        if (is.null(sol.matrix)) {
            admisc::stopError(
                "There are no solutions, given these constraints.",
                enter
            )
        }
        expressions <- expressions[[1]]
        if (nrow(unique(expressions)) != nrow(expressions)) {
            expressions <- unique(expressions)
            mtrx <- NULL
            sol.matrix <- NULL
        }
    }
    if (
        nrow(expressions) == 1 &&
        identical(
            unique(
                as.vector(expressions)
            ),
            0L
        )
    ) {
        admisc::stopError(
            paste(
                "All truth table configurations are used, all conditions are minimized.",
                "       Please check the truth table.",
                sep = "\n"
            ),
            enter
        )
    }
    if (FALSE) {
        if (!missing(indata)) {
            hastime <- logical(ncol(expressions))
            for (i in seq(ncol(expressions))) {
                if (any(is.element(indata[, i], c("-", "dc", "?")))) {
                    hastime[i] <- TRUE
                }
            }
            indata <- indata[, !hastime, drop = FALSE]
            expressions <- expressions[, !hastime, drop = FALSE]
            inputt <- inputt[, !hastime, drop = FALSE]
            relevant <- apply(expressions, 1, sum) > 0
            if (any(!relevant)) {
                sol.matrix <- NULL
                mtrx <- mtrx[relevant, , drop = FALSE]
                expressions <- expressions[relevant, , drop = FALSE]
            }
        }
    }
    PI <- admisc::writePrimeimp(
        impmat = expressions,
        mv = mv,
        collapse = collapse,
        curly = curly
    )
    rownames(expressions) <- PI
    if (pi.cons > 0 & outcome != "") {
        pofPI <- pof(
            paste(
                PI,
                collapse = " + "
            ),
            outcome = indata[, dots$outcome],
            data = indata,
            relation = "sufficiency",
            categorical = categorical
        )
        inclS <- pofPI$incl.cov[seq(length(PI)), 1]
        filterPI <- admisc::agteb(inclS, pi.cons)
        expressions <- expressions[filterPI, , drop = FALSE]
        PI <- PI[filterPI]
        mtrx <- makeChart(
            expressions,
            inputt,
            mv = mv,
            collapse = collapse,
            getSolution = TRUE,
            curly = curly
        )
        if (any(colSums(mtrx) == 0)) {
            admisc::stopError(
                "There are no solutions, given these constraints.",
                enter
            )
        }
    }
    if (is.null(mtrx)) {
        mtrx <- makeChart(
            expressions,
            inputt,
            mv = mv,
            collapse = collapse,
            getSolution = TRUE,
            curly = curly
        )
    }
    else {
        rownames(mtrx) <- PI
    }
    notempty <- apply(mtrx, 1, any)
        expressions <- expressions[notempty, , drop = FALSE]
        mtrx <- mtrx[notempty, , drop = FALSE]
    setColnames(mtrx, initial)
    reduced <- list(expressions = expressions, mtrx = mtrx)
    if (nrow(mtrx) > 0) {
        if (row.dom & is.null(sol.matrix)) {
            reduced.rows <- rowDominance(mtrx)
            if (length(reduced.rows) > 0) {
                reduced$mtrx <- mtrx[reduced.rows, , drop = FALSE]
                reduced$expressions <- expressions[reduced.rows, , drop = FALSE]
            }
        }
        mtrx <- reduced$mtrx
        setColnames(mtrx, initial)
        if (is.null(sol.matrix)) {
            sol.matrix <- solveChart(mtrx, all.sol = all.sol, ... = ...)
        }
        tokeep <- sort(unique(as.vector(unique(sol.matrix))))
        all.PIs <- rownames(mtrx)[tokeep]
        solm <- matrix(as.integer(sol.matrix), nrow = nrow(sol.matrix))
        sol.matrix[sol.matrix == 0] <- NA
        sol.matrix <- matrix(rownames(mtrx)[sol.matrix], nrow = nrow(sol.matrix))
        reduced$expressions <- reduced$expressions[tokeep, , drop = FALSE]
        solution.list <- writeSolution(sol.matrix, mtrx)
    }
    else {
        all.PIs <- NA
        solution.list <- NA
        solm <- NA
    }
    return(
        list(
            expressions = expressions,
            mtrx = mtrx,
            reduced = reduced,
            all.PIs = all.PIs,
            solution.list = solution.list,
            sol.matrix = solm,
            complex = complex
        )
    )
}
