verbose_message_factory <- function(verbose) {
    if (verbose) {
        function(...) {
            message(...)
        }
    } else {
        function(...) {
            invisible(NULL)
        }
    }
}

summ_beta <- function(fit_l, weights = rep(1 / length(fit_l), length(fit_l)), drop_0 = FALSE) {
    if (drop_0) {
        list(tibble::tibble(w = weights, Beta = fit_l) %>%
             tidyr::unnest(Beta) %>%
             dplyr::filter(Beta != 0) %>%
             dplyr::group_by(feature_name) %>%
             dplyr::summarise(
                           Beta_sd = sd(Beta),
                        Beta = stats::weighted.mean(Beta, w = w / sum(w))
                       ))
    }else{
        list(tibble::tibble(w = weights, Beta = fit_l) %>%
                tidyr::unnest(Beta) %>%
                dplyr::group_by(feature_name) %>%
                dplyr::summarise(
                        Beta_sd = sd(Beta),
                        Beta = stats::weighted.mean(Beta, w = w / sum(w))
                ))

    }
}


summarise_cv_lik <- function(fit,summarise_Beta=TRUE) {
    mgrps <- dplyr::groups(fit)
    if (!summarise_Beta) {
        fit %>%
            dplyr::group_by(group_l1, group_l2, .add = TRUE) %>%
            dplyr::summarize(
                       cv_sum = sum(cv_lik),
                       l0_mean = mean(l0n),
                       l0_sd = sd(l0n),
                       ) %>%
            dplyr::group_by(!!!mgrps)
    }else{
        fit %>%
            dplyr::group_by(group_l1, group_l2, .add = TRUE) %>%
            dplyr::summarize(
                       Beta = summ_beta(Beta,
                                        weights = cv_lik / sum(cv_lik),
                                        drop_0 = FALSE),
                       cv_sum = sum(cv_lik),
                       l0_mean = mean(l0n),
                       l0_sd = sd(l0n)
                   ) %>%
            dplyr::group_by(!!!mgrps)
    }
}

subset_long <- function(fit, long_df) {
        dplyr::semi_join(long_df,
                tidyr::unnest(fit, Beta) %>%
                        dplyr::filter(Beta != 0),
                by = "feature_name"
        )
}

create_X <- function(fit, long_df, BF_df) {
        rx <- dplyr::semi_join(long_df, tidyr::unnest(fit, Beta), by = "feature_name") %>%
                join_long_wide(BF_df)
}



pull_beta <- function(Beta){
    magrittr::set_names(
            c(
                    Beta$Beta[Beta$feature_name == "Intercept"],
                    Beta$Beta[Beta$feature_name != "Intercept"]
            ),
            c(
                    "Intercept",
                    Beta$feature_name[Beta$feature_name != "Intercept"]
            )
            )
}


predict_null_fgem <- function(BF, log = TRUE) {
        if (log) {
                eBF <- exp(BF)
        } else {
                eBF <- BF
        }
        rv <- gen_u(
                Beta = fgem:::fgem_null(eBF)[1],
                x = matrix(numeric(), nrow = length(BF), ncol = 0),
                B = BF,
                log = log
        )
}




join_long_wide <- function(long_df, wide_df, key="feature_name", value="value", by = "Gene") {

    ldv <- long_df[[value]]
    if (is.null(ldv))
        ldv <- rep(1.0, nrow(long_df))

    if( length(unique(wide_df[[by]])) < NROW(wide_df) ){
        u_f <- unique(long_df[[key]])
        ixm <- dplyr::mutate(long_df, tv = 1) %>%
            tidyr::spread(key = {{ key }}, value = "tv", fill = 0L)
        ix_df <- dplyr::mutate(
                            dplyr::select(ixm, {{ by }}),
                            X = as.matrix(dplyr::select(ixm, -{{ by }}))) %>%
            dplyr::inner_join(wide_df) %>%
                    dplyr::relocate(X, .after = dplyr::last_col())
        return(ix_df)
    }

    ixm <- trip2sparseMatrix(
        rowname_vals = long_df[[by]],
        colname_vals = long_df[[key]],
        values = ldv,
        total_rownames = unique(wide_df[[by]]),
        total_colnames = unique(long_df[[key]]),
        add_intercept = FALSE
    )
    dplyr::mutate(wide_df, X = ixm)
}

predict_long_df <- function(beta_df, x, log = FALSE) {
    beta_df <- beta_df %>% dplyr::filter( Beta != 0) %>%
        dplyr::arrange(feature_name != "Intercept")

    beta <- beta_df$Beta
    stopifnot(
        nrow(beta_df) > 0,
        beta_df$feature_name[1] == "Intercept"
    )
    if (NCOL(x) > 0) {
        x <- x[, beta_df$feature_name[beta_df$feature_name != "Intercept"], drop = FALSE]
    }
    gen_p(beta, x, log = log) # rownames(x))

}



##' @title Predict posterior from fgem fit
##'
##'
##' @param fit result of call to `fgem`
##' @param x annotation matrix (colnames should correspond to `feature_name` column of fit$Beta[[1]])
##' @param BF Bayes Factor
##' @param lambda value of lambda to use in prediction
##' @param log return log posterior
##' @return vector of length length(BF) with posterior probabilities
##' @export
predict_fgem <- function(fit, x, log = FALSE) {
        beta_l <- fit[["Beta"]]
        stopifnot(
                !is.null(beta_l),
                is.list(beta_l)
        )
        beta_df <- fit$Beta[[1]]
        predict_long_df(beta_df, x, log = log)
}



empty_x <- function(x, rownames = names(x)) {
        nr <- length(x)
        if (!is.null(rownames)) {
                nr <- min(nr, length(rownames))
        }
        matrix(0, nrow = nr, ncol = 0, dimnames = list(rownames, NULL))
}

predict_uniform_prior <- function(BF, Genes = names(BF)) {
    predict_fgem(fgem_null_fit(BF), empty_x(BF, Genes))
}

BF2posterior <- function(BF,Genes=names(BF), prior = predict_uniform_prior(BF,Genes)) {
    return((prior * BF) / ((prior * BF) + (1 - prior)))
}


log_BF2_logpost <- function(logBF, xb) {
        apply(xb - logBF, 2, log_1p_exp)
}




##' @title Convert a table-style sparse matrix representation to a sparse matrix
##'
##'
##'
##' @param rowname_vals vector with the rownames of the values
##' @param colname_vals vector with the colnames of the values
##' @param values the values themselves (default is 1.0)
##' @param total_rownames the universe of rownames, defaults to all the unique rowname_vals
##' @param total_colnames the universe of colnames, defaults to all the unique colname_vals
##' @param add_intercept whether to add an intercept term ("Intercept" will be added to total_colnames)
##' @return a sparse matrix of dim(length(total_rownames),length(total_colnames)) with (at most) length(values)
##' non-zero values
##' @export
trip2sparseMatrix <- function(rowname_vals,
                              colname_vals,
                              values = rep(1.0,length(colname_vals)),
                              total_rownames = unique(rowname_vals),
                              total_colnames = unique(colname_vals),
                              add_intercept = TRUE,
                              csr=FALSE) {

    stopifnot(
            all.equal(length(rowname_vals), length(colname_vals)),
            all.equal(length(rowname_vals), length(values))
        )

    bad_vals <- (!(rowname_vals %in% total_rownames)) | (!(colname_vals %in% total_colnames))

    rowname_vals <- rowname_vals[!bad_vals]
    colname_vals <- colname_vals[!bad_vals]
    values <- values[!bad_vals]

    if (add_intercept) {
        total_colnames <- unique(c(total_colnames, "Intercept"))
        rowname_vals <- c(rowname_vals, total_rownames)
        colname_vals <- c(colname_vals, rlang::rep_along(total_rownames, "Intercept"))
        values <- c(values, rlang::rep_along(total_rownames, 1.0))
    }

    row_id <- purrr::set_names(seq_along(total_rownames), total_rownames)
    col_id <- purrr::set_names(seq_along(total_colnames), total_colnames)
    ret_m <- Matrix::sparseMatrix(
            i = row_id[as.character(rowname_vals)],
            j = col_id[as.character(colname_vals)],
            dims = c(length(total_rownames), length(total_colnames)),
            x = values, dimnames = list(total_rownames, total_colnames)
            )
    ret_m
}

df2sparse <- function(long_df,rownames="Gene",colnames="feature_name",value="value",total_rownames=unique(long_df[[rownames]])){

    rnv <- long_df[[rownames]]
    cnv <- long_df[[colnames]]
    stopifnot(
            !is.null(rnv),
            !is.null(cnv)
    )

    ldv <- long_df[[value]]
    if (is.null(ldv))
        ldv <- rep(1.0, nrow(long_df))

    trip2sparseMatrix(rnv, cnv, ldv, total_rownames = total_rownames, add_intercept = FALSE)
}
