


FGEM_Logit <- function(Beta, x, B, ctrl = list()) {
    uvec <- gen_u(Beta, x, B)
    speedglm::speedglm.wfit(
            y = uvec,
            X = x,
            family = stats::quasibinomial(link = "logit")
    )[["coefficients"]]
  ## return(stats::coefficients(stats::glm.fit(
  ##         x = x,
  ##         y = uvec,
  ##         family = stats::quasibinomial(link = "logit"),
  ##         control = ctrl
  ## )))
}

## \if{html}{\out{
## <script id="MathJax-script" async
##    src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
## </script>
## }}

obj_factory <- function(feat_mat, BF, verbose) {
        tfm <- feat_mat
        tbf <- BF
        mvmessage <- verbose_message_factory(verbose)
        old_beta <- rep(0, NCOL(feat_mat))
        old_obj <- -Inf
        ret_f <- function(par, ...) {
            beta_diff <- sum((par - old_beta)^2)
            mvmessage("Beta: ", paste0(par, collapse = ", "))
            mvmessage("Beta_diff: ", beta_diff)
            ret <- (-FGEM_Logit_log_lik(Beta = par, x = tfm, B = tbf))
            obj_diff <- sum(abs(ret - old_obj)^2)
            mvmessage("obj_diff: ", obj_diff)
            old_obj <<- ret
            old_beta <<- par
            return(-ret)
        }
        return(ret_f)
}


EM_mat <- function(Beta0, feat_mat, BF, verbose=FALSE, em_methods=c("squarem"),  ...) {
    ctrl_l <- list()
    if (verbose) {
        ctrl_l[["trace"]] <-  TRUE
    }
    vmessage <- verbose_message_factory(verbose)
    nobj_f <- obj_factory(feat_mat, BF, verbose = verbose)
    vmessage("finding Beta using: ", paste0(em_methods, collapse = ","))
    cn <- colnames(feat_mat) %||% names(Beta0)
    stopifnot(!is.null(cn))
    opf <- turboEM::turboem(
            par = Beta0,
            fixptfn = FGEM_Logit,
            boundary = boundary,
            objfn = nobj_f,
            x = feat_mat,
            method = em_methods,
            control.run = list(convtype = "objfn", tol = 1e-7),
            B = BF
            )

    if (!is.null(dim(pars(opf)))){
        best_opf <- which.max(opf$value.objfn)
        objf <- -opf$value.objfn[best_opf]
        Beta <- c(pars(opf)[best_opf, ])
        converged <- is.na(opf$errors[best_opf])
    }else{
        objf <- opf$value.objfn
        Beta <- pars(opf)
        converged <- is.na(opf$errors)
    }
    if (!converged) {
            warning("with features: ", paste0(colnames(feat_mat), collapse = ","))
            error(opf)
    }
    opdf <- dplyr::tibble(
            data = list(tibble::tibble(
                    Beta = Beta,
                    feature_name = cn
            )),
            LogLik = objf,
            convergence = converged
    )
    return(opdf)
}


FGEM_Logit_log_lik <- function(Beta, x, B) {
        pvec <- gen_p(Beta = Beta, x = x)
        return(sum(log(pvec * B + (1 - pvec))))
}

gen_p <- function(Beta, x) {
    pvec <- 1 / (1 + exp(-(x %*% Beta)))
    return(as.vector(pvec))
}

gen_u <- function(Beta, x, B) {
        p <- gen_p(Beta, x)
        return(as.vector((p * B) / ((p * B) + (1 - p))))
}


#' Generate prior from fitted FGEM_df result and annotation dataframe
#' @param feat_df annotation dataframe (as described inf `FGEM_df`)
#' @param fgem_result_df result of running FGEM_df
gen_prior <- function(feat_df, fgem_result_df) {
    Beta0 <- t(as.matrix(tidyr::unnest(fgem_result_df) %>%
                         dplyr::select(Beta, feature_name) %>%
                         spread(feature_name, Beta)))
  # BF <- feat_df[["BF"]]
  feat_df <- mutate(feat_df, Intercept = 1)
  data_mat <- data.matrix(feat_df[, rownames(Beta0)])
  prior_n <- c(gen_p(Beta = Beta0, x = data_mat))
  return(mutate(feat_df, prior = prior_n) %>% dplyr::select(-Intercept))
}


#' Generate posterior from fitted FGEM_df result and annotation dataframe
#' @param feat_df annotation dataframe (as described inf `FGEM_df`)
#' @param fgem_result_df result of running FGEM_df
gen_posterior <- function(feat_df, fgem_result_df) {
    Beta0 <- t(as.matrix(tidyr::unnest(fgem_result_df) %>%
                         dplyr::select(Beta, feature_name) %>%
                         spread(feature_name, Beta)))
  BF <- feat_df[["BF"]]
  feat_df <- mutate(feat_df, Intercept = 1)
  data_mat <- data.matrix(feat_df[, rownames(Beta0)])
  post_n <- c(gen_u(Beta = Beta0, x = data_mat, B = BF))
  return(mutate(feat_df, posterior = post_n) %>% dplyr::select(-Intercept))
}

pmean <- function(Beta, feat_mat) {
        pvec <- 1 / (1 + exp(-(feat_mat %*% Beta)))
        return(mean(pvec))
}


fgem_null <- function(BF,prior=0.02, tmu = prior_mean(BF,prior)) {
    tsp <- Matrix::Matrix(rep(1,length(BF)), nrow = length(BF), ncol = 1)
    colnames(tsp) <- "Intercept"
    tbeta0 <- speedglm::speedglm.wfit(
                            y = as.vector(tmu),
                            X = tsp,
                            family = stats::quasibinomial(link = "logit")
                        )[["coefficients"]]
    res <- fgem:::EM_mat(tbeta0, tsp, BF)
    res$data[[1]]$Beta
}

prior_mean <- function(BF, prior = 0.02) {
        (prior * BF) / ((prior * BF) + (1 - prior))
}



guess_beta0 <- function(X,BF,prior=0.02,tmu = prior_mean(BF,prior)){
    speedglm::speedglm.wfit(
                  y = tmu,
                  X = X,
                  family = stats::quasibinomial(link = "logit")
              )[["coefficients"]]
}

##' @title Run FGEM
##'
##'
##' @param X matrix of gene-level annotations
##' @param BF vector of gene-level bayes factors
##' @param prior prior probability a gene is causal (before annotations)
##' @param tmu precomputed gene-level prior
##' @param null_beta model to test against for likelihood ratio test
##' @param verbose whether to use verbose output
##' @return tibble with fgem results
##' @export
fgem <- function(X,
                 BF,
                 prior=0.02,
                 tmu = prior_mean(BF, prior),
                 null_beta=fgem_null(BF, tmu),
                 verbose = FALSE) {
    stopifnot(NROW(X) == length(BF))
    vmessage <- verbose_message_factory(verbose)
    vmessage("starting feature:", col_name)
    Beta0 <- guess_beta0(X, BF, prior, tmu)
    vmessage("Beta0", paste0(Beta0, collapse = ","))
    ret_fgem <- FGEM(
            Beta0 = Beta0,
            feat_mat = X,
            BF = BF,
            null_beta = null_beta,
            verbose = verbose
    )
    return(ret_fgem)
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
                              add_intercept = TRUE) {

    stopifnot(
            all.equal(length(rowname_vals), length(colname_vals)),
            all.equal(length(rowname_vals), length(values))
        )

    bad_vals <- (!(rowname_vals %in% total_rownames)) | (!(colname_vals %in% total_colnames))
    if (any(bad_vals)) {
            warning("rowname_vals and/or colname_vals not found in total_rownames/total_colnames, dropping: ", sum(bad_vals))
    }


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

    Matrix::sparseMatrix(
            i = row_id[rowname_vals],
            j = col_id[colname_vals],
            dims = c(length(total_rownames), length(total_colnames)),
            x = values, dimnames = list(total_rownames, total_colnames)
    )
}



##' Run FGEM with a dataframe as input
##' This version of FGEM takes a single dataframe with
##' both gene-level bayes factors and gene-level features
##' and returns a dataframe with effect size estimates
##' for every feature (run as univariate)
##' @param feat_df dataframe with data. `feat_df`
##' must have a column called `BF` specifying bayes factors.
##'   Any additional columns (besides an optional gene-name
##' column that must be titled `Gene`) are treated as annotations.
##' @param prior_mean scalar between 0 and 1 specifying the starting
##' prior probability (this is for initializing the EM algorithm).
##' If unspecified, it defaults to 0.02
##' @param verbose Whether to print debug output
##' @param ...
##' @export
FGEM_marginal <- function(X,BF,prior_mean = 0.02, verbose = FALSE,parallel=FALSE, ...) {
    if(nrow(X) != length(BF)){
        stopifnot(
            !is.null(names(BF)),
            !is.null(rownames(X)),
            all(names(BF) %in% rownames(X))
        )
        X <- X[names(BF),]
    }
    vmessage <- verbose_message_factory(verbose)
    tmu <- (prior_mean * BF) / ((prior_mean * BF) + (1 - prior_mean))
    vmessage("setting null_beta...")
    sm <- Matrix::Matrix(rep(1,length(BF)),nrow=length(BF),ncol=1,sparse = TRUE)
    colnames(sm) <- "Intercept"
    null_Beta0 <- speedglm::speedglm.wfit(y = tmu,
                                         X = sm,
                                         intercept=FALSE,
                                         family = stats::quasibinomial(link = "logit"))[["coefficients"]]
    null_ret <- EM_mat(null_Beta0,
            sm,
            BF = BF,
            verbose = verbose
    )
    nb0 <- null_ret$data[[1]]$Beta
    if(parallel){
        Xl <- purrr::map(
                colnames(X),
                ~ X[, .x]
        )
        Beta0l <- purrr::map(
                Xl,
                ~ speedglm::speedglm.wfit(
                        y = as.vector(tmu),
                        X = cbind(.x, rep(1, NROW(.x))),
                        family = stats::quasibinomial(link = "logit")
                )[["coefficients"]]
        )
        ret_df <- purrr::pmap_dfr(list(x=Beta0l,
                                       y=Xl,
                                       z = colnames(X)
                                       ),
                                  function(x, y, z, BF, null_beta) {
                                      FGEM(x, magrittr::set_colnames(cbind(y, rep(1, NROW(y))), c(z, "Intercept")),
                                           BF = BF, null_beta = null_beta, verbose = FALSE
                                           )
                                  },BF = BF, null_beta = null_Beta0)
        return(ret_df)
    }
    else{
        purrr::map_dfr(colnames(X),
                fgem_col,
                X = X,
                null_Beta0 = nb0,
                Intercept_col = sm,
                tmu = tmu,
                BF = BF,
                vmessage = vmessage
                )
    }
}






##' Run FGEM with a dataframe as input
##' This version of FGEM takes a single dataframe with
##' both gene-level bayes factors and gene-level features
##' and returns a dataframe with effect size estimates
##' for every feature.
##' @param feat_df dataframe with data.  `feat_df` must have a column called `BF` specifying bayes factors.
##'   Any additional columns (besides an optional gene-name column that must be titled `Gene`) are treated as annotations.
##' @param prior_mean scalar between 0 and 1 specifying the starting
##' prior probability (this is for initializing the EM algorithm). If unspecified, it defaults to 0.02
##' @param Beta0 vector specifying starting guesses for effect sizes of annotations.
##' In almost all cases the default (NULL) should be kept
##' @param verbose Whether to print debug output
##' @param null_beta Instead of refitting the null model, you can pass the precomputed beta for the intercept-only model
##' @param ...
##' @export
FGEM_df <- function(feat_df, prior_mean = 0.02, Beta0 = NULL,verbose=FALSE,null_beta=NULL,...) {
    null_features <- "Intercept"
    vmessage <- verbose_message_factory(verbose)
    BF <- dplyr::pull(feat_df, BF)
    tmu <- (prior_mean * BF) / ((prior_mean * BF) + (1 - prior_mean))
    if (any(!null_features %in% colnames(feat_df))) {
        feat_df <- dplyr::mutate(feat_df, Intercept = 1)
        null_features <- "Intercept"
    }
    data_mat_df <- dplyr::select(feat_df, -Gene, -BF) %>% mutate(tmu = tmu)
    data_mat_df <- dplyr::mutate_all(data_mat_df, function(x) {
        if (all(is.na(x))) {
            return(rep(1, length(x)))
        } else {
            return(x)
        }
    })

    if (is.null(Beta0)) {
        vmessage("Settting Beta0...")
        Beta0 <- stats::coefficients(stats::glm(tmu ~ . + 0,
                                         data = data_mat_df,
                                         family = stats::quasibinomial(link = "logit")))
        vmessage(verbose, "Beta0: ", Beta0)
    }

    while (any(is.na(Beta0))) {
            bad_betas <- names(Beta0)[is.na(Beta0)]
            vmessage("Beta0 with NA:", paste0(bad_betas, collapse = ","))
            data_mat_df <- dplyr::select(data_mat_df, -one_of(bad_betas[1]))

            Beta0 <- stats::coefficients(stats::glm(tmu ~ . + 0,
                    data = data_mat_df,
                    family = stats::quasibinomial(link = "logit")
            ))
    }
    data_mat <- data.matrix(dplyr::select(data_mat_df, -tmu))
    vmessage("Starting FGEM:")
    ret <- FGEM(
            Beta0 = Beta0,
            feat_mat = data_mat,
            BF = BF,
            verbose = verbose,
            null_beta = null_beta,
            ...
    )
    return(ret)
}


##' @title Run FGEM with a matrix as input
##'
##'
##'This version of FGEM takes a single dataframe with
##' both gene-level bayes factors and gene-level features
##' and returns a dataframe with effect size estimates
##' for every feature.
##' @param Beta0 initial Beta0
##' @param feat_mat matrix of features
##' @param BF vector of bayes factors
##' @param verbose whether to give print intermediate results to the user
##' @param null_beta Instead of refitting the null model, you can pass the
##' precomputed beta for the intercept-only model
##' @param ... currently unused
##' @param null_features which features are "null" (deprecated)
##' @return results of model fit as a tibble
##' @export
FGEM <- function(Beta0, feat_mat, BF, verbose = FALSE, null_beta=NULL, ...) {
    null_features = "Intercept"
    nvmessage <- verbose_message_factory(verbose)
    nvmessage("Starting EM on full model ...")
    retdf <- EM_mat(Beta0, feat_mat, BF, verbose = verbose, ...) %>%
        dplyr::mutate(nac = NA)
    nvmessage("Starting EM on null model ...")
    if(is.null(null_beta)){
    retdf <- EM_mat(Beta0[null_features],
                    feat_mat[, null_features, drop = F],
                    BF = BF,verbose = verbose, ...) %>%
        dplyr::select(NullLogLik = LogLik) %>%
        dplyr::mutate(nac = NA) %>%
        dplyr::inner_join(retdf, by = "nac") %>%
        dplyr::select(-nac)
    }else{
        retdf <- mutate(retdf,
                NullLogLik = (FGEM_Logit_log_lik(
                        Beta = null_beta,
                        x = feat_mat[, null_features, drop = F],
                        B = BF
                ))
        ) %>% select(-nac)
    }
    retdf <- dplyr::mutate(retdf,
                           Chisq = -2 * (NullLogLik - LogLik),
                           pval = stats::pchisq(Chisq, df = ncol(feat_mat) - 1, lower.tail = F)
                           )
    return(retdf)
}
