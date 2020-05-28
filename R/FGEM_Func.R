FGEM_Logit <- function(Beta, x, B, ctrl = list()) {
        uvec <- gen_u(Beta, x, B)

        fastglm::fastglmPure(
                x = as.matrix(x),
                y = uvec,
                family = stats::quasibinomial(link = "logit")
        )[["coefficients"]]
}

obj_factory <- function(feat_mat, BF, verbose) {
        tfm <- feat_mat
        tbf <- BF
        mvmessage <- verbose_message_factory(verbose)
        ret_f <- function(par, ...) {
            ret <- (-FGEM_Logit_log_lik(Beta = par, x = tfm, B = tbf))
            return(-ret)
        }
        return(ret_f)
}

prep_fgem <- function(X, BF, l2) {
        renv <- rlang::env(X = X, BF = BF, prec = l2)
        rl <- make_env_obj(inherits(X, "dgCmatrix"))
        rl[["env"]] <- renv
        return(rl)
}



##' @title fit cross-validated FGEM model
##'
##'
##' @param X Feature matrix
##' @param BF vector of bayes factors
##' @param alpha alpha (as in glmnet). Alpha  is the (scalar) proportion of
##' `lambda` applied to l1 penalization, while `1-alpha` is applied to l2
##' @param lambda vector of shrinkage parameters
##' @return
##' @export
cv_fgem <- function(X,
                    BF,
                    alpha = 1,
                    lambda = c(0, 5 * 10^(-(seq(6, 1, length.out = 75))), 1, 2)){


    idf <- tibble::tibble(BF = BF, X)
    v <- 10
    cv_X <- rsample::vfold_cv(X)
    cv_idf <- rsample::vfold_cv(idf, v = v, strata = BF)
    tivx <- cv_idf$splits[[1]]
    furrr::future_imap_dfr(cv_idf$splits, function(tivx,y) {
        tiv_df <- tibble::as_tibble(tivx)
        tBF <- tiv_df$BF
        tX <- tiv_df$X
        res_d <- fgem_bfgs(tX, tBF, lambda = lambda, alpha = alpha)
        tBeta <- coeff_mat(res_d$Beta)
        civx <- as.data.frame(tivx, data = "assessment")
        cv_lik <- apply(tBeta, 2, FGEM_Logit_log_lik, x = civx$X, B = civx$BF)
        dplyr::mutate(res_d, cv_lik = cv_lik,cv_i=y)
    })
}


coeff_mat <- function(Beta_l) {
        Beta <- simplify2array(purrr::map(Beta_l, "Beta"), higher = TRUE)
        dimnames(Beta) <- list(Beta_l[[1]]$feature_name, paste0("s", seq_along(Beta_l) - 1))
        return(Beta)
}




##' @title fit fgem using lbfgs
##'
##'
##' @param X feature matrix
##' @param BF vector of bayes factors
##' @param Beta0 initial guess (defaults to 0)
##' @param verbose (whether to give verbose output)
##' @param alpha alpha (as in glmnet). Alpha  is the (scalar) proportion of
##' `lambda` applied to l1 penalization, while `1-alpha` is applied to l2
##' @param lambda vector of shrinkage parameters
##' @param ... currently unused
##' @return tibble with results of fit
##' @export
fgem_bfgs <- function(X,
                      BF,
                      Beta0 = rep(0, ncol(X) + 1),
                      verbose = FALSE,
                      alpha = 1,
                      lambda = c(0,5 * 10^(-(seq(6, 1, length.out = 75))),1,2), ...) {


    iseq <- seq_along(lambda)
    stopifnot(length(alpha) == 1)
    reg_resl <- list()
    for (i in iseq) {
        l <- lambda[i] * NROW(X)
        l2 <- (1 - alpha) * l
        l1 <- (alpha) * l
        tto <- prep_fgem(X, BF, l2)
        algo <- dplyr::if_else(l1!=0,"LBFGS_LINESEARCH_BACKTRACKING","LBFGS_LINESEARCH_DEFAULT")
        pta <- proc.time()
        lbr <- lbfgs::lbfgs(
                call_eval = tto$lik,
                call_grad = tto$grad,
                environment = tto$env,
                vars = c(rep(0, ncol(X) + 1)),
                orthantwise_c = l1,
                linesearch_algorithm = algo,
                orthantwise_start = 1,
                orthantwise_end = ncol(X) + 1,
                invisible = dplyr::if_else(verbose, 0L, 1L)
        )
        lbr$time <- (proc.time() - pta)["elapsed"]
        lbr$l0n <- sum(lbr$par != 0)
        lbr$l1n <- sum(abs(lbr$par))
        lbr$l2n <- sum(lbr$par^2)
        lbr$lambda <- l
        lbr$alpha <- alpha
        reg_resl[[i]] <- tibble::tibble_row(
                                     Beta = list(tibble::tibble(
                                                             Beta = lbr$par,
                                                             feature_name = c("Intercept", colnames(X))
                                                         )),
                                     l0n = lbr$l0n,
                                     l1n = lbr$l1n,
                                     l2n = lbr$l2n,
                                     lik = -lbr$value,
                                     lambda = l,
                                     l1 = l1,
                                     l2 = l2,
                                     alpha = alpha,
                                     time = lbr$time,
                                     convergence = lbr$convergence
                                 )
        if (lbr$l0n == 1) {
                break
        }
    }
    return(dplyr::bind_rows(reg_resl))
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

log1pexp <- function(x)
{
  indx <- .bincode(x,
                   c(-Inf, -37, 18, 33.3, Inf),
                   right = TRUE,
                   include.lowest = TRUE)

  kk <- which(indx==1)
  if( length(kk) ){  x[kk] <- exp(x[kk])  }

  kk <- which(indx==2)
  if( length(kk) ){  x[kk] <- log1p( exp(x[kk]) ) }

  kk <- which(indx==3)
  if( length(kk) ){  x[kk] <- x[kk] + exp(-x[kk]) }
  return(x)
}

logsum <- function(l1, l2) {
        pmax(l1, l2) + log1p(exp(-abs(l1 - l2)))
}

log_FGEM_log_lik <- function(Beta, x, logBF) {
        xb <- c((x %*% Beta[-1]) + Beta[1])
        sum(logsum(-xb, logBF) + plogis(xb, log = TRUE))
}


FGEM_Logit_log_lik <- function(Beta, x, B) {
        pvec <- gen_p(Beta = Beta, x = x)
        return(sum(log(pvec * B + (1 - pvec))))
}



gen_p <- function(Beta, x, log = FALSE) {

        stats::plogis((x %*% Beta[-1]) + Beta[1], log = log)
}

gen_u <- function(Beta, x, B, log=FALSE) {
    if(!log){
        p <- gen_p(Beta, x, log = log)
        return(as.vector((p * B) / ((p * B) + (1 - p))))
    }else{
        return(as.vector(-log1pexp(- (x %*% Beta[-1] + Beta[1]) - log(B))))
    }
}


log_log_u <- function(Beta, x, logBF) {
    c(-log1pexp(-(x %*% Beta[-1] + Beta[1]) - logBF))
}


##' @title Predict posterior from fgem fit
##'
##'
##' @param fit result of call to `fgem`
##' @param x annotation matrix (colnames should correspond to `feature_name` column of fit$data[[1]])
##' @param BF Bayes Factor
##' @param lambda value of lambda to use in prediction
##' @param log return log posterior
##' @return vector of length length(BF) with posterior probabilities
##' @export
predict_fgem <- function(fit, x, BF, lambda, log = FALSE) {
        beta_df <- dplyr::pull(dplyr::filter(fit, lambda == lambda), Beta)[[1]]
        beta <- magrittr::set_names(beta_df$Beta, beta_df$feature_name)
        beta <- beta[order(beta)]
        magrittr::set_names(gen_u(Beta = beta, x = x[, names(beta)], B = BF, log = log), rownames(x))
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
        pvec <- 1 / (1 + exp(-(feat_mat %*% Beta[-1] + Beta[1])))
        return(mean(pvec))
}


fgem_null <- function(BF, prior = 0.02, tmu = prior_mean(BF, prior)) {
        fgem_bfgs(matrix(0, length(BF), 0, dimnames = list(names(BF), character())), BF, lambda = 0)$Beta[[1]]$Beta
}

fgem_null_lik <- function(BF) {
        fgem_bfgs(matrix(0, length(BF), 0, dimnames = list(names(BF), character())), BF, lambda = 0)$lik
}


prior_mean <- function(BF, prior = 0.02, log = FALSE) {
        if (!log) {
                (prior * BF) / ((prior * BF) + (1 - prior))
        } else {
                -log1p((1 - prior) / (prior * BF))
        }
}



guess_beta0 <- function(X, BF, prior = 0.02, tmu = prior_mean(BF, prior),add_intercept=TRUE) {
    if (!"Intercept" %in% colnames(X)) {
        if(add_intercept){
            inter_x <- matrix(1, nrow = NROW(X), ncol = 1, dimnames = list(rownames(X), "Intercept"))
            x <- cbind(X, inter_x)
        }
    }
    fastglm::fastglmPure(
                 y = tmu,
                 x = X,
                 family = stats::quasibinomial(link = "logit")
             )[["coefficients"]]
}


##' @title Maximize FGEM marginalized likelihood
##'
##'
##' @param x gene-level feature set
##' @param BF gene-level bayes factors
##' @param prior (scalar?) prior probability z=1
##' @param tmu (gene-level) prior probability z=1 (takes precedence over `prior` if both are specified)
##' @param Beta0 initial guess for enrichment
##' @param verbose logical scalar.  Print out lots of debug info?
##' @param ...
##' @return dataframe with nested dataframe with multivariate effect-size estimates, and FGEM likelihood
##' @export
fgem_lik <- function(x,BF,
                     prior=0.02,
                     tmu = prior_mean(BF, prior),
                     add_intercept=TRUE,
                     Beta0 = guess_beta0(x, BF, prior, tmu, add_intercept),
                     verbose = FALSE,...){
    stopifnot(NROW(x) == length(BF))
    vmessage <- verbose_message_factory(verbose)
    vmessage("Beta0", paste0(Beta0, collapse = ","))
    if (!"Intercept" %in% colnames(x)) {
        if(add_intercept){
            inter_x <- matrix(1, nrow = NROW(x), ncol = 1, dimnames = list(rownames(x), "Intercept"))
            x <- cbind(x, inter_x)
        }
    }
    EM_mat(Beta0, x, BF, verbose = verbose, ...)
}






##' @title Maximize FGEM marginalized likelihood
##'
##'
##' @param x gene-level feature set
##' @param BF gene-level bayes factors
##' @param prior (scalar?) prior probability z=1
##' @param tmu (gene-level) prior probability z=1 (takes precedence over `prior` if both are specified)
##' @param Beta0 initial guess for enrichment
##' @param verbose logical scalar.  Print out lots of debug info?
##' @param jaccard_cutoff cutoff for removing correlated features
##' @param qval_cutoff adjusted p-value (by method="fdr") for the stopping criteria for forward selection
##' @return dataframe with nested dataframe with multivariate effect-size estimates, and FGEM likelihood
##' @export
forward_select_fgem_lik <- function(x,BF,
                                    prior=0.02,
                                    tmu = prior_mean(BF, prior),
                                    Beta0 = guess_beta0(x, BF, prior, tmu),
                                    verbose = FALSE,
                                    jaccard_cutoff=0.1,
                                    qval_cutoff=0.05,
                                    ...){
    stopifnot(NROW(x) == length(BF))
    vmessage <- verbose_message_factory(verbose)
    vmessage("Beta0", paste0(Beta0, collapse = ","))
    if (!"Intercept" %in% colnames(x)) {
        inter_x <- matrix(1, nrow = NROW(x), ncol = 1, dimnames = list(rownames(x), "Intercept"))
        x <- cbind(x, inter_x)
    }else{
        inter_x <- x[, "Intercept", drop = FALSE]
    }
    remaining_terms <- colnames(x)
    current_model <- c("Intercept")
    null_model <- c("Intercept")
    correlated_features <- character()
    null_results <- fgem_lik(x[, null_model, drop = FALSE], BF, prior, tmu, add_intercept = FALSE)
    remaining_terms <- remaining_terms[!remaining_terms %in% current_model]
    if (length(remaining_terms) > 0) {
            axm <- purrr::array_branch(x[, remaining_terms,drop=FALSE], 2)
            jacc_fun <- function(x, y) {
                    sum((x & y)) / sum(x | y)
            }
            feat_dist <- outer(axm, axm, FUN = Vectorize(jacc_fun))
    }
    while (length(remaining_terms) > 0) {
        NullLik <- null_results$LogLik
        fit_fun <- function(i, remaining_terms, current_model, BF, prior, tmu, NullLik, x) {
            lterm <- remaining_terms[i]
            dplyr::mutate(fgem::fgem_lik(x[,
                                           unique(c(lterm, current_model)),
                                           drop = FALSE
                                           ],
                                         BF,
                                         prior,
                                         tmu,
                                         add_intercept = FALSE
                                         ),
                          term = lterm,
                          Chisq = -2 * (NullLik - LogLik),
                          pval = stats::pchisq(Chisq, df = 1, lower.tail = F)
                          )
        }
        single_df <- purrr::map_dfr(
                seq_along(remaining_terms),
                fit_fun,
                remaining_terms = remaining_terms,
                current_model = current_model,
                BF = BF,
                prior = prior,
                tmu = tmu,
                NullLik = NullLik,
                x = x
                ) %>%
            mutate(qval = p.adjust(pval, method = "fdr"))
        bad_terms <- filter(single_df, qval >= qval_cutoff) %>% pull(term)
        single_df <- filter(single_df, qval < qval_cutoff)
        remaining_terms <- remaining_terms[!remaining_terms %in% bad_terms]
        if (nrow(single_df) > 0) {
            null_results <- filter(single_df, pval == min(pval)) %>%
                slice(1)
            next_term <- null_results$term
            current_model <- unique(c(current_model, next_term))
        }
        remaining_terms <- remaining_terms[!remaining_terms %in% current_model]
        not_int <- current_model[current_model != "Intercept"]
        t_feat_dist <- feat_dist[not_int,,drop=FALSE]
        correlated_features <- unique(c(
            correlated_features,
            colnames(t_feat_dist)[which(t_feat_dist > jaccard_cutoff, arr.ind = TRUE)[, "col"]]
        ))
        remaining_terms <- remaining_terms[!remaining_terms %in% correlated_features]
    }
    return(null_results)
}




##' @title Run FGEM
##'
##'
##' @param X matrix of gene-level annotations
##' @param BF vector of gene-level bayes factors
##' @param prior prior probability a gene is causal (before annotations)
##' @param tmu precomputed gene-level prior
##' @param null_beta model to test against for likelihood ratio test
##' @param add_intercept add intercept to X if it isn't there?
##' @param verbose whether to use verbose output
##' @return tibble with fgem results
##' @export
fgem <- function(x,
                 BF,
                 prior=0.02,
                 tmu = prior_mean(BF, prior),
                 null_beta=fgem_null(BF, tmu),
                 add_intercept=TRUE,
                 verbose = FALSE) {
    stopifnot(NROW(x) == length(BF))
    vmessage <- verbose_message_factory(verbose)


    if (!"Intercept" %in% colnames(x)) {
        stopifnot(add_intercept)
        inter_x <- matrix(1, nrow = NROW(x), ncol = 1, dimnames = list(rownames(x), "Intercept"))
        x <- cbind(x, inter_x)
    }else{
        inter_x <- x[, "Intercept"]
    }
    Beta0 <- guess_beta0(x, BF, prior, tmu)
    vmessage("Beta0", paste0(Beta0, collapse = ","))
    ret_fgem <- FGEM(
            Beta0 = Beta0,
            feat_mat = x,
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
                              add_intercept = TRUE,
                              csr=FALSE) {

    stopifnot(
            all.equal(length(rowname_vals), length(colname_vals)),
            all.equal(length(rowname_vals), length(values))
        )

    bad_vals <- (!(rowname_vals %in% total_rownames)) | (!(colname_vals %in% total_colnames))
    ## if (any(bad_vals)) {
    ##     warning("rowname_vals and/or colname_vals not found in total_rownames/total_colnames, dropping: ",
    ##             sum(bad_vals))
    ## }


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
            i = row_id[rowname_vals],
            j = col_id[colname_vals],
            dims = c(length(total_rownames), length(total_colnames)),
            x = values, dimnames = list(total_rownames, total_colnames)
    )
    if (csr) {
        rstan::extract_sparse_parts(ret_m)
    }
    else{
        as.matrix(ret_m)
    }
}



##' Run FGEM with a dataframe as input
##' This version of FGEM takes a single dataframe with
##' both gene-level bayes factors and gene-level features
##' and returns a dataframe with effect size estimates
##' for every feature (run as univariate)
##' @param X feature matrix (must have column names)
##' @param BF vector of bayes factors
##' @param verbose Whether to print debug output
##' @export
fgem_marginal <- function(X,BF,prior_mean = 0.02, ...) {
    null_lik <-fgem_null_lik(BF)
    results <- marginal_fgem_fit_bfgs(X,BF)
    Chisq <- -2 * (null_lik - (results$lik))
    tibble::tibble(
            pval = stats::pchisq(Chisq, df = 1, lower.tail = F),
            iter = results$iter,
            coeff = results$coeff[2, ],
            Intercept = results$coeff[1, ],
            feat_sum = results$feat_sum,
            feature_name = results$feature_name,
            qval = p.adjust(pval, method = "fdr")
    )
}



join_long_wide <- function(long_df, wide_df, key="feature_name", value="value", by = "Gene") {

    ldv <- long_df[[value]]
    if (is.null(ldv))
        ldv <- rep(1.0, nrow(long_df))


    ixm <- trip2sparseMatrix(
        rowname_vals = long_df[[by]],
        colname_vals = long_df[[key]],
        values = ldv,
        total_rownames = unique(wide_df[[by]]),
        total_colnames = unique(long_df[[key]]),
        add_intercept = FALSE
    )
    dplyr::inner_join(wide_df, tibble::as_tibble(ixm, rownames = by), by = by)
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
