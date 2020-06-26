prep_fgem <- function(X, BF, l2, log_BF=FALSE) {
        renv <- rlang::env(X = X, BF = BF, l2 = l2, log_BF = log_BF)
        rl <- make_env_obj(inherits(X, "dgCMatrix"))
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
                    lambda = c(2,1, 5 * 10^(-(seq(1, 6, length.out = 75))), 0),
                    stratify_BF=TRUE,
                    log_BF=FALSE,
                    ...){

    v <- 10

    if (!inherits(X, "dgCMatrix")) {
      fgls <- function(par, x, BF, l2 = 0, neg = FALSE, log_BF = FALSE){
        fgem_lik_stan(par,x,BF,l2,neg,log_BF)
      }
            idf <- tibble::tibble(BF = BF, X)
            strat_fun <- function(tivx, y) {
                tiv_df <- tibble::as_tibble(tivx)
                tBF <- tiv_df$BF
                tX <- tiv_df$X
                res_d <- fgem_bfgs(tX, tBF, lambda = lambda, alpha = alpha, log_BF = log_BF)
                tBeta <- coeff_mat(res_d$Beta)
                civx <- as.data.frame(tivx, data = "assessment")
                aX <- civx$X
                aBF <- civx$BF
                cv_lik <- apply(X=tBeta,MARGIN= 2,FUN= fgls, x = aX, BF = aBF, log_BF = log_BF)
                p(message = sprintf("cv-%d", as.integer(y)))
                dplyr::mutate(res_d,
                              cv_lik = cv_lik,
                              cv_i = y,
                              group_l1 = round(l1 / num_X, digits = 10),
                              group_l2 = round(l2 / num_X, digits = 10))
            }
    }else{
      idf <- dplyr::mutate(tibble::tibble(BF = BF), idx = 1:dplyr::n())
      fgls <- function(par, x, BF, l2 = 0, neg = FALSE, log_BF = FALSE){
        sp_fgem_lik_stan(par,x,BF,l2,neg,log_BF)
      }
      strat_fun <- function(tivx, y) {
        tiv_df <- tibble::as_tibble(tivx)
        tBF <- tiv_df$BF
        tX <- X[tiv_df$idx, ,drop=FALSE]
        res_d <- fgem_bfgs(tX, tBF, lambda = lambda, alpha = alpha, log_BF = log_BF)
        tBeta <- coeff_mat(res_d$Beta)
        civx <- as.data.frame(tivx, data = "assessment")
        aX <- X[civx$idx,,drop=FALSE]
        aBF <- civx$BF
        cv_lik <- apply(X=tBeta, MARGIN = 2,FUN = fgls , x = aX, BF = aBF, log_BF = log_BF)
        p(message = sprintf("cv-%d", as.integer(y)))
            dplyr::mutate(res_d,
                          cv_lik = cv_lik,
                          cv_i = y,
                          group_l1 = round(l1 / num_X, digits = 10),
                          group_l2 = round(l2 / num_X, digits = 10))
        }
    }
    if (stratify_BF) {
        cv_idf <- rsample::vfold_cv(idf, v = v, strata = BF)
    } else {
        cv_idf <- rsample::vfold_cv(idf, v = v)
    }
    p <- progressr::progressor(along = cv_idf$splits)
    furrr::future_imap_dfr(cv_idf$splits, strat_fun)
}




coeff_mat <- function(Beta_l) {
        Beta <- simplify2array(purrr::map(Beta_l, "Beta"), higher = TRUE)
        dimnames(Beta) <- list(Beta_l[[1]]$feature_name, paste0("s", seq_along(Beta_l) - 1))
        return(Beta)
}



fgem_elasticnet <- function(X, BF, Beta0=rep(0, NCOL(X) + 1), alpha=1,lambda=0, verbose=FALSE,log_BF=FALSE,grad=FALSE,hess=FALSE,check_conv=TRUE, ...) {

    if(!inherits(X,"dgCMatrix")){
        X <- as.matrix(X)
    }

    l <- lambda * NROW(X)
    l2 <- (1 - alpha) * l
    l1 <- (alpha) * l
    tto <- prep_fgem(X, BF, l2, log_BF)
    al <- rlang::list2(...)
    pta <- proc.time()
    lbr <- rlang::exec(lbfgs::lbfgs,
            !!!al,
            call_eval = tto$lik,
            call_grad = tto$grad,
            environment = tto$env,
            vars = c(rep(0, NCOL(X) + 1)),
            orthantwise_c = l1,
            orthantwise_start = 1,
            orthantwise_end = NCOL(X) + 1,
            invisible = dplyr::if_else(verbose, 0L, 1L)
            )
    lbr$time <- (proc.time() - pta)["elapsed"]
    lbr$l0n <- sum(lbr$par != 0)
    lbr$l1n <- sum(abs(lbr$par))
    lbr$l2n <- sum(lbr$par^2)
    lbr$lambda <- l
    lbr$alpha <- alpha
    gradient <- NULL
    if (lbr$convergence == 0 && check_conv) {
        gradient <- fgem_grad(par = lbr$par, X = X, BF = BF, l2 = l2, l1 = l1, log_BF = log_BF)
        sg <- sign(gradient)


    }


    cn <- c("Intercept", colnames(X) %||% paste0("V", seq_len(NCOL(X))))
    if (NCOL(X) == 0) {
            cn <- character("Intercept")
    }
    BL <- list(tibble::tibble(
            Beta = lbr$par,
            feature_name = cn,
    ))
    hessian <- NULL
    if (grad) {
      if(is.null(gradient)){
        gradient <- fgem_grad(par = lbr$par, X = X, BF = BF, l2 = l2, l1 = l1, log_BF = log_BF)
      }
      BL[[1]]$gradient <- gradient
    }
    if (hess) {
            hessian <- fgem_hess(par = lbr$par, X = X, BF = BF, l2 = l2, l1 = l1, log_BF = log_BF)
            BL[[1]]$hessian <- hessian
    }

    rdf <- tibble::tibble_row(
            Beta = BL,
            l0n = lbr$l0n,
            l1n = lbr$l1n,
            l2n = lbr$l2n,
            lik = -lbr$value,
            lambda = l,
            l1 = l1,
            l2 = l2,
            alpha = alpha,
            time = lbr$time,
            num_X = NROW(X),
            convergence = lbr$convergence
    )

    return(rdf)
}


##' @title Helper to get resolved futures
##'
##'
##' @param x an object
##' @return if x inherits Future, it returns
##' the value of x if it is resolved, or NULL
##' if is not.  If x is not a Future, it returns
##' x
##'
resolved_or_null <- function(x) {
        if (inherits(x, "Future")) {
                if (future::resolved(x)) {
                        return(future::value(x))
                } else {
                        return(NULL)
                }
        } else {
                return(x)
        }
}


##' @title retrieve only resolved results from list of futures
##'
##'
##' @param x a list of futures
##' @param atleast the min number of resolved futures to return
##' @return a list of values of length <= length(x)
##'
retrieve_resolved <- function(x, atleast = 1) {
        x <- wait_on_n(x, atleast)
        purrr::map(x, resolved_or_null)
}



##' @title Helper to get resolved futures
##'
##'
##' @param x an object
##' @return if x inherits Future, it returns
##' the value of x if it is resolved, or itself
##' if is not.  If x is not a Future, it returns
##' itself
partial_resolve <- function(x) {
        if (inherits(x, "Future")) {
                if (future::resolved(x)) {
                        return(future::value(x))
                } else {
                        return(x)
                }
        } else {
                return(x)
        }
}

value_wrapper <- function(x) {
        if (inherits(x, "Future")) {
                      return(future::value(x))
              }
        return(x)
}

num_resolved <- function(x) {
    sum(purrr::map_lgl(x, future::resolved))
    
}

##' @title take a list of futures and block until at least `n` future is resolved
##'
##'
##' @param x container of futures
##' @return x
wait_on_n <- function(x, n=1,pause = 0.02) {
    n <- pmin(n, length(x))
    rate <- purrr::rate_delay(pause)
    ready_n <- num_resolved(x)
    while (ready_n < n) {
        purrr::rate_sleep(rate)
        ready_n <- num_resolved(x)
    }
    return(x)
}

last_completed <- function(x) {
        rr <- wait_on_n(x, 1) %>%
                retrieve_resolved()
        return(rr[[length(rr)]])
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
                      Beta0 = c(0, rep(0, ncol(X))),
                      verbose = FALSE,
                      alpha = 1,
                      lambda = c(2,1, 5 * 10^(-(seq(1, 6, length.out = 75))), 0),
                      add_lambda=FALSE,
                      progress = FALSE,
                      log_BF=FALSE,
                      ...){

    lambda <- sort(lambda, decreasing = TRUE)
    iseq <- seq_along(lambda)
    stopifnot(length(alpha) == 1)
    arl <- rlang::list2(...)
    reg_resl <- list()
    old_max_lambda <- lambda[1]
    max_lambda <- old_max_lambda
    reg_resl <- purrr::prepend(reg_resl, list(fgem_elasticnet(X,
            BF,
            Beta0,
            alpha,
            max_lambda,
            verbose,
            log_BF=log_BF,
            ... = ...
    )))
    Beta0 <- reg_resl[[length(reg_resl)]]$Beta[[1]]$Beta

    tl0n <- reg_resl[[length(reg_resl)]]$l0n
    if(add_lambda){
        while (tl0n > 1) {
                max_lambda <- max_lambda * 2
                lambda <- purrr::prepend(lambda, max_lambda)
                reg_resl <- purrr::prepend(reg_resl, list(fgem_elasticnet(X,
                        BF,
                        Beta0,
                        alpha,
                        max_lambda,
                        verbose,
                        log_BF=log_BF,
                        ... = ...
                )))
                Beta0 <- reg_resl[[length(reg_resl)]]$Beta[[1]]$Beta
                tl0n <- reg_resl[[length(reg_resl)]]$l0n
        }
    }else{
        if (tl0n > 1 && alpha>0) {
            warning("largest lambda provided does not shrink to 0, rerun with `add_lambda=TRUE` to ensure complete regularization path")
        }
    }
    rest_lambda <- lambda[lambda < old_max_lambda]
    if (progress) {
            pb <- utils::txtProgressBar(style=3)
    }
    pc <- 0
    for (i in seq_along(rest_lambda)) {
        l <- rest_lambda[i]

            Beta0 <- last_completed(reg_resl)$Beta[[1]]$Beta

            reg_resl <- append(reg_resl, future::future(fgem_elasticnet(X,
                    BF,
                    Beta0,
                    alpha,
                    l,
                    verbose,
                    log_BF=log_BF,
                    ... = ...
                    )))
        opc <- pc
        npc <- num_resolved(reg_resl)
        nc <- npc-pc
        pc <- npc

        if (progress) {
                utils::setTxtProgressBar(pb, npc / length(rest_lambda))
        }
        Sys.sleep(.15)
    }
    return(purrr::map_dfr(reg_resl, value_wrapper))
}

logsum <- function(l1, l2) {
        pmax(l1, l2) + log1p(exp(-abs(l1 - l2)))
}

logspace_sub <- function (logx, logy)
{
    # log(exp(logx) - exp(logy)), avoiding unnecessary floating point error
    dxy <- logx - logy
    # log(2) looks like best breakpoint
    logx + ifelse(dxy < 0.693147180559945, log(-expm1(-dxy)), log1p(-exp(-dxy)))
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
    if (NCOL(Beta) == 1) {
        xb <- (x %*% Beta[-1]) + Beta[1]
        if (is.null(ncol(xb))) {
            xb <- as.matrix(xb)
        }
    }else{
        xb <- t(t(x %*% Beta[-1, ]) + Beta[1, ])
    }
    apply(xb, 2, stats::plogis, log = log)
}

gen_u <- function(Beta, x, B, log = FALSE) {
    if(!log){
        p <- gen_p(Beta, x, log = log)
        return(as.vector((p * B) / ((p * B) + (1 - p))))
    }else{
        return(-log_1p_exp(-(x %*% Beta[-1] + Beta[1]) - log(B)))
    }
}


log_log_u <- function(Beta, x, logBF) {
    c(-log_1p_exp(-(x %*% Beta[-1] + Beta[1]) - logBF))
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

fgem_null <- function(BF, log_BF = FALSE, ...) {
    stats::optimize(fgem_lik, lower = -10, upper = 100,  X = matrix(0, nrow = length(BF), ncol = 0), BF = BF, log_BF = log_BF, maximum = TRUE)$maximum
}


fgem_null_lik <- function(BF,log_BF=FALSE) {
    stats::optimize(fgem_lik, lower = -10, upper = 100,  X = matrix(0, nrow = length(BF), ncol = 0), BF = BF,log_BF = log_BF, maximum = TRUE)$objective

}

fgem_null_fit <- function(BF, log_BF = FALSE, ...) {
    pta <- proc.time()
    lbf <- stats::optimize(fgem_lik_stan,
                    lower = -10,
                    upper = 100, X = matrix(0, nrow = length(BF), ncol = 0), BF = BF,
                    log_BF=log_BF,
                    maximum = TRUE)
    ptb <- (proc.time() - pta)["elapsed"]
    return(tibble::tibble_row(
                       Beta = list(tibble::tibble(
                                               Beta = lbf$maximum,
                                               feature_name = "Intercept"
                                           )),
                       l0n = 1,
                       l1n = abs(lbf$maximum),
                       l2n = lbf$maximum^2,
                       lik = lbf$objective,
                       lambda = Inf,
                       l1 = Inf,
                       l2 = 0 * Inf,
                       alpha = 1,
                       time = ptb,
                       num_X = length(BF),
                       convergence = 0
                   ))
}


prior_mean <- function(BF, prior = 0.02, log = FALSE) {
    if (!log) {
        (prior * BF) / ((prior * BF) + (1 - prior))
    } else {
        -log1p((1 - prior) / (prior * BF))
    }
}


##' @title Run FGEM
##'
##'
##' @param X matrix of gene-level annotations
##' @param BF vector of gene-level bayes factors
##' @param verbose whether to use verbose output
##' @return tibble with fgem results
##' @export
fgem <- function(x,
                 BF,
                 alpha=1,
                 lambda = c(2,1, 5 * 10^(-(seq(1, 6, length.out = 75))), 0),
                 verbose = FALSE,
                 ...) {
    stopifnot(NROW(x) == length(BF))
    vmessage <- verbose_message_factory(verbose)

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





##' Run FGEM with a dataframe as input
##' This version of FGEM takes a single dataframe with
##' both gene-level bayes factors and gene-level features
##' and returns a dataframe with effect size estimates
##' for every feature (run as univariate)
##' @param X feature matrix (must have column names)
##' @param BF vector of bayes factors
##' @param verbose Whether to print debug output
##' @export
fgem_marginal <- function(X,BF,prior_mean = 0.02,epsilon=1e-06,max_iter=150,parallel=FALSE,log_BF=FALSE, ...) {
    null_lik <-fgem_null_lik(BF,log_BF=log_BF)
    if (parallel) {
#        X <- as.matrix(X)
        results <- future.apply::future_apply(X, 2, fgem_elasticnet, BF = BF, log_BF = log_BF, ... = ...)
    } else {
        results <- apply(X, 2, fgem_elasticnet, BF = BF, log_BF = log_BF, ... = ...)
    }
    if (!is.null(colnames(X))) {
            results <- purrr::map2_dfr(results, colnames(X), function(x, y) {
                    x$Beta[[1]]$feature_name[2] <- y
                    return(x)
            })
    } else {
            results <- dplyr::bind_rows(results)
    }
    Chisq <- -2 * (null_lik - (results$lik))
    dplyr::mutate(results,
            pval = stats::pchisq(Chisq, df = 1, lower.tail = F),
            qval = p.adjust(dplyr::if_else(convergence == 0, pval, 1), method = "fdr"),
            sum_X = colSums(X)
    )
}
