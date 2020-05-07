


FGEM_Logit <- function(Beta, x, B, ctrl = list()) {
  uvec <- gen_u(Beta, x, B)
  if(inherits(x,"dgCMatrix")){
      return(speedglm::speedglm.wfit(
              y = as.vector(uvec),
              X = x,
              family = quasibinomial(link = "logit")
      )[["coefficients"]])
  }
  return(stats::coefficients(stats::glm.fit(
          x = x,
          y = uvec,
          family = quasibinomial(link = "logit"),
          control = ctrl
  )))
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
    cn <- colnames(feat_mat)
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
    return(pvec)
}

gen_u <- function(Beta, x, B) {
        p <- gen_p(Beta, x)
        return((p * B) / ((p * B) + (1 - p)))
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


fgem_col <- function(col_name,X,
                     Intercept_col,
                     tmu,
                     BF,
                     null_Beta0,
                     vmessage=message_factory(FALSE)) {
    vmessage("starting feature:", col_name)
    Xm <- magrittr::set_colnames(cbind(X[, col_name], Intercept_col), c(col_name, "Intercept"))
    stopifnot(ncol(Xm) == 2)
    Beta0 <- speedglm::speedglm.wfit(y = as.vector(tmu),
                                   X = Xm,
                                   family = quasibinomial(link = "logit"))[["coefficients"]]
    vmessage("Beta0", paste0(Beta0, collapse = ","))
    ret_fgem <- FGEM(
            Beta0 = Beta0,
            feat_mat = Xm,
            BF = BF,
            null_beta = null_Beta0,
            verbose = FALSE
    )
    return(ret_fgem)
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
    sm <- Matrix(rep(1,length(BF)),nrow=length(BF),ncol=1,sparse = TRUE)
    null_Beta0<- speedglm::speedglm.wfit(y = tmu,
                                         X = sm,
                                         intercept=FALSE,
                                         family = quasibinomial(link = "logit"))[["coefficients"]]
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
                        family = quasibinomial(link = "logit")
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
                                         family = quasibinomial(link = "logit")))
        vmessage(verbose, "Beta0: ", Beta0)
    }

    while (any(is.na(Beta0))) {
            bad_betas <- names(Beta0)[is.na(Beta0)]
            vmessage("Beta0 with NA:", paste0(bad_betas, collapse = ","))
            data_mat_df <- dplyr::select(data_mat_df, -one_of(bad_betas[1]))

            Beta0 <- stats::coefficients(stats::glm(tmu ~ . + 0,
                    data = data_mat_df,
                    family = quasibinomial(link = "logit")
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
