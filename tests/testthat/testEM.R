context("simulation")

test_that("FGEM works on simulated data",{
  g <- 2000
  sz <- 200
  X <- cbind(rbinom(n = g,size=1,prob=0.5),
             rbinom(n = g,size=1,prob=0.1),
             rbinom(n = g,size=1,prob=0.4))
  storage.mode(X) <- "numeric"
  eff <- c(-3,0,5,5)  
  fn <- c("Intercept","V1","V2","V3")
  true_value <- tibble::tibble(feature_name=fn,true_Beta=eff)  
  pr <- fgem:::gen_p(eff,X)
  y <- rbinom(g,size=1,prob=pr)
  q <- dplyr::if_else(y==1,0.5,runif(g))
  xg <- rbinom(g,size=sz,prob = q)
  pxz1 <- dbinom(x = xg,size = sz,prob=0.5)
  pxz2 <- beta((xg+1),sz-xg+1)*choose(sz,xg)
  bf <- log(pxz1)-log(pxz2)
  lambda <- c(2, 1, 5 * 10^(-(seq(1, 
    9, length.out = 75))), 0)
  tr <- fgem:::fgem_bfgs(X,bf,log_BF = TRUE,alpha=0.9,lambda=lambda)
  progressr::with_progress({
  cvtr <- fgem:::cv_fgem(X,bf,alpha=0.9,lambda=lambda,log_BF = TRUE)
  })
  summ_cv <- fgem:::summarise_cv_lik(cvtr)

comp_df <-   tibble(l1=tr$lambda,Beta=map2(rev(summ_cv$Beta),tr$Beta,function(cv_df,ncv_df){
    bind_rows(transmute(cv_df,Beta=Beta,feature_name=feature_name,model="cv"),
               transmute(ncv_df,Beta=Beta,feature_name=feature_name,model="ncv"))
  }),lik=summ_cv$cv_sum) %>% unnest(Beta)
  
  unnest(summ_cv,Beta) %>% 
    filter(group_l1 > 0) %>% 
    inner_join(true_value) %>% 
    ggplot(aes(x = group_l1, y = Beta, col = feature_name))+
    geom_line()+
    scale_x_log10()+
    geom_hline(aes(yintercept = true_Beta, col = feature_name, linetype = feature_name))
  filter(summ_cv, group_l1 > 0) %>% 
    ggplot(aes(x=group_l1,y=cv_sum))+geom_point()+scale_x_log10()
  
  unnest(tr,Beta) %>% 
    inner_join(true_value) %>% 
    ggplot(aes(x=lambda,y=Beta,col=feature_name))+
    geom_line()+
    scale_x_log10()+
    geom_hline(aes(yintercept=true_Beta,col=feature_name,linetype=feature_name))
  preds <- purrr::map(tr$Beta,~fgem:::gen_u(Beta = .x$Beta,x = X,B = bf))

tr <- mutate(tr,Beta=purrr::map(Beta,
                           ~mutate(.x,grad=fgem:::fgem_grad_stan(Beta,X,bf))))
  
  hess_mats <- purrr::map(tr$Beta,
                          ~fgem:::fgem_hess_stan(.x$Beta,X,bf))
})




