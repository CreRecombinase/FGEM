#include "fgem.h"
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(unwindProtect)]]
//[[Rcpp::depends(BH)]]

typedef Rcpp::NumericVector (*efuncPtr)(SEXP,SEXP);


//' Gradient for FGEM likelihood
//' This is an attempt to use stan's AD features to optimize the FGEM likelihood
//' @export
//[[Rcpp::export]]
Eigen::ArrayXd fgem_grad_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double l2=0.0,const double l1=0.0,const bool log_BF=false){
  using Mt = Eigen::MatrixXd;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  const size_t p=par.size();
  double fx=0;
  if(l2<0.0 or l1<0.0){
    Rcpp::stop("penalties (l2/L2 and l1) must be non-negative");
  }

  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }
  if(l1>0.0){
    if(!log_BF){
      fgem_lik_l1<Mt> f(X,BF,l2,l1);
      stan::math::gradient(f,param,fx,grad_fx);
    }else{
      log_fgem_lik_l1<Mt> f(X,BF,l2,l1);
      stan::math::gradient(f,param,fx,grad_fx);
    }
    for(int i=1; i<p; i++){
      if(par[i]>0){
        grad_fx[i]+=l1;
      }
      if(par[i]<0){
        grad_fx[i]-=l1;
      }
      if(par[i]==0){
        if((grad_fx[i]+l1) < 0){
          grad_fx[i]+=l1;
        }else{
          if((grad_fx[i]-l1) > 0){
          grad_fx[i]-=l1;
          }else{
            grad_fx[i]=0;
          }
        }
      }
    }
  }else{
    if(!log_BF){
      fgem_lik<Mt> f(X,BF,l2);
      stan::math::gradient(f,param,fx,grad_fx);
    }else{
      log_fgem_lik<Mt> f(X,BF,l2);
      stan::math::gradient(f,param,fx,grad_fx);
    }
  }
  return(grad_fx.array());

}


//' Gradient for FGEM likelihood
//' This is an attempt to use stan's AD features to optimize the FGEM likelihood
//' @export
//[[Rcpp::export]]
Eigen::ArrayXd sp_fgem_grad_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::SparseMatrix<double>> X, const Eigen::Map<Eigen::ArrayXd> BF,const double l2=0.0,const double l1=0.0,const bool log_BF=false){
  using Mt = Eigen::SparseMatrix<double>;
Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;

  if(l2<0 or l1<0){
    Rcpp::stop("penalties (l2/L2 and l1) must be non-negative");
  }
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }


  double fx=0;
  const size_t p=par.size();
  if(l1!=0){
    if(!log_BF){
      fgem_lik_l1<Mt> f(X,BF,l2,l1);
      stan::math::gradient(f,param,fx,grad_fx);
    }else{
      log_fgem_lik_l1<Mt> f(X,BF,l2,l1);
      stan::math::gradient(f,param,fx,grad_fx);
    }
    for(int i=1; i<p; i++){
      if(par[i]>0){
        grad_fx[i]+=l1;
      }
      if(par[i]<0){
        grad_fx[i]-=l1;
      }
      if(par[i]==0){
        if((grad_fx[i]+l1) < 0){
          grad_fx[i]+=l1;
        }else{
          if((grad_fx[i]-l1) > 0){
            grad_fx[i]-=l1;
          }else{
            grad_fx[i]=0;
          }
        }
      }
    }
  }else{
    if(!log_BF){
      fgem_lik<Mt> f(X,BF,l2);
      stan::math::gradient(f,param,fx,grad_fx);
    }else{
      log_fgem_lik<Mt> f(X,BF,l2);
      stan::math::gradient(f,param,fx,grad_fx);
    }
  }
  return(grad_fx.array());
}




//' evd_dnorm_hess_stan
//'
//' This is an attempt to use stan's AD features to calculate a hessian
//' for the RSSp likelihood
//'
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd sp_fgem_hess_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::SparseMatrix<double>> X, const Eigen::Map<Eigen::ArrayXd> BF,const double l2=0.0,const double l1=0.0,const bool neg=false,const bool log_BF=false){

  using Mt = Eigen::SparseMatrix<double>;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_fx;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }
  if(l2<0 or l1<0){
    Rcpp::stop("penalties (l2/L2 and l1) must be non-negative");
  }

  double fx=0;

  if(l1!=0){
    if(!log_BF){
      fgem_lik_l1<Mt> f(X,BF,l2,l1);
      stan::math::hessian(f,param,fx,grad_fx,H_fx);
    }else{
      log_fgem_lik_l1<Mt> f(X,BF,l2,l1);
      stan::math::hessian(f,param,fx,grad_fx,H_fx);
    }
  }else{
    if(!log_BF){
      fgem_lik<Mt> f(X,BF,l2);
      stan::math::hessian(f,param,fx,grad_fx,H_fx);
    }else{
      log_fgem_lik<Mt> f(X,BF,l2);
      stan::math::hessian(f,param,fx,grad_fx,H_fx);
    }
  }

  return (H_fx);
}




//' evd_dnorm_hess_stan
//'
//' This is an attempt to use stan's AD features to calculate a hessian
//' for the RSSp likelihood
//'
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd fgem_hess_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double l2=0.0,const double l1 = 0.0,const bool log_BF=false){
    using Mt = Eigen::MatrixXd;

  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_fx;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }
  if(l2<0 or l1<0){
    Rcpp::stop("penalties (l2/L2 and l1) must be non-negative");
  }

  double fx=0;
  if(l1!=0){
    if(!log_BF){
      fgem_lik_l1<Mt> f(X,BF,l2,l1);
      stan::math::hessian(f,param,fx,grad_fx,H_fx);
    }else{
      log_fgem_lik_l1<Mt> f(X,BF,l2,l1);
      stan::math::hessian(f,param,fx,grad_fx,H_fx);
    }
  }else{
    if(!log_BF){
      fgem_lik<Mt> f(X,BF,l2);
      stan::math::hessian(f,param,fx,grad_fx,H_fx);
    }else{
      log_fgem_lik<Mt> f(X,BF,l2);
      stan::math::hessian(f,param,fx,grad_fx,H_fx);
    }
  }

  return (H_fx);
}




//'
//' This is an attempt to use stan's AD features to calculate a gradient
//' for the FGEM likelihood
//'
//' @export
//[[Rcpp::export]]
double sp_fgem_lik_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::SparseMatrix<double>> X, const Eigen::Map<Eigen::ArrayXd> BF,const double l2=0.0,const double l1=0.0,const bool log_BF=false){
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  using Mt = Eigen::SparseMatrix<double>;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }
  if(l2<0 or l1<0){
    Rcpp::stop("penalties (l2/L2 and l1) must be non-negative");
  }
  if(l1!=0){
    if(log_BF){
      return log_fgem_lik_l1<Mt>(X,BF,l2,l1)(param);
    }
    return fgem_lik_l1<Mt>(X,BF,l2,l1)(param);
  }
  if(log_BF){
    return log_fgem_lik<Mt>(X,BF,l2)(param);
  }
  return fgem_lik<Mt>(X,BF,l2)(param);
}

//[[Rcpp::export]]
double l1_norm_penalty(Eigen::Map<Eigen::ArrayXd> par, const double l1){
  return(l1_penalty<double>(par,l1));
}

//[[Rcpp::export]]
double l2_norm_penalty(Eigen::Map<Eigen::ArrayXd> par, const double l2){
  return(l2_penalty<double>(par,l2));
}

//[[Rcpp::export]]
double fgem_lik_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double l2=0.0,const double l1=0.0,const bool log_BF=false){
    using Mt = Eigen::MatrixXd;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);

  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }
  if(l2<0 or l1<0){
    Rcpp::stop("penalties (l2/L2 and l1) must be non-negative");
  }
  if(l1!=0){
    if(log_BF){
      return log_fgem_lik_l1<Mt>(X,BF,l2,l1)(param);
    }
    return fgem_lik_l1<Mt>(X,BF,l2,l1)(param);
  }
  if(log_BF){
    return log_fgem_lik<Mt>(X,BF,l2)(param);
  }
  return fgem_lik<Mt>(X,BF,l2)(param);
}



Rcpp::NumericVector sp_grad_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto l2 = Rcpp::as<double>(e["l2"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(X);

  auto log_sexp = Rcpp::as<SEXP>(e["log_BF"]);
  bool logs = Rcpp::as<bool>(log_sexp);

  return Rcpp::wrap(sp_fgem_grad_stan(par,sX,BF,l2,0.0,logs));
}


Rcpp::NumericVector grad_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto l2 = Rcpp::as<double>(e["l2"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  auto log_sexp = Rcpp::as<SEXP>(e["log_BF"]);
  bool logs = Rcpp::as<bool>(log_sexp);
  return Rcpp::wrap(fgem_grad_stan(par,sX,BF,l2,0.0,logs));
}



Rcpp::NumericVector lik_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto l2 = Rcpp::as<double>(e["l2"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd >>(X);

  auto log_sexp = Rcpp::as<SEXP>(e["log_BF"]);
  bool logs = Rcpp::as<bool>(log_sexp);
  return Rcpp::wrap(fgem_lik_stan(par,sX,BF,l2,0.0,logs));
}


Rcpp::NumericVector sp_lik_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto l2 = Rcpp::as<double>(e["l2"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(X);

  auto log_sexp = Rcpp::as<SEXP>(e["log_BF"]);
  bool logs = Rcpp::as<bool>(log_sexp);

  return Rcpp::wrap(sp_fgem_lik_stan(par,sX,BF,l2,0.0,logs));
}


//[[Rcpp::export]]
Rcpp::List make_env_obj(bool sparse=false,const bool log=false){
  using namespace Rcpp;
  if(!log){
    if(!sparse){
      return List::create(_["lik"]=XPtr<efuncPtr>(new efuncPtr(&lik_ext)),
                          _["grad"]=XPtr<efuncPtr>(new efuncPtr(&grad_ext)));
    }else{
      return List::create(_["lik"]=XPtr<efuncPtr>(new efuncPtr(&sp_lik_ext)),
                          _["grad"]=XPtr<efuncPtr>(new efuncPtr(&sp_grad_ext)));
    }
  }else{
    if(!sparse){
      return List::create(_["lik"]=XPtr<efuncPtr>(new efuncPtr(&lik_ext)),
                          _["grad"]=XPtr<efuncPtr>(new efuncPtr(&grad_ext)));
    }else{
      return List::create(_["lik"]=XPtr<efuncPtr>(new efuncPtr(&sp_lik_ext)),
                          _["grad"]=XPtr<efuncPtr>(new efuncPtr(&sp_grad_ext)));
    }
  }
}




//[[Rcpp::export]]
Rcpp::NumericVector log_1p_exp(SEXP x){

  SEXP rx = Rcpp::clone(x);
  double *xb= REAL(rx);
  const size_t p = LENGTH(rx);
  double *xe= xb+p;
  std::transform(xb,xe,xb,Rf_log1pexp);
  return rx;
  // return Rcpp::NumericVector::import_transform(x.begin(),x.end(),[](const double d){
  //                                                                  return(Rf_log1pexp(d));
  //                                                                });
}
