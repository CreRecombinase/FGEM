#include "fgem.h"
#include <progress.hpp>
#include "progress_bar.hpp"
//[[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(unwindProtect)]]
//[[Rcpp::depends(BH)]]

typedef Rcpp::NumericVector (*efuncPtr)(SEXP,SEXP);


//' Gradient for FGEM likelihood
//' This is an attempt to use stan's AD features to optimize the FGEM likelihood
//' @export
//[[Rcpp::export]]
Eigen::ArrayXd fgem_grad_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false,const bool log_BF=false){
  using Mt = Eigen::MatrixXd;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  double fx=0;

  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }



  if(!log_BF){
    if(!neg){
      fgem_lik<Mt,1> f(X,BF,prec);
      stan::math::gradient(f,param,fx,grad_fx);
    }else{
      fgem_lik<Mt,-1> f(X,BF,prec);
      stan::math::gradient(f,param,fx,grad_fx);
    }
  }else{
    if(!neg){
      log_fgem_lik<Mt,1> f(X,BF,prec);
      stan::math::gradient(f,param,fx,grad_fx);
    }else{
      log_fgem_lik<Mt,-1> f(X,BF,prec);
      stan::math::gradient(f,param,fx,grad_fx);
    }
  }
  return(grad_fx.array());
}


//' Gradient for FGEM likelihood
//' This is an attempt to use stan's AD features to optimize the FGEM likelihood
//' @export
//[[Rcpp::export]]
Eigen::ArrayXd sp_fgem_grad_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::SparseMatrix<double>> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false,const bool log_BF=false){
  using Mt = Eigen::SparseMatrix<double>;
Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;

  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }


  double fx=0;
  if(!log_BF){
    if(!neg){
      fgem_lik<Mt,1> f(X,BF,prec);
      stan::math::gradient(f,param,fx,grad_fx);
    }else{
      fgem_lik<Mt,-1> f(X,BF,prec);
      stan::math::gradient(f,param,fx,grad_fx);
    }
  }else{
    if(!neg){
      log_fgem_lik<Mt,1> f(X,BF,prec);
      stan::math::gradient(f,param,fx,grad_fx);
    }else{
      log_fgem_lik<Mt,-1> f(X,BF,prec);
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
Eigen::MatrixXd sp_fgem_hess_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::SparseMatrix<double>> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false,const bool log_BF=false){

  using Mt = Eigen::SparseMatrix<double>;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_fx;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }

  double fx=0;
  if(!log_BF){
    if(!neg){
      fgem_lik<Mt,1> f(X,BF,prec);
      stan::math::hessian(f, param, fx, grad_fx, H_fx);
    }else{
      fgem_lik<Mt,-1> f(X,BF,prec);
      stan::math::hessian(f, param, fx, grad_fx, H_fx);
    }
  }else{
    if(!neg){
      log_fgem_lik<Mt,1> f(X,BF,prec);
      stan::math::hessian(f, param, fx, grad_fx, H_fx);
    }else{
      log_fgem_lik<Mt,-1> f(X,BF,prec);
      stan::math::hessian(f, param, fx, grad_fx, H_fx);
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
Eigen::MatrixXd fgem_hess_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false,const bool log_BF=false){
    using Mt = Eigen::MatrixXd;

  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_fx;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }

  double fx=0;
  if(!log_BF){
    if(!neg){
      fgem_lik<Mt,1> f(X,BF,prec);
      stan::math::hessian(f, param, fx, grad_fx, H_fx);
    }else{
      fgem_lik<Mt,-1> f(X,BF,prec);
      stan::math::hessian(f, param, fx, grad_fx, H_fx);
    }
  }else{
    if(!neg){
      log_fgem_lik<Mt,1> f(X,BF,prec);
      stan::math::hessian(f, param, fx, grad_fx, H_fx);
    }else{
      log_fgem_lik<Mt,-1> f(X,BF,prec);
      stan::math::hessian(f, param, fx, grad_fx, H_fx);
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
double sp_fgem_lik_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::SparseMatrix<double>> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false,const bool log_BF=false){
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  using Mt = Eigen::SparseMatrix<double>;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }

  if(log_BF){
    if(!neg)
      return log_fgem_lik<Mt,1>(X,BF,prec)(param);
    return log_fgem_lik<Mt,-1>(X,BF,prec)(param);
  }
  if(!neg)
    return fgem_lik<Mt,1>(X,BF,prec)(param);
  return fgem_lik<Mt,-1>(X,BF,prec)(param);
}

//[[Rcpp::export]]
double fgem_lik_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false,const bool log_BF=false){
    using Mt = Eigen::MatrixXd;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);

  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(X.rows()!=BF.size()){
    Rcpp::stop("BF must be of length: NROW(x) ("+std::to_string(X.rows())+"), and not: "+std::to_string(BF.size()));
  }
  if(log_BF){
    if(!neg)
      return log_fgem_lik<Mt,1>(X,BF,prec)(param);
    return log_fgem_lik<Mt,-1>(X,BF,prec)(param);
  }
  if(!neg)
    return fgem_lik<Mt,1>(X,BF,prec)(param);
  return fgem_lik<Mt,-1>(X,BF,prec)(param);
}

 
//' fit every column of X with the BFGS algorithm
//'
//' fit fgem using c++ BFGS algorithm
//'
//' @export
//[[Rcpp::export]]
Rcpp::List marginal_fgem_fit_bfgs(Rcpp::NumericMatrix X, const Eigen::Map<Eigen::ArrayXd> BF,double prec=0.0,const double epsilon=1e-6, const int max_iter=100,const bool progress=true,const bool log_BF=false){
  LBFGSpp::LBFGSBParam<double> param;
  param.epsilon = epsilon;
  param.max_iterations = max_iter;
  // param.past=1;
  // param.delta=1e-10;
  //  LBFGSpp::LineSearchMoreThuente<double>

  //  LBFGSpp::LBFGSSolver<double,LBFGSpp::LineSearchMoreThuente> solver(param);
  LBFGSpp::LBFGSBSolver<double> solver(param);

  const size_t g= BF.size();
  using namespace Rcpp;
  StringVector feature_name = colnames(X);
  auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  if(sX.rows()!=g){
      Rcpp::stop("rows of X must be equal to length of BF");
  }
  Eigen::ArrayXd sum_feat =  sX.colwise().sum();
  double* cr=sX.data();
  Eigen::Map<Eigen::MatrixXd> sc(cr,g,1);
  const size_t p=sX.cols();
  Progress pp(p, progress);                  
  auto iv = seq_len(p);
  Rcpp::NumericMatrix coeffMat(2,p);

  
  colnames(coeffMat) = colnames(X);
  Rcpp::NumericVector lik_vec(p);
  Rcpp::StringVector message(p);
  
  Rcpp::IntegerVector iter_vec(p);
  auto coefM = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(coeffMat);
  const char* msg = "converged";
  auto msg_exp = Rf_mkChar(msg);
  
  Eigen::VectorXd lb = Eigen::VectorXd::Constant(2, -150.0);
  Eigen::VectorXd ub = Eigen::VectorXd::Constant(2, 150.0);
  Eigen::VectorXd ptx = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd grd = Eigen::VectorXd::Zero(2);
  for(auto i:iv){
    if ( Progress::check_abort() )
      break;
    pp.increment();
    double* crp=cr+(i-1)*g;
    Eigen::Map<Eigen::MatrixXd> sc(crp,g,1);
    auto rmsg_exp=msg_exp;
    //    fun.update_X(Eigen::Map<Eigen::MatrixXd>(&cr,g,1));
    int niter=0;
    double fx=0;
    try{
      if(log_BF){
        auto fun = log_fgem_bfg<Eigen::MatrixXd>(log_fgem_lik<Eigen::MatrixXd,-1>(sc,BF,prec));
        niter = solver.minimize(fun,ptx,fx,lb,ub);
      }     else{
        auto fun = fgem_bfg<Eigen::MatrixXd>(fgem_lik<Eigen::MatrixXd,-1>(sc,BF,prec));
        niter = solver.minimize(fun,ptx,fx,lb,ub);

      }
    } catch(const std::runtime_error& error){
      niter =NA_INTEGER;
      fx = NA_REAL;
      rmsg_exp=Rf_mkChar(error.what());
    }

    iter_vec[i-1]=niter;
    lik_vec[i-1]=-fx;
    coefM.col(i-1)=ptx;
    message[i-1]=rmsg_exp;
  }
  
  return List::create(_["iter"]=iter_vec,
                      _["lik"]=lik_vec,
                      _["coeff"]=coeffMat,
                      _["message"]=message,
                      _["feat_sum"]=Rcpp::wrap(sum_feat),
                      _["feature_name"]=feature_name);
      

}


Rcpp::NumericVector sp_grad_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto prec = Rcpp::as<double>(e["prec"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(X);

  auto log_sexp = Rcpp::as<SEXP>(e["log_BF"]);
  bool logs = Rcpp::as<bool>(log_sexp);

  return Rcpp::wrap(sp_fgem_grad_stan(par,sX,BF,prec,true,logs));
}


Rcpp::NumericVector grad_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto prec = Rcpp::as<double>(e["prec"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  auto log_sexp = Rcpp::as<SEXP>(e["log_BF"]);
  bool logs = Rcpp::as<bool>(log_sexp);
  return Rcpp::wrap(fgem_grad_stan(par,sX,BF,prec,true,logs));
}



Rcpp::NumericVector lik_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto prec = Rcpp::as<double>(e["prec"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd >>(X);

  auto log_sexp = Rcpp::as<SEXP>(e["log_BF"]);
  bool logs = Rcpp::as<bool>(log_sexp);
  return Rcpp::wrap(fgem_lik_stan(par,sX,BF,prec,true,logs));
}


Rcpp::NumericVector sp_lik_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto prec = Rcpp::as<double>(e["prec"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(X);

  auto log_sexp = Rcpp::as<SEXP>(e["log_BF"]);
  bool logs = Rcpp::as<bool>(log_sexp);

  return Rcpp::wrap(sp_fgem_lik_stan(par,sX,BF,prec,true,logs));
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
