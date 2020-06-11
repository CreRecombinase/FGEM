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
Eigen::ArrayXd fgem_grad_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false){
  using Mt = Eigen::MatrixXd;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  double fx=0;
  if(!neg){
    fgem_lik<Mt,1> f(X,BF,prec);
    stan::math::gradient(f,param,fx,grad_fx);
  }else{
    fgem_lik<Mt,-1> f(X,BF,prec);
    stan::math::gradient(f,param,fx,grad_fx);
  }
  return(grad_fx.array());
}

// 
// 
// SEXP fgem_deviance(SEXP y, SEXP mu, SEXP wt)
// {
//     int i, n = LENGTH(y), lmu = LENGTH(mu), lwt = LENGTH(wt), nprot = 1;
//     SEXP ans;
//     double mui, yi, *rmu, *ry, *rwt, *rans;
// 
//     if (!isReal(y)) {y = PROTECT(coerceVector(y, REALSXP)); nprot++;}
//     ry = REAL(y);
//     ans = PROTECT(shallow_duplicate(y));
//     rans = REAL(ans);
//     if (!isReal(mu)) {mu = PROTECT(coerceVector(mu, REALSXP)); nprot++;}
//     if (!isReal(wt)) {wt = PROTECT(coerceVector(wt, REALSXP)); nprot++;}
//     rmu = REAL(mu);
//     rwt = REAL(wt);
//     if (lmu != n && lmu != 1)
// 	error(_("argument %s must be a numeric vector of length 1 or length %d"),
// 	      "mu", n);
//     if (lwt != n && lwt != 1)
// 	error(_("argument %s must be a numeric vector of length 1 or length %d"),
// 	      "wt", n);
//     /* Written separately to avoid an optimization bug on Solaris cc */
//     if(lmu > 1) {
// 	for (i = 0; i < n; i++) {
// 	    mui = rmu[i];
// 	    yi = ry[i];
// 	    rans[i] =  rwt[lwt > 1 ? i : 0] *
//               mui*yi + (1-mui);
// 	}
//     } else {
// 	mui = rmu[0];
// 	for (i = 0; i < n; i++) {
// 	    yi = ry[i];
//             rans[i] =  log(rwt[lwt > 1 ? i : 0] *
//                            (mui*yi + (1-mui)));
// 	}
//     }
// 
//     UNPROTECT(nprot);
//     return ans;
// }
//   
    


//' Gradient for FGEM likelihood
//' This is an attempt to use stan's AD features to optimize the FGEM likelihood
//' @export
//[[Rcpp::export]]
Eigen::ArrayXd sp_fgem_grad_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::SparseMatrix<double>> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false){
  using Mt = Eigen::SparseMatrix<double>;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  double fx=0;
  if(!neg){
    fgem_lik<Mt,1> f(X,BF,prec);
    stan::math::gradient(f,param,fx,grad_fx);
  }else{
    fgem_lik<Mt,-1> f(X,BF,prec);
    stan::math::gradient(f,param,fx,grad_fx);
  }
  return(grad_fx.array());
}



//'
//' This is an attempt to use stan's AD features to calculate a gradient
//' for the FGEM likelihood
//'
//' @export
//[[Rcpp::export]]
double sp_fgem_lik_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::SparseMatrix<double>> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false){
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  using Mt = Eigen::SparseMatrix<double>;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(!neg){
    return fgem_lik<Mt,1>(X,BF,prec)(param);
  }else{
    return fgem_lik<Mt,-1>(X,BF,prec)(param);
  }
}

//[[Rcpp::export]]
double fgem_lik_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false,const bool log_BF=false){
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  using Mt = Eigen::MatrixXd;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(log_BF){
    if(!neg)
      return log_fgem_lik<Mt,1>(X,BF,prec)(param);
    return log_fgem_lik<Mt,-1>(X,BF,prec)(param);
  }
  if(!neg)
    return fgem_lik<Mt,1>(X,BF,prec)(param);
  return fgem_lik<Mt,-1>(X,BF,prec)(param);
}



//' evd_dnorm_hess_stan
//' 
//' This is an attempt to use stan's AD features to calculate a hessian 
//' for the RSSp likelihood
//'
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd fgem_hess_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0){

  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  fgem_lik<Eigen::MatrixXd,1> f(X,BF,prec);
  double fx=0;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_fx;
  stan::math::hessian(f, param, fx, grad_fx, H_fx);
  return (H_fx);
}




//' fit fgem with bfgs
//'
//' fit fgem using c++ BFGS algorithm
//'
//' @export
//[[Rcpp::export]]
Rcpp::List fgem_fit_bfgs(const Eigen::ArrayXd par,SEXP X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const double epsilon=1e-6, const int max_iter=100){
  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = epsilon;
  param.max_iterations = max_iter;
  //  LBFGSpp::LineSearchMoreThuente<double>
  LBFGSpp::LBFGSSolver<double> solver(param);
  Eigen::VectorXd x = par;
  double fx;
  int niter=0;
  if(Rf_inherits(X,"dgCMatrix")){
  // Create solver and function object
    auto sX = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(X);
    if(sX.cols()!=par.size()-1)
      Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(sX.cols()+1)+"), and not: "+std::to_string(par.size()));
    fgem_bfg<Eigen::SparseMatrix<double>> fun(fgem_lik<Eigen::SparseMatrix<double>,-1>(sX,BF,prec));
    niter = solver.minimize(fun,x,fx);
  }else{
    auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
    if(sX.cols()!=par.size()-1)
      Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(sX.cols()+1)+"), and not: "+std::to_string(par.size()));
    fgem_bfg<Eigen::MatrixXd> fun(fgem_lik<Eigen::MatrixXd,-1>(sX,BF,prec));
    niter = solver.minimize(fun,x,fx);
  }
  // x will be overwritten to be the best point found
  using namespace Rcpp;
  return List::create(_["iter"]=Rcpp::wrap(niter),
                      _["par"]=Rcpp::wrap(x),
                      _["obj"]=Rcpp::wrap(fx));
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
  return Rcpp::wrap(sp_fgem_grad_stan(par,sX,BF,prec,true));
}


Rcpp::NumericVector grad_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto prec = Rcpp::as<double>(e["prec"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  return Rcpp::wrap(fgem_grad_stan(par,sX,BF,prec,true));
}



Rcpp::NumericVector lik_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto prec = Rcpp::as<double>(e["prec"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd >>(X);
  return Rcpp::wrap(fgem_lik_stan(par,sX,BF,prec,true));
}


Rcpp::NumericVector sp_lik_ext(SEXP xs, SEXP env){
  auto par =  Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  auto BF = Rcpp::as<Eigen::Map<Eigen::ArrayXd >>(e["BF"]);
  auto prec = Rcpp::as<double>(e["prec"]);
  auto X =Rcpp::as<SEXP>(e["X"]);
  auto sX = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(X);
  return Rcpp::wrap(sp_fgem_lik_stan(par,sX,BF,prec,true));
}


//[[Rcpp::export]]
Rcpp::List make_env_obj(bool sparse=false){
  using namespace Rcpp;
  if(!sparse){
    return List::create(_["lik"]=XPtr<efuncPtr>(new efuncPtr(&lik_ext)),
                        _["grad"]=XPtr<efuncPtr>(new efuncPtr(&grad_ext)));
  }else{
    return List::create(_["lik"]=XPtr<efuncPtr>(new efuncPtr(&sp_lik_ext)),
                        _["grad"]=XPtr<efuncPtr>(new efuncPtr(&sp_grad_ext)));
  }
}
