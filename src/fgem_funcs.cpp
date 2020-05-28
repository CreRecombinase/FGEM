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

double fgem_lik_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const bool neg=false){
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  using Mt = Eigen::MatrixXd;
  if(X.cols()!=par.size()-1)
    Rcpp::stop("par must be of length: NCOL(x)+1 ("+std::to_string(X.cols()+1)+"), and not: "+std::to_string(par.size()));
  if(!neg){
    return fgem_lik<Mt,1>(X,BF,prec)(param);
  }else{
    return fgem_lik<Mt,-1>(X,BF,prec)(param);
  }
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
Rcpp::List marginal_fgem_fit_bfgs(Rcpp::NumericMatrix X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0,const double epsilon=1e-6, const int max_iter=100,const bool progress=true){
  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = epsilon;
  param.max_iterations = max_iter; 
  //  LBFGSpp::LineSearchMoreThuente<double>

  LBFGSpp::LBFGSSolver<double> solver(param);
  double fx;
  const size_t g= BF.size();
  using namespace Rcpp;
  StringVector feature_name = colnames(X);
  auto sX = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  if(sX.rows()!=g){
      Rcpp::stop("rows of X must be equal to length of BF");
  }
  Eigen::ArrayXd sum_feat =  sX.colwise().sum();
  double& cr=sX.coeffRef(0,0);
  Eigen::Map<Eigen::MatrixXd> sc(&cr,g,1);
  const size_t p=sX.cols();
  Progress pp(p, progress);                  
  auto iv = seq_len(p);
  Rcpp::NumericMatrix coeffMat(2,p);
  colnames(coeffMat) = colnames(X);
  Rcpp::NumericVector lik_vec(p);
  Rcpp::IntegerVector iter_vec(p);
  Rcpp::IntegerVector obj_vec(p);
  auto coefM = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(coeffMat);
      
      
  fgem_bfg<Eigen::MatrixXd> fun(fgem_lik<Eigen::MatrixXd,-1>(sc,BF,prec));
  Eigen::VectorXd ptx = Eigen::VectorXd::Zero(2);
  for(auto i:iv){
    if ( Progress::check_abort() )
      break;
    pp.increment();
    cr=sX.coeffRef(0,i-1);
    fun.update_X(Eigen::Map<Eigen::MatrixXd>(&cr,g,1));
    int niter=0;
    niter = solver.minimize(fun,ptx,fx);
    iter_vec[i]=niter;
    obj_vec[i]=-fx;
    coefM.col(i)=ptx;        
  }
  
  return List::create(_["iter"]=iter_vec,
                      _["lik"]=obj_vec,
                      _["coeff"]=coeffMat,
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
