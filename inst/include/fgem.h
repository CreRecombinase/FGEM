#pragma once
#include <RcppEigen.h>
#include <stan/math.hpp>
#include <LBFGS.h>
#include <LBFGSB.h>
#include <stan/math/mix/mat.hpp>
#include <cmath>
#include <cstddef>
#include <limits>

  
template<typename T,typename U,int neg=-1>
inline T FGEM_log_lik(const  Eigen::Array<T,Eigen::Dynamic,1> &Beta, const Eigen::Map<U> X, const Eigen::Map<Eigen::ArrayXd> BF) noexcept{

  // arma::vec pvec = 1/(1+arma::exp(-(x*Beta)));
  // Rcpp::Rcout<<"pvec:"<<pvec<<std::endl;


  const size_t p= Beta.size();
  auto pvec =  stan::math::inv_logit((X*(Beta.tail(p-1).matrix())).array()+Beta[0]);
  return ((pvec.array()*BF)+(1.0-pvec.array())).log().sum()*neg;

}


template<typename T,typename U,int neg=-1>
inline T FGEM_log_lik_l2(const  Eigen::Array<T,Eigen::Dynamic,1> &Beta, const Eigen::Map<U> X, const Eigen::Map<Eigen::ArrayXd> BF,const double prec=0.0) noexcept{

  // arma::vec pvec = 1/(1+arma::exp(-(x*Beta)));
  // Rcpp::Rcout<<"pvec:"<<pvec<<std::endl;

//  arma::vec uvec = (pvec%B)/((pvec%B)+(1-pvec));
  const size_t p= Beta.size();
  auto pvec =  stan::math::inv_logit((X*(Beta.tail(p-1).matrix())).array()+Beta[0]);
  return (((pvec.array()*BF)+(1.0-pvec.array())).log().sum() - 0.5 * (Beta.array().tail(p-1).square()*prec).sum()) *neg;
}


template<typename U,int neg=-1>
struct fgem_lik {
  Eigen::Map<U> X;
  Eigen::Map<Eigen::ArrayXd> BF;
  double prec;
  fgem_lik(const Eigen::Map<U> X_,  const Eigen::Map<Eigen::ArrayXd> BF_,const double prec_=0.0): X(X_),BF(BF_),prec(prec_){}
  void update_X(Eigen::Map<U> &&x){
    X=x;
  }
  
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &par) const noexcept{
      Eigen::Array<T,Eigen::Dynamic,1> tcvec(par);
      return FGEM_log_lik_l2<T,U,neg>(tcvec,X,BF,prec);
  }
  
};


template<typename U>
class fgem_bfg
{
private:
  fgem_lik<U,-1> f;
  mutable Eigen::VectorXd tgrad;
public:
  fgem_bfg(fgem_lik<U,-1> && f_) : f(f_),tgrad() {}
  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
  {
    double fx=0;
    stan::math::gradient(f,x,fx,grad);
    return fx;
  }
  void update_X(Eigen::Map<U> &&x){
    f.update_X(std::forward<Eigen::Map<U>>(x));
  }
  Rcpp::NumericVector grad(SEXP xs) const noexcept{
    auto x =  Rcpp::as<Eigen::Map<Eigen::VectorXd>>(xs);
    tgrad.resize(x.size());
    double fx=0;
    stan::math::gradient(f,x,fx,tgrad);
    return Rcpp::wrap(tgrad);
  }
  Rcpp::NumericVector lik(SEXP xs) const noexcept{
    auto fr = f(Rcpp::as<Eigen::VectorXd>(xs));
    return Rcpp::NumericVector::create(fr);
  }
};
