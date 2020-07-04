#pragma once
#include <stan/math.hpp>
#include <RcppEigen.h>
#include <LBFGS.h>
#include <LBFGSB.h>
#include <stan/math/mix/mat.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>

#if __has_include(<execution>)
#define have_transform_reduce 1
#else
#define have_transform_reduce 0
#endif



template<typename T,typename U>
inline auto l1_penalty(const Eigen::Matrix<T,Eigen::Dynamic,1> &Beta,const U l1) {
  const size_t p=Beta.size();
  return l1*(Beta.segment(1,Beta.size()-1).template lpNorm<1>());
}

template<typename T,typename U>
inline auto l2_penalty(const Eigen::Matrix<T,Eigen::Dynamic,1> &Beta,const U l2) {
  return (l2/2.0)*(Beta.segment(1,Beta.size()-1).squaredNorm());
}

template<typename T, typename U>
inline auto FGEM_log_lik(const  Eigen::Array<T,Eigen::Dynamic,1> &Beta, const Eigen::Map<U> X, const Eigen::Map<Eigen::ArrayXd> BF) noexcept{
  const size_t p= Beta.size();
  auto pvec =  stan::math::inv_logit((X*(Beta.tail(p-1).matrix())).array()+Beta[0]);
  return (((pvec.array()*BF)+(1.0-pvec.array())).log().sum());
}


template<typename T,typename U>
inline T log_FGEM_log_lik(const  Eigen::Array<T,Eigen::Dynamic,1> &Beta, const Eigen::Map<U> X, const Eigen::Map<Eigen::ArrayXd> log_BF) noexcept{

  const size_t p= Beta.size();
  T init =  0.0;
  Eigen::Array<T,Eigen::Dynamic,1> xb = (X*(Beta.tail(p-1).matrix())).array()+Beta[0];

  //  return (stan::math::log1p_exp(+log_BF)-stan::math::log1p_exp(xb)).sum();

  #if have_transform_reduce
  return  std::inner_product(xb.data(),xb.data()+xb.size(),log_BF.data(),init,std::plus<>(),
                                    [](T a, double b) -> T{
                                      return stan::math::log1p_exp(a+b)-stan::math::log1p_exp(a);
                                    });
  #else
  return  std::inner_product(xb.data(),xb.data()+xb.size(),log_BF.data(),init,std::plus<>(),
                                    [](T a, double b) -> T{
                                      return stan::math::log1p_exp(a+b)-stan::math::log1p_exp(a);

                                    });
  #endif

}


template <typename T, typename U,auto TF>
inline T l2_penalized_lik(const  Eigen::Array<T,Eigen::Dynamic,1> &Beta, const Eigen::Map<U> X, const Eigen::Map<Eigen::ArrayXd> BF,const T l2) noexcept{
  return -TF(Beta,X,BF) + l2_penalty<T,T>(Beta,l2);
}

template <typename T, typename U,auto TF>
inline T l1_l2_penalized_lik(const  Eigen::Array<T,Eigen::Dynamic,1> &Beta, const Eigen::Map<U> X, const Eigen::Map<Eigen::ArrayXd> BF,const T l2,const T l1) noexcept{
  return -TF(Beta,X,BF) + l2_penalty<T,T>(Beta,l2)+l1_penalty<T,T>(Beta,l1);
}

// template <typename T, typename U,auto TF>
// inline T ap_l1_l2_penalized_lik(const  Eigen::Array<T,Eigen::Dynamic,1> &Beta, const Eigen::Map<U> X, const Eigen::Map<Eigen::ArrayXd> BF) noexcept{
//   return l1_l2_penalized_lik(Beta.tail(Beta.size()-2),X,BF,Beta[0],Beta[1]);
// }



template<typename U>
struct fgem_lik {
  Eigen::Map<U> X;
  Eigen::Map<Eigen::ArrayXd> BF;
  double l2;
  fgem_lik(const Eigen::Map<U> X_,  const Eigen::Map<Eigen::ArrayXd> BF_,const double l2_=0.0): X(X_),BF(BF_),l2(l2_){}
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &par) const noexcept{
      Eigen::Array<T,Eigen::Dynamic,1> tcvec(par);
      return l2_penalized_lik<T,U,FGEM_log_lik<T,U> >(tcvec,X,BF,l2);
  }
};

template<typename U>
struct log_fgem_lik {
  Eigen::Map<U> X;
  Eigen::Map<Eigen::ArrayXd> log_BF;
  double l2;
  log_fgem_lik(const Eigen::Map<U> X_,  const Eigen::Map<Eigen::ArrayXd> log_BF_,const double l2_=0.0): X(X_),log_BF(log_BF_),l2(l2_){}
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &par) const noexcept{
      Eigen::Array<T,Eigen::Dynamic,1> tcvec(par);
      return l2_penalized_lik<T,U,log_FGEM_log_lik<T,U> >(tcvec,X,log_BF,l2);
  }
};


template<typename U>
struct fgem_lik_l1 {
  Eigen::Map<U> X;
  Eigen::Map<Eigen::ArrayXd> BF;
  double l2;
  double l1;
  fgem_lik_l1(const Eigen::Map<U> X_,  const Eigen::Map<Eigen::ArrayXd> BF_,const double l2_,const double l1_): X(X_),BF(BF_),l2(l2_),l1(l1_){}
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &par) const noexcept{
      Eigen::Array<T,Eigen::Dynamic,1> tcvec(par);
      return l1_l2_penalized_lik<T,U,FGEM_log_lik<T,U> >(tcvec,X,BF,l2,l1);
  }

};

template<typename U>
struct log_fgem_lik_l1 {
  Eigen::Map<U> X;
  Eigen::Map<Eigen::ArrayXd> log_BF;
  double l2;
  double l1;
  log_fgem_lik_l1(const Eigen::Map<U> X_,  const Eigen::Map<Eigen::ArrayXd> log_BF_,const double l2_,const double l1_): X(X_),log_BF(log_BF_),l2(l2_),l1(l1_){}
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &par) const noexcept{
      Eigen::Array<T,Eigen::Dynamic,1> tcvec(par);
      return l1_l2_penalized_lik<T,U,log_FGEM_log_lik<T,U> >(tcvec,X,log_BF,l2,l1);
  }
};



// template<typename U>
// struct log_fgem_lik_l1g {
//   Eigen::Map<U> X;
//   Eigen::Map<Eigen::ArrayXd> log_BF;
//   log_fgem_lik_l1g(const Eigen::Map<U> X_,  const Eigen::Map<Eigen::ArrayXd> log_BF_): X(X_),log_BF(log_BF_){}
//   template <typename T>
//   T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &par) const noexcept{
//       Eigen::Array<T,Eigen::Dynamic,1> tcvec(par);
//       return ap_l1_l2_penalized_lik<T,U,log_FGEM_log_lik<T,U> >(tcvec,X,log_BF);
//   }
// };



// template<typename U,template<typename> typename L = fgem_lik>
// class fgem_bfg
// {
// private:
//   L<U> f;
//   mutable Eigen::VectorXd tgrad;
// public:
//   fgem_bfg(L<U> && f_) : f(f_),tgrad() {}
//   double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
//   {
//     double fx=0;
//     stan::math::gradient(f,x,fx,grad);
//     return fx;
//   }
//   Rcpp::NumericVector grad(SEXP xs) const noexcept{
//     auto x =  Rcpp::as<Eigen::Map<Eigen::VectorXd>>(xs);
//     tgrad.resize(x.size());
//     double fx=0;
//     stan::math::gradient(f,x,fx,tgrad);
//     return Rcpp::wrap(tgrad);
//   }
//   Rcpp::NumericVector lik(SEXP xs) const noexcept{
//     auto fr = f(Rcpp::as<Eigen::VectorXd>(xs));
//     return Rcpp::NumericVector::create(fr);
//   }
// };
