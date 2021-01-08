#include "n_dimensional_integral.h"

Integral::Integral( std::map<std::string, std::pair<double, double>> limits, std::map<std::string, Integrand> integrands ){
    limits_ = limits;
    dimension_ = limits.size();
    integrands_ = integrands;
    length_ = integrands_.size();    
    integrator_ = new Vegas<Integral>(*this, length_);
}

Integral::~Integral()
{
    delete integrator_;
    std::cout << "integrator is deleted" << std::endl;
}

double Integral::ExecuteVegas(int entry, int itmx, int ncall,int nprn) {

  iterations_ += itmx;
  calls_=ncall;

  integrator_->setitmx(iterations_);
  integrator_->setncall(calls_);
 

  integrator_->setnprn(nprn);

  switch(entry){
   case 0: integrator_->vegas();
           break;
   case 1: integrator_->vegas1();
           break;
   case 2: integrator_->vegas2();
           break;
   default: std::cout << "wrong vegas entry level! must be 0,1 or 2" << std::endl;
            exit(0);
            break;
   }       
   
   int i = 0;
   for ( auto it = integrands_.begin(); it != integrands_.end(); ++it ) {
     std::string name = it->first;
    integral_[name] = integrator_->avgi_[i];
    error_[name] = integrator_->sd_[i];
    chi2a_[name] = integrator_->chi2a_[i];
    ++i;
  }                    

  
  return integrator_->avgi;
}


double Integral::f(double x[],double wgt, double res[]) {
  /*
   *  Differential integrand
   */

  std::map<std::string, double> variables;
  double jac = 1.0;
  int index = 0;
  for ( std::map<std::string, std::pair<double, double>>::iterator it = limits_.begin(); it != limits_.end(); ++it ){
    std::string name = it->first;
    std::pair<double, double> limit = it->second;
    double lower_limit = limit.first;
    double upper_limit = limit.second;
    variables[name] = x[index]*(upper_limit - lower_limit) + lower_limit;
    jac *= (upper_limit - lower_limit);
    ++index;
  }

  int i = 0;
  int number_of_calls = integrator_->getcalls();
  double weight = wgt / number_of_calls / iterations_;
  for ( auto it = integrands_.begin(); it != integrands_.end(); ++it ) {
      //std::cout << weight << std::endl;
      res[i] = it->second(variables, weight) * jac;
      ++i;
  }
  
  // std::cout << result.at(0) << std::endl;

  return(res[0]);
}