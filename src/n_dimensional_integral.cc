#include "n_dimensional_integral.h"

Integral::Integral(std::vector<Integral*> integrals) {
    if (integrals.size() != 2) {
        std::cerr << "combined integration with " << integrals.size()
            << " has not been implemented yet..." << std::endl;
        exit(1);
    }
    performing_combined_integration_ = true;
    limits1_ = integrals.at(0)->limits_;
    limits2_ = integrals.at(1)->limits_;
    integrands1_ = integrals.at(0)->integrands_;
    integrands2_ = integrals.at(1)->integrands_;
    dimension1_ = integrals.at(0)->dimension_;
    dimension2_ = integrals.at(1)->dimension_;
    dimension_ = dimension1_ + dimension2_ + 1;
    constants1_ = integrals.at(0)->constants_;
    constants2_ = integrals.at(1)->constants_;
    length_ = 1;
    integrator_ = new Vegas<Integral>(*this, length_);
}

Integral::Integral( 
    IntegrationLimitsMap limits, std::map<std::string, Integrand> integrands, IntegrationVariablesMap constants ){
    limits_ = limits;
    constants_ = constants;
    dimension_ = limits.size() - constants.size();
    integrands_ = integrands;
    length_ = integrands_.size();    
    integrator_ = new Vegas<Integral>(*this, length_);
}

Integral::~Integral()
{
    delete integrator_;
}

void Integral::ExecuteVegas(
    int entry, int itmx, int ncall,int nprn) {

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
  
  if (performing_combined_integration_) {
    results_["default"] = std::make_tuple( integrator_->avgi_[0],
                                            integrator_->sd_[0],
                                            integrator_->chi2a_[0] );
  }


   int i = 0;
   for ( auto it = integrands_.begin(); it != integrands_.end(); ++it ) {
     std::string name = it->first;
     results_[name] = std::make_tuple(integrator_->avgi_[i],
                                      integrator_->sd_[i],
                                      integrator_->chi2a_[i]);
     ++i;
  }                    
}


double Integral::f(double x[],double wgt, double res[]) {

    if (performing_combined_integration_) {
        return CombinedIntegrationF(x, wgt, res);
    }

  double jac = 1.0;
  

  IntegrationVariablesMap variables = MapToHyperCube(limits_, x, jac, constants_);

  int i = 0;
  int number_of_calls = integrator_->getcalls();
  double weight = wgt / number_of_calls / iterations_;
  for ( auto it = integrands_.begin(); it != integrands_.end(); ++it ) {
      res[i] = it->second(variables, weight) * jac;
      ++i;
  }
  
  // std::cout << result.at(0) << std::endl;

  return(res[0]);
}

double Integral::CombinedIntegrationF(double x[5], double wgt, double res[]) {
    double jac1 = 1.0, jac2 = 1.0;
    std::vector<double> x1, x2;
    for (int i = 0; i < dimension1_; ++i) {
        x1.push_back(x[i]);
    }
    for (int i = dimension1_; i < dimension1_ + dimension2_; ++i) {
        x2.push_back(x[i]);
    }

    IntegrationVariablesMap variables1 = MapToHyperCube(limits1_, &x1[0], jac1, constants1_);
    IntegrationVariablesMap variables2 = MapToHyperCube(limits2_, &x2[0], jac2, constants2_);

    //std::cout << "variables1: ";
    //for (auto& l : variables1) {
    //    std::cout << l.first << "=" << l.second << ", ";
    //}
    //std::cout << std::endl;
    //std::cout << "variables2: ";
    //for (auto& l : variables2) {
    //    std::cout << l.first << "="<< l.second << ", ";
    //}
    //std::cout << std::endl;
    bool performing_first_integration = x[dimension_ - 1] < 0.1;

    int number_of_calls = integrator_->getcalls();
    double weight = wgt / number_of_calls / iterations_;
    if (performing_first_integration) {
        res[0] = integrands1_["default"](variables1, weight) * jac1 * 10.0;
    }
    else {
        res[0] = integrands2_["default"](variables2, weight) * jac2 * 10.0 / 9.0;
    }
    return res[0];
}

IntegrationVariablesMap MapToHyperCube(
    IntegrationLimitsMap limits, double x[], double& jacobian, IntegrationVariablesMap constants)
{
    IntegrationVariablesMap variables;
    double jac = jacobian;
    int index = 0;
    for (auto it = limits.begin(); it != limits.end(); ++it) {
        std::string name = it->first;
        std::pair<double, double> limit = it->second;
        double lower_limit = limit.first;
        double upper_limit = limit.second;
        jac *= (upper_limit - lower_limit);
        if (constants.count(name) > 0) {
            variables[name] = constants[name];
            continue;
        }
        variables[name] = x[index] * (upper_limit - lower_limit) + lower_limit;
        ++index;
    }
    jacobian = jac;
    return variables;
}
