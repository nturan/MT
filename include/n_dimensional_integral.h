#pragma once

#include "ndimvegas.h"
#include <vector>
#include <array>
#include <map>
#include <complex>
#include <functional>

typedef std::function<double(std::map<std::string,double>,double&)> Integrand;
typedef std::map<std::string, double> IntegrationVariablesMap;
typedef std::map<std::string, std::pair<double, double>> IntegrationLimitsMap;
typedef std::map<std::string, 
    std::tuple<double, double, double>> IntegrationResultsMap;

IntegrationVariablesMap MapToHyperCube(
    IntegrationLimitsMap limits, double x[], double& jac,
    IntegrationVariablesMap constants = IntegrationVariablesMap{} );

class Integral{
public:
  Integral( IntegrationLimitsMap limits, 
            std::map<std::string, Integrand> integrands,
      IntegrationVariablesMap constants = IntegrationVariablesMap{});
  Integral(std::vector<Integral*> integrals, 
      std::pair<double, double> ratio = { 0.5, 0.5 });
  ~Integral();
  std::tuple<double, double, double> ExecuteVegas( 
      int entry, int itmx, int ncall, int nprn=0 );
    
  IntegrationResultsMap results_, results1_, results2_;
  std::map<std::string, Integrand> integrands_, integrands1_, integrands2_;
  IntegrationLimitsMap limits_, limits1_, limits2_;
  IntegrationVariablesMap constants_, constants1_, constants2_;

protected:
  int dimension_, dimension1_, dimension2_, 
      calls_ = 100000, iterations_ = 0, length_;
  std::pair<double, double> ratio_;
  bool performing_combined_integration_ = false;
  double f ( double x[], double wgt, double res[] );
  double CombinedIntegrationF(double x[5], double wgt, double res[]);
  int getDimension() { return dimension_; }
  int VegasIterations() { return iterations_; }  
  int VegasPoints() { return calls_; }
  double VegasAcc() { return 1.0E-10; }  
  friend class Vegas<Integral>;
  Vegas<Integral> *integrator_;  
};