#pragma once

#include "vegas/ndimvegas.h"
#include <vector>
#include <array>
#include <map>
#include <complex>
#include <functional>

typedef std::function<double(std::map<std::string,double>,double&)> Integrand;

class Integral{
public:
  Integral( std::map<std::string, std::pair<double, double>>, std::map<std::string, Integrand> integrands );
  ~Integral();
  void ExecuteVegas( int entry, int itmx, int ncall, int nprn=0 );
    
  std::map<std::string, std::tuple<double, double, double>> results_;
  std::map<std::string, Integrand> integrands_;

protected:
  std::map<std::string, std::pair<double, double>> limits_;
  Integrand integrand_;
  int dimension_, calls_ = 100000, iterations_ = 0, length_;
  double f ( double x[], double wgt, double res[] );
  int getDimension() { return dimension_; } 
  int VegasIterations() { return iterations_; }  
  int VegasPoints() { return calls_; }
  double VegasAcc() { return 1.0E-10; }  
  friend class Vegas<Integral>;
  Vegas<Integral> *integrator_;  
};