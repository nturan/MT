#include "matrix_element_method.h"

void EvaluateLoEvents(EventGenerator* eg, Parameters* p) {
	double weight, error;
	double x[4];
	double jac = 1.0;
	double wgt = 1.0;
	std::vector<Histogram*>* empty_histogram_set = new std::vector<Histogram*>();
	IntegrationVariablesMap v = MapToHyperCube(lo_variables, x, jac);
	for (auto& e : eg->GetEvents()) {
		std::tie(v["k1p"], v["phi1"], v["theta1"], v["theta2"], weight, error) = e;
		double integrand = lo::Hadronic(v, wgt, p, empty_histogram_set) * jac;
		std::cout << integrand << " " << error << std::endl;
	}
	delete empty_histogram_set;
}

void EvaluateNloEvents(EventGenerator* eg, Parameters* p) {
	double weight, error;
	double x[4];
	double jac = 1.0;
	double wgt = 1.0;
	std::vector<Histogram*>* empty_histogram_set = new std::vector<Histogram*>();
	IntegrationVariablesMap v = MapToHyperCube(lo_variables, x, jac);

	using namespace std::placeholders;
	Integrand z_integrand = std::bind(&nlo::HadronicZ, _1, _2, p, empty_histogram_set);
	Integrand j_integrand = std::bind(&nlo::Hadronic3, _1, _2, p, empty_histogram_set); 
	std::map<std::string, Integrand> z_integrands = { {"default", z_integrand} };
	std::map<std::string, Integrand> j_integrands = { {"default", j_integrand} };
	Integral* z_integral = new Integral(nlo_2_variables, z_integrands, v);
	Integral* j_integral = new Integral(nlo_3_variables, j_integrands, v);
	Integral* nlo_integral = new Integral(std::vector<Integral*>{z_integral, j_integral});
	for (auto& e : eg->GetEvents()) {
		std::tie(v["k1p"], v["phi1"], v["theta1"], v["theta2"], weight, error) = e;
		auto [integrand, sigma, chi] = NloHadronic(v, wgt, p, empty_histogram_set, nlo_integral);
		std::cout << integrand << " " << sigma << std::endl;
	}
	delete empty_histogram_set;
	delete nlo_integral;
}