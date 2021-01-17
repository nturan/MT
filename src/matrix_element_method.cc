#include "matrix_element_method.h"

void EvaluateEvents(EventGenerator* eg, Parameters* p) {
	double weight;
	double x[4];
	double jac = 1.0;
	double wgt = 1.0;
	double m = p->GetTopQuarkMass();
	std::vector<Histogram*>* empty_histogram_set = new std::vector<Histogram*>();
	IntegrationVariablesMap v = MapToHyperCube(lo_variables, x, jac);
	for (auto& e : eg->GetEvents()) {
		std::tie(v["E1"], v["phi1"], v["theta1"], v["theta2"], weight) = e;
		double integrand = lo::Hadronic(v, wgt, p, empty_histogram_set) * jac;
		std::cout << m << " " << v["E1"] << " " << v["phi1"] << " " << v["theta1"] << " " << v["theta2"] << " " << integrand << std::endl;
	}
	delete empty_histogram_set;
}