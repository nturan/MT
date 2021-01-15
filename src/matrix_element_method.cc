#include "matrix_element_method.h"

void EvaluateEvents(EventGenerator* eg, Parameters* p) {
	double weight;
	double x[4];
	double jac = 1.0;
	double wgt = 1.0;
	double m = p->GetTopQuarkMass();
	std::string file_ending = ".mdat";
	std::ostringstream name;
	name.precision(2);
	name << std::fixed << m;
	std::ofstream data_file(name.str() + file_ending);
	std::vector<Histogram*>* empty_histogram_set = new std::vector<Histogram*>();
	IntegrationVariablesMap v = MapToHyperCube(lo_variables, x, jac);
	for (auto& e : eg->GetEvents()) {
		std::tie(v["E1"], v["phi1"], v["theta1"], v["theta2"], weight) = e;
		double integrand = lo::Hadronic(v, wgt, p, empty_histogram_set) * jac;
		data_file << v["E1"] << " " << v["phi1"] << " " 
			      << v["theta1"] << " " << v["theta2"] << " " << integrand << std::endl;
	}
	data_file.close();
	delete empty_histogram_set;
}