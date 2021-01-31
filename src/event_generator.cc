#include "event_generator.h"

double LoMaxIntegrandFinderIntegral(
	std::map<std::string, double> v, double& wgt, Parameters* p, 
	std::vector<Histogram*>* histograms, double *max_integrand){
	double result = lo::Hadronic(v, wgt, p, histograms);
	if (result > *max_integrand) { *max_integrand = result; }
	return result;
}

std::tuple<double, double, double> NloHadronic(
	std::map<std::string, double> v, double& wgt, Parameters* p,
	std::vector<Histogram*>* histograms, Integral* nlo_integral) {
	
	double x[4], jac = 1.0;
	MapToHyperCube(lo_variables, x, jac);
	double dLO = lo::Hadronic(v, wgt, p, histograms) * jac;
	double dNLOconst = nlo::HadronicConstZ(v, wgt, p, histograms) * jac;
	double accuracy = 1.0E-2;
	double calls = 10000;
	nlo_integral->constants1_ = v;
	nlo_integral->constants2_ = v;
	auto [dNLOint, integration_error, chi] =
		nlo_integral->ExecuteVegas(0, 30, (int)calls, 1);
	double dNLO = dNLOint + dNLOconst + dLO;
	/* since other contributions are practically without error */
	double error = integration_error; 
	double relative_nlo_correction = std::abs( (dNLO - dLO) / dLO );
	double relative_nlo_error = std::abs(error / dNLO);
	while ( relative_nlo_error / relative_nlo_correction > accuracy ) {
		calls *= 1.5;
		auto [dNLOint, integration_error, chi] = 
			nlo_integral->ExecuteVegas(2, 1, (int)calls, 1);
		dNLO = dNLOint + dNLOconst + dLO;
		error = integration_error;
		relative_nlo_correction = std::abs((dNLO - dLO) / dLO);
		relative_nlo_error = std::abs(error / dNLO);
		std::cout << "     relative nlo error: " 
			<< relative_nlo_error << std::endl;
		std::cout << "relative nlo correction: " 
			<< relative_nlo_correction << std::endl;
	}

	return std::tie(dNLO, error, chi);
}

EventGenerator::EventGenerator(){}

EventGenerator::EventGenerator(std::string file_name) {
	events.clear();
	std::ifstream in_file(file_name);
	while (!in_file.eof()) {
		double E1, phi1, theta1, theta2, weight, error;
		in_file >> E1;
		in_file >> phi1;
		in_file >> theta1;
		in_file >> theta2;
		in_file >> weight;
		in_file >> error;
		events.push_back(std::tie(E1, phi1, theta1, theta2, weight, error));
	}
}


EventGenerator::EventGenerator(
	unsigned int number, std::string perturbation_order, Parameters* p){
	if (perturbation_order == "lo") {
		GenerateLOEvents(number, p);
	}
	else if (perturbation_order == "nlo") {
		GenerateNLOEvents(number, p);
	}
}

void EventGenerator::GenerateLOEvents(unsigned int number, Parameters* p) {
	/* Determine max integrand */
	using namespace std::placeholders;
	double max_integrand = 0.0;
	std::vector<Histogram*>* empty_histogram_set = 
		new std::vector<Histogram*>();
	Integrand new_integrand = 
		std::bind(&LoMaxIntegrandFinderIntegral, _1, _2, p, 
		empty_histogram_set, &max_integrand);
	std::map<std::string, Integrand> integrands = { 
		{"max finder", new_integrand} };
	Integral* max_integrand_finder_integral = 
		new Integral(lo_variables, integrands);
	std::cout << "Now executing max integrand finder integral..." << std::endl;
	max_integrand_finder_integral->ExecuteVegas(0, 10, 100000, 0);
	//here i want to calculate jac by hand and multiple max_integrand with it.
	double x[4], jac=1.0;
	MapToHyperCube(lo_variables, x, jac);
	max_integrand *= jac;
	std::cout << "DONE! Integrand maximum is: " << max_integrand << std::endl;

	int attempt_counter = 0;
	events.clear();
	while ( events.size() < number ) {
		++attempt_counter;
		double x[4];
		ranlxd(x, 4);
		double jac = 1.0;
		double wgt = 1.0;
		IntegrationVariablesMap v = MapToHyperCube(lo_variables, x, jac);
		double integrand = lo::Hadronic(v, wgt, p, empty_histogram_set)*jac;
		double y[1];
		ranlxd(y, 1);
		double rho = y[0] * max_integrand * 1.1;
		if (integrand > rho) {
			double error = 0.0;
			events.push_back(
				std::tie(v["k1p"], v["phi1"], v["theta1"], v["theta2"], 
					integrand, error));
		}
	}
	double efficiency = 100.0 * ((double)number) / attempt_counter;
	std::cout << "Event generation efficiency is: " << efficiency 
		<< "%." << std::endl;
	delete max_integrand_finder_integral;
	delete empty_histogram_set;
}

void EventGenerator::GenerateNLOEvents(unsigned int number, Parameters* p) {

	/*
	  Initialize integrands together with other parts of calculation
	
	*/

	std::vector<Histogram*>* empty_histogram_set = 
		new std::vector<Histogram*>();
	int attempt_counter = 0;
	std::cout << 
		"Generating LO events to sample from their distribution" << std::endl;
	GenerateLOEvents(number*3, p);
	std::vector<Event> lo_events = events;
	events.clear();


	using namespace std::placeholders;
	Integrand z_integrand = std::bind(
		&nlo::HadronicZ, _1, _2, p, empty_histogram_set);
	Integrand j_integrand = std::bind(
		&nlo::Hadronic3, _1, _2, p, empty_histogram_set);

	std::map<std::string, Integrand> z_integrands = { {"default", z_integrand} };
	std::map<std::string, Integrand> j_integrands = { {"default", j_integrand} };

	IntegrationVariablesMap v{ 
		{"k1p", 0.0}, 
		{"phi1", 0.0}, 
		{"theta1", 0.0}, 
		{"theta2", 0.0} };

	double x[4], jac=1.0;
	MapToHyperCube(lo_variables, x, jac);

	Integral* z_integral = new Integral(nlo_2_variables, z_integrands, v);
	Integral* j_integral = new Integral(nlo_3_variables, j_integrands, v);
	Integral* nlo_integral = 
		new Integral(std::vector<Integral*>{z_integral, j_integral}, 
			std::make_pair(0.01, 0.99));
	for (auto e = lo_events.begin(); e != lo_events.end(); ++e) {
		++attempt_counter;
		double lo_weight, lo_error;
		double wgt = 1.0;
		std::tie(v["k1p"], v["phi1"], v["theta1"], v["theta2"], 
			lo_weight, lo_error) = *e;

//		double dNLO, err, chi;
		auto [dNLO, err, chi] = 
			NloHadronic(v, wgt, p, empty_histogram_set, nlo_integral);
		
		double y[1];
		ranlxd(y, 1);
		double rho = y[0] * lo_weight * 2.0;
		if (dNLO > rho) {
			events.push_back(
				std::tie(v["k1p"], v["phi1"], v["theta1"], v["theta2"], 
				dNLO, err));
		}
		if (events.size() >= number) {
			break;
		}
	}

	delete nlo_integral;
	delete z_integral;
	delete j_integral;
	delete empty_histogram_set;
	double efficiency = 100.0 * ((double)number) / attempt_counter;
	std::cout << "Event generation efficiency is: " 
		<< efficiency << "%." << std::endl;
}

void EventGenerator::Print() {
	std::cout << "#START: Events" << std::endl;
	for (auto e = events.begin(); e != events.end(); ++e) {
		double E1, phi1, theta1, theta2, weight, error;
		std::tie(E1, phi1, theta1, theta2, weight, error) = *e;
		std::cout << E1 << " " 
			<< phi1 << " " 
			<< theta1 << " " 
			<< theta2 << " " 
			<< weight  << " " << error << std::endl;
	}
	std::cout << "#END: Events" << std::endl;
}

std::vector<Event> EventGenerator::GetEvents() {
	return events;
}


