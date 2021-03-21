
#include <sstream>
#include <iostream>
#include <algorithm>
#include "SysInfo.h"
#include "cross_sections.h"
#include "cxxopts.hpp"
#include "event_generator.h"
#include "matrix_element_method.h"



int main( int argc, char* argv[] ) {

	SysInfo sys_info;
	cxxopts::Options options(
		"ttbar@NLO", 
		"Hadronic ttbar production cross section at NLO accuracy");
	options.add_options()
		("ecms", "Collider energy", 
			cxxopts::value<double>()->default_value("13000.0"))
		("mtop", "Top quark mass", 
			cxxopts::value<double>()->default_value("173.2"))
		("mur", "Renormalization scale", 
			cxxopts::value<double>()->default_value("173.2"))
		("muf", "Factorization scale", 
			cxxopts::value<double>()->default_value("173.2"))
		("dynamic_scale", "Dynamic Scale", 
			cxxopts::value<bool>()->default_value("false"))
		("xmin", "Phase space slicing cut parameter", 
			cxxopts::value<double>()->default_value("0.0001"))
		("calculateXS", "Calculate cross section: lo, nlo or lo+nlo", 
			cxxopts::value<std::vector<std::string>>()->default_value(""))
		("ch", "Parton channels", 
			cxxopts::value<std::vector<std::string>>()->default_value("all"))
		("pdf", "LHAPDF set name", 
			cxxopts::value<std::string>()->default_value("CT10nlo"))
		("histogram", "Histogram creation strings", 
			cxxopts::value<std::vector<std::string>>()->default_value(""))
		("calls", "Number of calls to integrand", 
			cxxopts::value<int>()->default_value("100000"))
		("iterations", "Number of vegas iterations", 
			cxxopts::value<int>()->default_value("10"))
		("lo_events", "Number of LO Events", 
			cxxopts::value<int>()->default_value("0"))
		("nlo_events", "Number of NLO Events", 
			cxxopts::value<int>()->default_value("0"))
		("evaluate_lo_events", "Evaluate LO events", 
			cxxopts::value<std::string>()->default_value(""))
		("evaluate_nlo_events", "Evaluate NLO events", 
			cxxopts::value<std::string>()->default_value(""))
		("test", "Perform test routines", 
			cxxopts::value<bool>()->default_value("false"))
		("seed", "Seed for random number generator", 
			cxxopts::value<unsigned int>()->default_value("0"))
		("h,help", "Print usage")
		;
	auto result = options.parse(argc, argv);
	if (result.count("help")){
		std::cout << options.help() << std::endl;
		exit(0);
	}

	double ecms = result["ecms"].as<double>();
	double m = result["mtop"].as<double>();
	double mur = result["mur"].as<double>();
	double muf = result["muf"].as<double>();
	bool scale_is_dynamic = result["dynamic_scale"].as<bool>();
	double xmin = result["xmin"].as<double>();
	int calls = result["calls"].as<int>();
	int iterations = result["iterations"].as<int>();
	std::vector<std::string> perturbation_order = 
		result["calculateXS"].as<std::vector<std::string>>();
	std::vector<std::string> channels = 
		result["ch"].as<std::vector<std::string>>();
	std::string pdf_name = result["pdf"].as<std::string>();
	std::vector<std::string> histogram_strings = 
		result["histogram"].as<std::vector<std::string>>();

	int number_of_lo_events = result["lo_events"].as<int>();
	int number_of_nlo_events = result["nlo_events"].as<int>();
	std::string evaluate_lo_events =
		result["evaluate_lo_events"].as<std::string>();
	std::string evaluate_nlo_events =
		result["evaluate_nlo_events"].as<std::string>();
	bool running_test = result["test"].as<bool>();
	unsigned int iseed = result["seed"].as<unsigned int>();


	/* Initializing extern libraries */
	const int Nf = NF - 1;
	// const int FDH = 0;  /* Four-Dimensional-Helicity scheme                */
	const int HV_CDR = 1;  /* 't Hooft-Veltman scheme                         */
	bsyppttinit_(&m, &Nf, &HV_CDR);
	coupqcd_.gg[0] = -1.0, coupqcd_.gg[1] = -1.0, coupqcd_.g = 1.0; // 1.0 for gs 
	fermions_.fmass[10] = m;
	iseed = iseed == 0 ? (unsigned int)time(0): iseed;
	const int lxlev = 1;
	std::cout << "ranlxs is initialized with seed=" << iseed 
		<< " and lxlevel=" << lxlev << std::endl;
	rlxd_init(lxlev, iseed);

	parameter_sets.insert(
		{"default", 
		new Parameters (
			pdf_name, ecms, mur, muf, m, xmin, channels) });
	if (scale_is_dynamic) {
		parameter_sets.insert(
			{ "scale=2.0", 
			new Parameters(
				pdf_name, ecms, 2.0*mur, 2.0*muf, m, xmin, channels) });
		parameter_sets.insert(
			{ "scale=0.5", 
			new Parameters(
				pdf_name, ecms, 0.5*mur, 0.5*muf, m, xmin, channels) });
	}

	InitializeIntegrands(histogram_strings, ecms, m);
	if ( perturbation_order.size() != 0) {
		ExecuteIntegralsAndPrintResults(perturbation_order, iterations, calls);
	}

	if (running_test) {
		RunTestFunction(ecms, m);
		return 0;
	}



	if ( number_of_lo_events>0 ) {
		EventGenerator* eg = new EventGenerator(
			number_of_lo_events, "lo", parameter_sets["default"]);
		eg->Print();
		delete eg;
	}
	if (number_of_nlo_events > 0) {
		EventGenerator* eg = new EventGenerator(
			number_of_nlo_events, "nlo", parameter_sets["default"]);
		eg->Print();
		delete eg;
	}
	if ( evaluate_lo_events != "") {

		std::vector<Event> events = EventGenerator(evaluate_lo_events).GetEvents();
		EvaluateLoEvents(events, parameter_sets["default"]);
		}
	if (evaluate_nlo_events != "") {
		std::vector<Event> events = EventGenerator(evaluate_nlo_events).GetEvents();
		EvaluateNloEvents(events, parameter_sets["default"]);
	}
	

	for (const auto& it: histogram_sets) {
		for ( const auto& ij: *(it.second) ) {
			delete ij;
		}
		it.second->clear();
		delete it.second;
	}
	histogram_sets.clear();

	for (const auto& it: integrals) {
		for (const auto& ij: *(it.second)) {
			delete ij;
		}
		it.second->clear();
		delete it.second;
	}
	integrals.clear();


	for (const auto& it: parameter_sets) {
		delete it.second;
	}
	parameter_sets.clear();

	return 0;
}
