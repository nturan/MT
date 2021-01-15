
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
	cxxopts::Options options("MyProgram", "One line description of MyProgram");
	options.add_options()
		("ecms", "Collider energy", cxxopts::value<double>()->default_value("13000.0"))
		("mtop", "Top quark mass", cxxopts::value<double>()->default_value("173.2"))
		("mur", "Renormalization scale", cxxopts::value<double>()->default_value("173.2"))
		("muf", "Factorization scale", cxxopts::value<double>()->default_value("173.2"))
		("dynamic_scale", "Dynamic Scale", cxxopts::value<bool>()->default_value("false"))
		("xmin", "Phase space slicing cut parameter", cxxopts::value<double>()->default_value("0.0005"))
		("po", "Perturbation order", cxxopts::value<std::vector<std::string>>()->default_value("lo"))
		("ch", "Parton channels", cxxopts::value<std::vector<std::string>>()->default_value("all"))
		("pdf", "LHAPDF set name", cxxopts::value<std::string>()->default_value("CT10nlo"))
		("histogram", "Histogram creation strings", cxxopts::value<std::vector<std::string>>()->default_value(""))
		("calls", "Number of calls to integrand", cxxopts::value<int>()->default_value("100000"))
		("iterations", "Number of vegas iterations", cxxopts::value<int>()->default_value("10"))
		("lo_events", "Number of LO Events", cxxopts::value<int>()->default_value("0"))
		("nlo_events", "Number of NLO Events", cxxopts::value<int>()->default_value("0"))
		("events_file", "Events from file", cxxopts::value<std::string>()->default_value(""))
		("evaluate_events", "Evaluate events", cxxopts::value<bool>()->default_value("false"))
		("h,help", "Print usage")
		;
	auto result = options.parse(argc, argv);
	if (result.count("help"))
	{
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
	std::vector<std::string> perturbation_order = result["po"].as<std::vector<std::string>>();
	std::vector<std::string> channels = result["ch"].as<std::vector<std::string>>();
	std::string pdf_name = result["pdf"].as<std::string>();
	std::vector<std::string> histogram_strings = result["histogram"].as<std::vector<std::string>>();

	int number_of_lo_events = result["lo_events"].as<int>();
	int number_of_nlo_events = result["nlo_events"].as<int>();
	std::string events_file = result["events_file"].as<std::string>();
	bool evaluating_events = result["evaluate_events"].as<bool>();




	/* Initializing extern libraries */
	const int Nf = NF - 1;
	// const int FDH = 0;  /* Four-Dimensional-Helicity scheme                    */
	const int HV_CDR = 1;  /* 't Hooft-Veltman scheme                           */
	bsyppttinit_(&m, &Nf, &HV_CDR);
	coupqcd_.gg[0] = -1.0, coupqcd_.gg[1] = -1.0, coupqcd_.g = 1.0; // 1.0 for gs 
	fermions_.fmass[10] = m;
	unsigned int iseed = 3720758; // I can also take the cpu time here.
	iseed = (unsigned int)time(0);
	const int lxlev = 1;
	rlxd_init(lxlev, iseed);

	parameter_sets.insert(
		{"default", new Parameters (pdf_name, ecms, mur, muf, m, xmin, channels) });
	if (scale_is_dynamic) {
		parameter_sets.insert(
			{ "scale=2.0", new Parameters(pdf_name, ecms, 2.0*mur, 2.0*muf, m, xmin, channels) });
		parameter_sets.insert(
			{ "scale=0.5", new Parameters(pdf_name, ecms, 0.5*mur, 0.5*muf, m, xmin, channels) });
	}

	InitializeIntegrands(histogram_strings, ecms, m);
	ExecuteIntegralsAndPrintResults(perturbation_order, iterations, calls);

	EventGenerator* eg;

	if ( number_of_lo_events>0 ) {
		eg = new EventGenerator(number_of_lo_events, "lo", parameter_sets["default"]);
		eg->Print();
	}
	if (number_of_nlo_events > 0) {
		eg = new EventGenerator(number_of_nlo_events, "nlo", parameter_sets["default"]);
		eg->Print();
	}
	if (events_file != "") {
		eg = new EventGenerator(events_file);
	}

	if (evaluating_events) {
		Parameters* p[20];
		for (int i = 0; i < 20; ++i) {
			p[i] = new Parameters(pdf_name, ecms, mur, muf, 160.0 + i * 1.25, xmin, channels);
			EvaluateEvents(eg, p[i]);
			delete p[i];
		}
	}
	
	delete eg;
	
	

	for (auto it = parameter_sets.begin(); it != parameter_sets.end(); ++it) {
		delete it->second;
	}
    parameter_sets.clear();
	for (auto it = histogram_sets.begin(); it != histogram_sets.end(); ++it) {
		for (auto ij = it->second->begin(); ij != it->second->end(); ++ij) {
			delete *ij;
		}
		it->second->clear();
		delete it->second;
	}
	histogram_sets.clear();

	for (auto it = integrals.begin(); it != integrals.end(); ++it) {
		for (auto ij = it->second->begin(); ij != it->second->end(); ++ij) {
			delete *ij;
		}
		it->second->clear();
		delete it->second;
	}
	integrals.clear();

	return 0;
}
