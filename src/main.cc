
#include <sstream>
#include <iostream>
#include <algorithm>
#include "SysInfo.h"
#include "cross_sections.h"
#include "n_dimensional_integral.h"
#include "cxxopts.hpp"



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
	std::vector<std::string> po = result["po"].as<std::vector<std::string>>();
	std::vector<std::string> channels = result["ch"].as<std::vector<std::string>>();
	std::string pdf_name = result["pdf"].as<std::string>();
	std::map<std::string, Parameters*> parameter_sets;
	std::map<std::string, std::vector<Histogram*>*> histogram_sets;
	std::vector<std::string> histogram_strings = {
		"M(top atop) 50 300.0 1000.0" };/* ,
		"PT(top) 50 0.0 500.0",
		"PT(atop) 50 0.0 500.0",
		"Y(top) 50 -5.0 5.0",
		"E(top) 50 0.0 1000.0",
		"PHI(top) 50 -4.0 4.0",
		"ETA(top) 50 -10.0 10.0",
		"ETA(atop) 50 -10.0 10.0" };*/



	/* Initializing extern libraries */
	const int Nf = NF - 1;
	// const int FDH = 0;  /* Four-Dimensional-Helicity scheme                    */
	const int HV_CDR = 1;  /* 't Hooft-Veltman scheme                           */
	bsyppttinit_(&m, &Nf, &HV_CDR);
	coupqcd_.gg[0] = -1.0, coupqcd_.gg[1] = -1.0, coupqcd_.g = 1.0; // 1.0 for gs 
	fermions_.fmass[10] = m;
	unsigned int iseed = 3720758; // I can also take the cpu time here.
	const int lxlev = 1;
	rlxd_init(lxlev, iseed);
	



	std::map<std::string, std::pair<double, double>> lo_variables{
		{"E1",     std::make_pair(m, ecms/2.0)},
		{"phi1",   std::make_pair(-M_PI, M_PI)},
		{"theta1", std::make_pair(0.0, M_PI)},
		{"theta2", std::make_pair(0.0, M_PI)}
	};

	std::map<std::string, std::pair<double, double>> nlo_2_variables{
		{"E1",     std::make_pair(m, ecms / 2.0)},
		{"phi1",   std::make_pair(-M_PI, M_PI)},
		{"theta1", std::make_pair(0.0, M_PI)},
		{"theta2", std::make_pair(0.0, M_PI)},
		{"z",      std::make_pair(0.0,1.0)}
	};

	std::map<std::string, std::pair<double, double>> nlo_3_variables{
		{"E1",     std::make_pair(m, ecms / 2.0)},
		{"phi1",   std::make_pair(-M_PI, M_PI)},
		{"theta1", std::make_pair(0.0, M_PI)},
		{"theta2", std::make_pair(0.0, M_PI)},
		{"E3",     std::make_pair(0.0,ecms)},
		{"phi3",   std::make_pair(-M_PI,M_PI)},
		{"theta3", std::make_pair(0.0, M_PI)}
	};


	std::map<std::string, Integrand> lo_integrands, nlo_2_integrands, nlo_3_integrands;


	using namespace std::placeholders;
	//std::bind always copies the argument.  to force it pass by reference i used std::ref
	parameter_sets.insert(
		{"scale=1.0", new Parameters (pdf_name, ecms, mur, muf, m, xmin, channels) });
	if (scale_is_dynamic) {
		parameter_sets.insert(
			{ "scale=2.0", new Parameters(pdf_name, ecms, 2.0*mur, 2.0*muf, m, xmin, channels) });
		parameter_sets.insert(
			{ "scale=0.5", new Parameters(pdf_name, ecms, 0.5*mur, 0.5*muf, m, xmin, channels) });
	}

	for (auto it = parameter_sets.begin(); it != parameter_sets.end(); ++it) {

		std::vector<Histogram*> *hists = new std::vector<Histogram*>();
		for (auto ij = histogram_strings.begin(); ij != histogram_strings.end(); ++ij) {
			hists->push_back(new Histogram(*ij, it->second));
		}

		histogram_sets.insert({ it->first, hists });
		lo_integrands[it->first]    = std::bind(&lo::Hadronic, _1, _2,   it->second, hists);
		nlo_2_integrands[it->first] = std::bind(&nlo::Hadronic2, _1, _2, it->second, hists);
		nlo_3_integrands[it->first] = std::bind(&nlo::Hadronic3, _1, _2, it->second, hists);
	}
	

	std::map<std::string, std::vector<Integral*>*> integrals = { 
		{"lo",  new std::vector<Integral*>{ new Integral(lo_variables,    lo_integrands) }},
		{"nlo", new std::vector<Integral*>{ new Integral(nlo_2_variables, nlo_2_integrands),
	                                        new Integral(nlo_3_variables, nlo_3_integrands)}}
};

	for (auto it = po.begin(); it != po.end(); ++it) {
		for (auto ij = integrals.find(*it)->second->begin(); ij != integrals.find(*it)->second->end(); ++ij) {
			(*ij)->ExecuteVegas(1, 10, 1000, 1);
		}
	}

	for (auto it = histogram_sets.begin(); it != histogram_sets.end(); ++it) {
		std::cout << "Histogram for parameter set: " << it->first << std::endl;
		for (auto ij = it->second->begin(); ij != it->second->end(); ++ij) {
			(*ij)->Print();
		}
	}

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
