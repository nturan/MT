
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
	std::map<std::string, Parameters> parameters;
	std::map<std::string, std::vector<Histogram*>*> histograms;




	const int Nf = NF - 1;
	// const int FDH = 0;  /* Four-Dimensional-Helicity scheme                    */
	const int HV_CDR = 1;  /* 't Hooft-Veltman scheme                           */
	bsyppttinit_(&m, &Nf, &HV_CDR);

	coupqcd_.gg[0] = -1.0, coupqcd_.gg[1] = -1.0, coupqcd_.g = 1.0; // 1.0 for gs 
	fermions_.fmass[10] = m;
	

	//for (auto it = parameters.begin(); it != parameters.end(); ++it) {
	//	histograms[it->first] = std::vector<Histogram>{
	//		Histogram("M(top atop) 50 300.0 1000.0", it->second),
	//		Histogram("PT(top) 50 0.0 500.0", it->second),
	//		Histogram("PT(atop) 50 0.0 500.0", it->second),
	//		Histogram("Y(top) 50 -5.0 5.0", it->second),
	//		Histogram("E(top) 50 0.0 1000.0", it->second),
	//		Histogram("PHI(top) 50 -4.0 4.0", it->second),
	//		Histogram("ETA(top) 50 -10.0 10.0", it->second),
	//		Histogram("ETA(atop) 50 -10.0 10.0", it->second) };
	//}



	unsigned int iseed = 3720758;
	const int lxlev = 1;
	rlxd_init(lxlev, iseed);
	std::map<std::string, std::pair<double, double>> leading_order_variables{
		{"E1",     std::make_pair(m, ecms/2.0)},
		{"phi1",   std::make_pair(-M_PI, M_PI)},
		{"theta1", std::make_pair(0.0, M_PI)},
		{"theta2", std::make_pair(0.0, M_PI)}
	};

	std::map<std::string, std::pair<double, double>> next_to_leading_order_twobody_variables{
		{"E1",     std::make_pair(m, ecms / 2.0)},
		{"phi1",   std::make_pair(-M_PI, M_PI)},
		{"theta1", std::make_pair(0.0, M_PI)},
		{"theta2", std::make_pair(0.0, M_PI)},
		{"z",      std::make_pair(0.0,1.0)}
	};

	std::map<std::string, std::pair<double, double>> next_to_leading_order_threebody_variables{
		{"E1",     std::make_pair(m, ecms / 2.0)},
		{"phi1",   std::make_pair(-M_PI, M_PI)},
		{"theta1", std::make_pair(0.0, M_PI)},
		{"theta2", std::make_pair(0.0, M_PI)},
		{"E3",     std::make_pair(0.0,ecms)},
		{"phi3",   std::make_pair(-M_PI,M_PI)},
		{"theta3", std::make_pair(0.0, M_PI)}
	};


	std::map<std::string, Integrand> 
		leading_order_integrands, 
		next_to_leading_order_twobody_integrands, 
		next_to_leading_order_threebody_integrands;


	using namespace std::placeholders;
	//std::bind always copies the argument.  to force it pass by reference i used std::ref
	

	Parameters p1(pdf_name, ecms, mur, muf, m, xmin, channels);
	Parameters p2(pdf_name, ecms, 2*mur, 2*muf, m, xmin, channels);
	Parameters p3(pdf_name, ecms, 0.5 * mur, 0.5 * muf, m, xmin, channels);


	parameters.insert(std::pair<std::string, Parameters>("one", p1));
	parameters.insert(std::pair<std::string, Parameters>("two", p2));
	parameters.insert(std::pair<std::string, Parameters>("three", p3));

	for (auto it = parameters.begin(); it != parameters.end(); ++it) {

		std::vector<Histogram*> *hists = new std::vector<Histogram*>();
		hists->push_back(new Histogram("M(top atop) 50 300.0 1000.0", it->second));
		hists->push_back(new Histogram("PT(top) 50 0.0 500.0", it->second));

		histograms.insert(std::pair<std::string, std::vector<Histogram*>*>(it->first, hists));
		leading_order_integrands[it->first] = std::bind(&sigma::leading_order::Hadronic, _1, _2, std::ref(it->second), hists);
		next_to_leading_order_twobody_integrands[it->first] = std::bind(&sigma::next_to_leading_order::Hadronic2, _1, _2, std::ref(it->second), hists);
		next_to_leading_order_threebody_integrands[it->first] = std::bind(&sigma::next_to_leading_order::Hadronic3, _1, _2, std::ref(it->second), hists);
	}
	

	std::map<std::string, std::vector<Integral*>> integrals;
	Integral *leading_order_integral = new Integral(leading_order_variables, leading_order_integrands);
	Integral *next_to_leading_order_twobody_integral = new Integral(next_to_leading_order_twobody_variables, next_to_leading_order_twobody_integrands);
	Integral *next_to_leading_order_threebody_integral = new Integral(next_to_leading_order_threebody_variables, next_to_leading_order_threebody_integrands);

	integrals.insert(std::pair<std::string, std::vector<Integral*>>("lo",  std::vector<Integral*>{ leading_order_integral }));
	integrals.insert(std::pair<std::string, std::vector<Integral*>>("nlo", std::vector<Integral*>{ next_to_leading_order_twobody_integral,
		                                                                                           next_to_leading_order_threebody_integral}));

	for (auto it = po.begin(); it != po.end(); ++it) {
		for (auto ij = integrals.find(*it)->second.begin(); ij != integrals.find(*it)->second.end(); ++ij) {
			(*ij)->ExecuteVegas(1, 10, 1000, 1);
		}
	}

	for (auto it = histograms.begin(); it != histograms.end(); ++it) {
		std::cout << "Histogram for parameter set: " << it->first << std::endl;
		for (auto ij = it->second->begin(); ij != it->second->end(); ++ij) {
			(*ij)->Print();
			delete *ij;
		}
		it->second->clear();
		delete it->second;
	}
	parameters.clear();
	histograms.clear();
	delete leading_order_integral;
	delete next_to_leading_order_twobody_integral;
	delete next_to_leading_order_threebody_integral;

	return 0;
}
