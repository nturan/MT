
#include <sstream>
#include <iostream>
#include <algorithm>
#include "SysInfo.h"
#include "cross_sections.h"
#include "n_dimensional_integral.h"
#include "cxxopts.hpp"

void InitializeChannels() {

	using namespace std::placeholders;
	sigma::leading_order::partonic::Born["gg"] = std::bind(&sigma::leading_order::partonic::gg, _1, _2);
	sigma::leading_order::partonic::Born["qqb"] = std::bind(&sigma::leading_order::partonic::qqb, _1, _2);
	sigma::leading_order::partonic::Born["qbq"] = std::bind(&sigma::leading_order::partonic::qbq, _1, _2);
	sigma::leading_order::partonic::Born["qg"] = std::bind(&NoContribution, _1, _2);
	sigma::leading_order::partonic::Born["gq"] = std::bind(&NoContribution, _1, _2);
	sigma::leading_order::partonic::Born["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::leading_order::partonic::Born["gqb"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["gg"] = std::bind(&sigma::next_to_leading_order::partonic::hard::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::hard::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::hard::qbq, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["qg"] = std::bind(&sigma::next_to_leading_order::partonic::hard::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["gq"] = std::bind(&sigma::next_to_leading_order::partonic::hard::gq, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["qbg"] = std::bind(&sigma::next_to_leading_order::partonic::hard::qbg, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["gqb"] = std::bind(&sigma::next_to_leading_order::partonic::hard::gqb, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["gg"] = std::bind(&sigma::next_to_leading_order::partonic::soft::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::soft::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::soft::qbq, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["qg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["gq"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["gqb"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["gg"] = std::bind(&sigma::next_to_leading_order::partonic::virt::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::virt::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::virt::qbq, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["qg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["gq"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["gqb"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["gg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["qg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["gq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_right_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["qbg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["gqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_right_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["gg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["qg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_right_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["gq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["qbg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_right_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["gqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["gg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_1::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_1::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_1::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["qg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["gq"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["gqb"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["gg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_0::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_0::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_0::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["qg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["gq"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["gqb"] = std::bind(&NoContribution, _1, _2);


}


int main( int argc, char* argv[] ) {

	SysInfo sys_info;
	cxxopts::Options options("MyProgram", "One line description of MyProgram");
	options.add_options()
		("ecms", "Collider energy", cxxopts::value<double>()->default_value("13000.0"))
		("mtop", "Top quark mass", cxxopts::value<double>()->default_value("173.2"))
		("mur", "Renormalization scale", cxxopts::value<double>()->default_value("173.2"))
		("muf", "Factorization scale", cxxopts::value<double>()->default_value("173.2"))
		("xmin", "Phase space slicing cut parameter", cxxopts::value<double>()->default_value("0.0005"))
		("po", "Perturbation order", cxxopts::value<std::string>()->default_value("lo"))
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
	double xmin = result["xmin"].as<double>();
	std::string po = result["po"].as<std::string>();
	std::vector<std::string> channels = result["ch"].as<std::vector<std::string>>();
	std::string pdf_name = result["pdf"].as<std::string>();
	Parameters parameters(pdf_name);
	parameters.SetColliderEnergy(ecms);
	parameters.SetTopQuarkMass(m);
	parameters.SetFactorizationScale(muf);
	parameters.SetRenormalizationScale(mur);
	parameters.SetCutParameter(xmin);
	if (std::find(channels.begin(), channels.end(), "all") != channels.end()) {
		parameters.channels.push_back("gg");
		parameters.channels.push_back("qqb");
		parameters.channels.push_back("qbq");
		parameters.channels.push_back("qg");
		parameters.channels.push_back("gq");
		parameters.channels.push_back("qbg");
		parameters.channels.push_back("gqb");
	}
	else {
		for (auto it = channels.begin(); it != channels.end(); ++it) {
			parameters.channels.push_back(*it);
		}
	}
	
	InitializeChannels();

//	std::vector<Histogram> histograms = { 
//		Histogram("M[top atop]", std::make_pair(300.0, 1000.0), 50, &calcInvariantMass, std::vector<std::string>{"top", "atop"}, parameters),
//	    Histogram("PT[top]", std::make_pair(0.0, 500.0), 50, &calcTransversalMomentum, std::vector<std::string>{"top"}, parameters),
//	    Histogram("PT[atop]", std::make_pair(0.0, 500.0), 50, &calcTransversalMomentum, std::vector<std::string>{"atop"}, parameters),
//	    Histogram("Y[top]", std::make_pair(-5.0, 5.0), 50, &calcRapidity, std::vector<std::string>{"top"}, parameters),
//	    Histogram("E[top]", std::make_pair(0.0, 1000.0), 50, &calcEnergy, std::vector<std::string>{"top"}, parameters),
//	    Histogram("PHI[top]", std::make_pair(-4.0, 4.0), 50, &calcAzimuth, std::vector<std::string>{"top"}, parameters),
//	    Histogram("ETA[top]", std::make_pair(-10.0, 10.0), 50, &calcPseudoRapidity, std::vector<std::string>{"top"}, parameters),
//		Histogram("ETA[atop]", std::make_pair(-10.0, 10.0), 50, &calcPseudoRapidity, std::vector<std::string>{"atop"}, parameters) };
	std::vector<Histogram> histograms = {
		Histogram("M(top atop) 50 300.0 1000.0", parameters),
		Histogram("PT(top) 20 0.0 500.0", parameters)
	};

	const int Nf = NF - 1;
	// const int FDH = 0;  /* Four-Dimensional-Helicity scheme                    */
	const int HV_CDR = 1;  /* 't Hooft-Veltman scheme                           */
	bsyppttinit_(&m, &Nf, &HV_CDR);


	std::map<std::string, Integrand> my_integrands;
	std::map<std::string, std::pair<double, double>> my_variables;
	my_variables["E1"]     = std::pair<double, double>{ m, ecms/2.0 };
	my_variables["phi1"]   = std::pair<double, double>{ -M_PI, M_PI };
	my_variables["theta1"] = std::pair<double, double>{ 0.0, M_PI };
	my_variables["theta2"] = std::pair<double, double>{ 0.0, M_PI };
//	my_variables["E3"] = std::pair<double, double>{ 0.0, ecms };
//	my_variables["phi3"] = std::pair<double, double>{ -M_PI, M_PI };
//	my_variables["theta3"] = std::pair<double, double>{ 0.0, M_PI };

	using namespace std::placeholders;
	//std::bind always copies the argument.  to force it pass by reference i used std::ref
	my_integrands["mur=m"] = std::bind(&sigma::leading_order::Hadronic, _1, _2, parameters, std::ref(histograms));
	coupqcd_.gg[0] = -1.0, coupqcd_.gg[1] = -1.0, coupqcd_.g = 1.0; // 1.0 for gs 
	fermions_.fmass[10] = 173.2;

	unsigned int iseed = 3720758;
	const int lxlev = 1;
	rlxd_init(lxlev, iseed);

	std::cout << "chosen perturbation order: " << po << std::endl;
	Integral my_integral(my_variables, my_integrands);
	my_integral.ExecuteVegas(1, 10, 100000, 1);
	for (auto it = histograms.begin(); it != histograms.end(); ++it) {
		(*it).Print();
	}


	return 0;
}
