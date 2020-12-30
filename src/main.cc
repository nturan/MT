
#include <iostream>
#include "SysInfo.h"
#include "cross_sections.h"
#include "n_dimensional_integral.h"


int main()
{
	SysInfo sys_info;
	Parameters parameters("CT10nlo");
	parameters.SetColliderEnergy(13000);
	parameters.SetTopQuarkMass(173.2);
	parameters.SetFactorizationScale(173.2);
	parameters.SetRenormalizationScale(173.2);
	parameters.channels.push_back("gg" );
	parameters.channels.push_back("qqb");
	parameters.channels.push_back("qbq");
	parameters.channels.push_back("qg" );
	parameters.channels.push_back("gq" );
	parameters.channels.push_back("qbg");
	parameters.channels.push_back("gqb");



	using namespace std::placeholders;
	sigma::lo::partonic::Born["gg"]  = std::bind(&sigma::lo::partonic::gg,  _1, _2);
	sigma::lo::partonic::Born["qqb"] = std::bind(&sigma::lo::partonic::qqb, _1, _2);
	sigma::lo::partonic::Born["qbq"] = std::bind(&sigma::lo::partonic::qbq, _1, _2);
	sigma::lo::partonic::Born["qg"]  = std::bind(&sigma::lo::partonic::NoContribution, _1, _2);
	sigma::lo::partonic::Born["gq"]  = std::bind(&sigma::lo::partonic::NoContribution, _1, _2);
	sigma::lo::partonic::Born["qbg"] = std::bind(&sigma::lo::partonic::NoContribution, _1, _2);
	sigma::lo::partonic::Born["gqb"] = std::bind(&sigma::lo::partonic::NoContribution, _1, _2);

	std::map<std::string, Integrand> my_integrands;
	std::map<std::string, std::pair<double, double>> my_variables;
	my_variables["E1"]     = std::pair<double, double>{ 173.2, 13000.0 / 2.0 };
	my_variables["phi1"]   = std::pair<double, double>{ -M_PI, M_PI };
	my_variables["theta1"] = std::pair<double, double>{ 0.0, M_PI };
	my_variables["theta2"] = std::pair<double, double>{ 0.0, M_PI };


	my_integrands["mur=m"] = std::bind(&sigma::lo::Hadronic, _1, _2, parameters);
	coupqcd_.gg[0] = -1.0, coupqcd_.gg[1] = -1.0, coupqcd_.g = 1.0; // 1.0 for gs 
	fermions_.fmass[10] = 173.2;

	Integral my_integral(my_variables, my_integrands);
	my_integral.ExecuteVegas(1, 10, 100000, 1);

	return 0;
}
