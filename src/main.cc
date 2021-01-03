
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
	parameters.SetCutParameter(0.0005);
	parameters.channels.push_back("gg" );
	parameters.channels.push_back("qqb");
	parameters.channels.push_back("qbq");
	parameters.channels.push_back("qg" );//correct here
	parameters.channels.push_back("gq" );//correct here
	parameters.channels.push_back("qbg");//correct here
	parameters.channels.push_back("gqb");//correct here



	using namespace std::placeholders;
	sigma::leading_order::partonic::Born["gg"]  = std::bind(&sigma::leading_order::partonic::gg,  _1, _2);
	sigma::leading_order::partonic::Born["qqb"] = std::bind(&sigma::leading_order::partonic::qqb, _1, _2);
	sigma::leading_order::partonic::Born["qbq"] = std::bind(&sigma::leading_order::partonic::qbq, _1, _2);
	sigma::leading_order::partonic::Born["qg"]  = std::bind(&NoContribution, _1, _2);
	sigma::leading_order::partonic::Born["gq"]  = std::bind(&NoContribution, _1, _2);
	sigma::leading_order::partonic::Born["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::leading_order::partonic::Born["gqb"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["gg"]  = std::bind(&sigma::next_to_leading_order::partonic::hard::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::hard::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::hard::qbq, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["qg"]  = std::bind(&sigma::next_to_leading_order::partonic::hard::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["gq"]  = std::bind(&sigma::next_to_leading_order::partonic::hard::gq, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["qbg"] = std::bind(&sigma::next_to_leading_order::partonic::hard::qbg, _1, _2);
	sigma::next_to_leading_order::partonic::Hard["gqb"] = std::bind(&sigma::next_to_leading_order::partonic::hard::gqb, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["gg"]  = std::bind(&sigma::next_to_leading_order::partonic::soft::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::soft::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::soft::qbq, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["qg"]  = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["gq"]  = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Soft["gqb"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["gg"]  = std::bind(&sigma::next_to_leading_order::partonic::virt::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::virt::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::virt::qbq, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["qg"]  = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["gq"]  = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Virt["gqb"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["gg"]  = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["qg"]  = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["gq"]  = std::bind(&sigma::next_to_leading_order::partonic::coll_right_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["qbg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_left_z["gqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_right_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["gg"]  = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["qg"]  = std::bind(&sigma::next_to_leading_order::partonic::coll_right_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["gq"]  = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["qbg"] = std::bind(&sigma::next_to_leading_order::partonic::coll_right_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_right_z["gqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_left_z::qg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["gg"]  = std::bind(&sigma::next_to_leading_order::partonic::coll_1::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_1::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_1::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["qg"]  = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["gq"]  = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_1["gqb"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["gg"]  = std::bind(&sigma::next_to_leading_order::partonic::coll_0::gg, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["qqb"] = std::bind(&sigma::next_to_leading_order::partonic::coll_0::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["qbq"] = std::bind(&sigma::next_to_leading_order::partonic::coll_0::qqb, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["qg"]  = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["gq"]  = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["qbg"] = std::bind(&NoContribution, _1, _2);
	sigma::next_to_leading_order::partonic::Coll_0["gqb"] = std::bind(&NoContribution, _1, _2);
	std::map<std::string, Integrand> my_integrands;
	std::map<std::string, std::pair<double, double>> my_variables;
	my_variables["E1"]     = std::pair<double, double>{ 173.2, 13000.0 / 2.0 };
	my_variables["phi1"]   = std::pair<double, double>{ -M_PI, M_PI };
	my_variables["theta1"] = std::pair<double, double>{ 0.0, M_PI };
	my_variables["theta2"] = std::pair<double, double>{ 0.0, M_PI };
	my_variables["E3"] = std::pair<double, double>{ 0.0, 13000.0 };
	my_variables["phi3"] = std::pair<double, double>{ -M_PI, M_PI };
	my_variables["theta3"] = std::pair<double, double>{ 0.0, M_PI };
	const int Nf = NF - 1;
	// const int FDH = 0;  /* Four-Dimensional-Helicity scheme                    */
	const int HV_CDR = 1;  /* 't Hooft-Veltman scheme                           */
	double m = parameters.GetTopQuarkMass();	
	bsyppttinit_(&m, &Nf, &HV_CDR);

	my_integrands["mur=m"] = std::bind(&sigma::next_to_leading_order::Hadronic3, _1, _2, parameters);
	coupqcd_.gg[0] = -1.0, coupqcd_.gg[1] = -1.0, coupqcd_.g = 1.0; // 1.0 for gs 
	fermions_.fmass[10] = 173.2;

	unsigned int iseed = 3720758;
	const int lxlev = 1;
	rlxd_init(lxlev, iseed);

	Integral my_integral(my_variables, my_integrands);
	my_integral.ExecuteVegas(1, 10, 100000, 1);

	return 0;
}
