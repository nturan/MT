// MT.cpp : Defines the entry point for the application.
//
#include <iostream>
#include "SysInfo.h"
#include "parton_distribution_function.h"
#include "cross_sections.h"
#include "qcd_parameters.h"


int main()
{
	SysInfo sys_info;
	PartonDistributionFunction* pdf = PartonDistributionFunction::GetInstance();
	Parameters* parameters = Parameters::GetInstance();
	parameters->SetColliderEnergy(13000);
	parameters->SetTopQuarkMass(173.2);
	parameters->SetFactorizationScale(173.2);
	parameters->SetRenormalizationScale(173.2);
	parameters->channels.push_back("gg" );
	parameters->channels.push_back("qqb");
	parameters->channels.push_back("qbq");
	parameters->channels.push_back("qg" );
	parameters->channels.push_back("gq" );
	parameters->channels.push_back("qbg");
	parameters->channels.push_back("gqb");

	
	pdf->InitPdfSet("CT10nlo");
	return 0;
}
