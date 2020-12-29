
#include <iostream>
#include "SysInfo.h"
#include "cross_sections.h"
#include "parameters.h"


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
	std::cout << parameters.Fs["gg"](0.3, 0.5) << std::endl;
	

	return 0;
}
