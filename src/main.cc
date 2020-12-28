// MT.cpp : Defines the entry point for the application.
//
#include <iostream>
#include "SysInfo.h"
#include "parton_distribution_function.h"
#include "cross_sections.h"


int main()
{
	SysInfo sys_info;
	PartonDistributionFunction* pdf = PartonDistributionFunction::GetInstance();
	
	pdf->InitPdfSet("CT10nlo");
	std::cout << "Hello World!\n" << pdf->Fs["gg"](0.4, 0.4, 173.2 * 173.2);
	return 0;
}
