// MT.cpp : Defines the entry point for the application.
//
#include <iostream>
#include "SysInfo.h"
#include "parton_distribution_function.h"


int main()
{
	SysInfo sys_info;
	PartonDistFunc* pdf = PartonDistFunc::GetInstance();
	pdf->InitPdfSet("CT10nlo");
	std::cout << "Hello World!\n";
	return 0;
}
