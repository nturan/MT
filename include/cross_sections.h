#pragma once

#include <map>
#include "amps.h"
#include "parameters.h"
#include "phase_space_generator.h"

typedef double Partonic(PhaseSpaceGenerator,Parameters&);
typedef std::map<std::string, std::function<double(PhaseSpaceGenerator, Parameters&)>> PartonicContribution;
namespace sigma::lo {
	double Hadronic(std::map<std::string, double> v, double* wgt, Parameters& p);
	namespace partonic {
		extern PartonicContribution Born;
		Partonic gg, qqb, qbq, NoContribution; 
	}
}

