#pragma once

#include "amps.h"
#include <map>
#include "parameters.h"
#include "phase_space_generator.h"

typedef double Partonic(PhaseSpaceGenerator);

namespace sigma::lo {
	double Hadronic(std::map<std::string, double> v, double *wgt, Parameters &p);
	namespace partonic {
		std::map<std::string, std::function<double(PhaseSpaceGenerator)> > Born;
		Partonic gg, qqb, qbq;
	}
}



