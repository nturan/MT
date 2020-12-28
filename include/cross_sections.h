#pragma once

#include "amps.h"
#include <map>
#include "qcd_parameters.h"
#include "parton_distribution_function.h"
#include "phase_space_generator.h"

namespace sigma::lo {
	double Hadronic(std::map<std::string, double> v, double *wgt);
	namespace partonic {
		std::map<std::string,
			std::function<double( PhaseSpaceGenerator PS )>> Born;
	}
}

