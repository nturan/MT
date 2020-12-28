#pragma once

#include "amps.h"
#include <map>
#include "qcd_parameters.h"
#include "parton_distribution_function.h"

namespace sigma::lo {
	double Hadronic(std::map<std::string, double> v, double *wgt, qcd::Parameters &param);
	namespace partonic {
		std::map<std::string,
			std::function<double(std::map<std::string, double> v,
				double* wgt,
				qcd::Parameters& param)>> Born;
	}
}