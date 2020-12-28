#include "qcd_parameters.h"

qcd::Parameters* qcd::Parameters::GetInstance() {
	if (has_instance) {
		return parameters;
	}
	else {
		parameters = new Parameters();
		has_instance = true;
	}
}