#include "cross_sections.h"

double sigma::lo::Hadronic(std::map<std::string, double> var, double* wgt, qcd::Parameters& param) {
	double res = 0.0;
	PartonDistributionFunction *pdf = PartonDistributionFunction::GetInstance();
	for (auto it = param.channels.begin(); it != param.channels.end(); ++it) {
		std::string channel = *it;
		res += pdf->Fs[channel](x1, x2, param.muf2_) * sigma::lo::partonic::Born[channel](var, wgt, param);
	}
}
