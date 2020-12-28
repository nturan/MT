#include "cross_sections.h"

double sigma::lo::Hadronic( std::map<std::string, double> var, double* wgt ) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var);
	PartonDistributionFunction *pdf = PartonDistributionFunction::GetInstance();
	for (auto it = param.channels.begin(); it != param.channels.end(); ++it) {
		std::string channel = *it;
		res += pdf->Fs[channel](PS.x1_, PS.x2_, param.muf2_)
			* sigma::lo::partonic::Born[channel]( PS );
	}
}
