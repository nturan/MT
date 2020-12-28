#include "cross_sections.h"

double sigma::lo::Hadronic( std::map<std::string, double> var, double* wgt, double scale_factor ) {
	double res = 0.0;
	PartonDistributionFunction* pdf = PartonDistributionFunction::GetInstance();
	Parameters* parameter = Parameters::GetInstance();
	PhaseSpaceGenerator PS(var);
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double muf2 = std::pow(parameter->GetFactorizationScale(), 2);
	for (auto it = parameter->channels.begin(); it != parameter->channels.end(); ++it) {
		std::string channel = *it;
		res += pdf->Fs[channel](x1, x2, muf2) * sigma::lo::partonic::Born[channel]( PS );
	}
}

double sigma::lo::partonic::gg(PhaseSpaceGenerator PS) {
	Parameters* parameter = Parameters::GetInstance();
	double s	  = PS.s_;
	double dGamma = PS.dGamma_;
	double coup   = parameter->GetAlphaS();
	return 0.5 / s * dGamma * std::pow(coup, 2) * sgg_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_);
}
