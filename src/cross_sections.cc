#include "cross_sections.h"

PartonicContribution sigma::lo::partonic::Born;


double sigma::lo::Hadronic( std::map<std::string, double> var, double* wgt, Parameters &p) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var, p);
	if (PS.dGamma_ == 0) {
		return 0.0;
	}
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	for (auto it = p.channels.begin(); it != p.channels.end(); ++it) {
		std::string channel = *it;
		res += p.Fs[channel](x1, x2) * sigma::lo::partonic::Born[channel]( PS, p );
	}
	return res;
}

double sigma::lo::partonic::gg(PhaseSpaceGenerator PS, Parameters & p) {
	double s	  = PS.s_;
	double dGamma = PS.dGamma_;
	double coup   = PS.p_->GetSquaredGs();
	double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * sgg_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
	
	return result;
}

double sigma::lo::partonic::qqb(PhaseSpaceGenerator PS, Parameters & p) {
	double s      = PS.s_;
	double dGamma = PS.dGamma_;
	double coup   = PS.p_->GetSquaredGs();
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double sigma::lo::partonic::qbq(PhaseSpaceGenerator PS, Parameters & p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = PS.p_->GetSquaredGs();
	double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * suub_ttb_(PS.p2_, PS.p1_, PS.k1_, PS.k2_) * kPbarn;
	//std::cout << result << std::endl;
	return result;
}

double sigma::lo::partonic::NoContribution(PhaseSpaceGenerator, Parameters & p) { return 0.0; }