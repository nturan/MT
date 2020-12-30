#include "cross_sections.h"

PartonicContribution sigma::leading_order::partonic::Born;
PartonicContribution sigma::next_to_leading_order::partonic::Soft;
PartonicContribution sigma::next_to_leading_order::partonic::Virt;
PartonicContribution sigma::next_to_leading_order::partonic::Coll_left_z;
PartonicContribution sigma::next_to_leading_order::partonic::Coll_right_z;
PartonicContribution sigma::next_to_leading_order::partonic::Coll_1;
PartonicContribution sigma::next_to_leading_order::partonic::Coll_0;


double sigma::leading_order::Hadronic( std::map<std::string, double> var, double* wgt, Parameters &p) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var, p);
	if (PS.dGamma_ == 0) {
		return 0.0;
	}
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	for (auto it = p.channels.begin(); it != p.channels.end(); ++it) {
		std::string channel = *it;
		res += p.Fs[channel](x1, x2) * sigma::leading_order::partonic::Born[channel]( PS, p );
	}
	return res;
}

double sigma::leading_order::partonic::gg(PhaseSpaceGenerator PS, Parameters & p) {
	double s	  = PS.s_;
	double dGamma = PS.dGamma_;
	double coup   = PS.p_->GetSquaredGs();
	double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * sgg_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
	
	return result;
}

double sigma::leading_order::partonic::qqb(PhaseSpaceGenerator PS, Parameters & p) {
	double s      = PS.s_;
	double dGamma = PS.dGamma_;
	double coup   = PS.p_->GetSquaredGs();
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double sigma::leading_order::partonic::qbq(PhaseSpaceGenerator PS, Parameters & p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = PS.p_->GetSquaredGs();
	double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * suub_ttb_(PS.p2_, PS.p1_, PS.k1_, PS.k2_) * kPbarn;
	//std::cout << result << std::endl;
	return result;
}

double NoContribution(PhaseSpaceGenerator, Parameters & p) { return 0.0; }

double sigma::next_to_leading_order::partonic::soft::gg(PhaseSpaceGenerator PS, Parameters& p) {
	double s = PS.s_;
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double m = p.GetTopQuarkMass();
	double mur2 = p.GetSquaredRenormalizationScale();
	double xmin = p.GetCutParameter();
	double m2 = m * m;
	double beta = (x2 - x1) / (x1 + x2);
	double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
	double b = std::sqrt(1.0 - 4.0 * m2 / s), b2 = b * b;
	double y = 1.0 / b * (2.0 * PS.k1_[3] / std::sqrt(s) / gamma + beta), y2 = y * y;
	double x = (1.0 - b) / (1.0 + b);

	divTerm Ceps = { 0.0, 0.0, 1.0,
							  ln(mur2 / s / xmin / xmin),
							  0.5 * ln2(mur2 / s / xmin / xmin) };


	divTerm SA = { Nc,
				  CF,
				  -Nc / 6.0 * M_PI * M_PI - CF / b * ln(x), 0.0, 0.0 };

	divTerm SB = { 0.0,

  (1.0 + b2) / 2.0 / b / Nc * (Nc * Nc * (1.0 - b2 * y2) + 2.0) * ln(x),

  (1.0 + b2) / 2.0 / b / Nc * (Nc * Nc * (1.0 - b2 * y2) + 2.0) *
							  (-2.0 * Li2(1.0 - x) - 0.5 * ln2(x))
		  + Nc / 2.0 * (4.0 - Nc * Nc * (1.0 + b2 * y2)) * ln2(x), 0.0, 0.0 };

	function<double(double)> S = [b](double y) {
		return 2.0 * Li2(-b * (1.0 - y) / (1.0 - b))
			- 2.0 * Li2(-b * (1.0 + y) / (1.0 - b * y))
			+ ln2((1.0 - b * y) / (1.0 - b));
	};

	function<divTerm(double)> SC = [s, b, b2, S](double y) {
		return (divTerm) {
			0.0,
				Nc / 2.0 * (4.0 - Nc * Nc * std::pow(1.0 + b * y, 2))
				* ln(std::pow(1.0 - b * y, 2) / (1.0 - b2)),
				-Nc / 2.0 * (4.0 - Nc * Nc * std::pow(1.0 + b * y, 2)) * S(y), 0.0, 0.0
		};
	};


	double Mgg = 0.0, Mtgg = 0.0, Agg = 0.0;

	Agg =
		1.0 + 2.0 * b2 * (1.0 - y2) - std::pow(b, 4) * (1.0 + std::pow(1.0 - y2, 2));

	Mtgg = 4.0 * M_PI * M_PI * (Nc * Nc - 1.0)
		/ Nc / std::pow(1.0 - b2 * y2, 2) * Agg;

	Mgg = 2.0 * (Nc * Nc * (1.0 + b2 * y2) - 2.0) * Mtgg;

	double Phigg = 1.0 / 64.0;
	divTerm softGG = Phigg * 1.0 / M_PI * Ceps *
		(SA * Mgg + Mtgg * (SB + SC(y) + SC(-y)));

	return softGG.eps0;
}

double sigma::next_to_leading_order::partonic::soft::qqb(PhaseSpaceGenerator PS, Parameters& p) {
	double s = PS.s_;
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double m = p.GetTopQuarkMass();
	double mur2 = p.GetSquaredRenormalizationScale();
	double xmin = p.GetCutParameter();
	double m2 = m * m;
	double beta = (x2 - x1) / (x1 + x2);
	double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
	double b = std::sqrt(1.0 - 4.0 * m2 / s);
	double y = 1.0 / b * (2.0 * PS.k1_[3] / std::sqrt(s) / gamma + beta);
	double x = (1.0 - b) / (1.0 + b);

	function<double(double)> S = [b](double y) {
		return 2.0 * Li2(-b * (1.0 - y) / (1.0 - b))
			- 2.0 * Li2(-b * (1.0 + y) / (1.0 - b * y))
			+ ln2((1.0 - b * y) / (1.0 - b));
	};

	double eps_2 = CF;
	double eps_1 = CF
		- (1.0 + b * b) / 4.0 / b / Nc * ln(x)
		- 2.0 / Nc * ln((1.0 + b * y) / (1.0 - b * y))
		- Nc / 2.0 * ln((1.0 - b * y) * (1.0 - b * y) / (1.0 - b * b));
	double eps0 = -CF * M_PI * M_PI / 6.0
		+ (1.0 + b * b) / 2.0 / Nc / b * (Li2(1.0 - x) + 0.25 * ln2(x))
		- CF / b * ln(x)
		- Nc / 4.0 * ln2(x)
		+ (Nc * Nc - 2.0) / 2.0 / Nc * S(y)
		+ 1.0 / Nc * S(-y);

	divTerm Ceps = { 0.0, 0.0, 1.0,
								ln(mur2 / s / xmin / xmin),
								0.5 * ln2(mur2 / s / xmin / xmin) };
	double bornQQ = suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_);

	divTerm result = bornQQ * 1.0 / M_PI * Ceps * (divTerm) {
		eps_2,
			eps_1,
			eps0, 0.0, 0.0};
	return result.eps0;
}

double sigma::next_to_leading_order::partonic::soft::qbq(PhaseSpaceGenerator PS, Parameters& p) {
	double s = PS.s_;
	double x1 = PS.x2_;
	double x2 = PS.x1_;
	double m = p.GetTopQuarkMass();
	double mur2 = p.GetSquaredRenormalizationScale();
	double xmin = p.GetCutParameter();
	double m2 = m * m;
	double beta = (x2 - x1) / (x1 + x2);
	double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
	double b = std::sqrt(1.0 - 4.0 * m2 / s);
	double y = 1.0 / b * (2.0 * PS.k1_[3] / std::sqrt(s) / gamma + beta);
	double x = (1.0 - b) / (1.0 + b);

	function<double(double)> S = [b](double y) {
		return 2.0 * Li2(-b * (1.0 - y) / (1.0 - b))
			- 2.0 * Li2(-b * (1.0 + y) / (1.0 - b * y))
			+ ln2((1.0 - b * y) / (1.0 - b));
	};

	double eps_2 = CF;
	double eps_1 = CF
		- (1.0 + b * b) / 4.0 / b / Nc * ln(x)
		- 2.0 / Nc * ln((1.0 + b * y) / (1.0 - b * y))
		- Nc / 2.0 * ln((1.0 - b * y) * (1.0 - b * y) / (1.0 - b * b));
	double eps0 = -CF * M_PI * M_PI / 6.0
		+ (1.0 + b * b) / 2.0 / Nc / b * (Li2(1.0 - x) + 0.25 * ln2(x))
		- CF / b * ln(x)
		- Nc / 4.0 * ln2(x)
		+ (Nc * Nc - 2.0) / 2.0 / Nc * S(y)
		+ 1.0 / Nc * S(-y);

	divTerm Ceps = { 0.0, 0.0, 1.0,
								ln(mur2 / s / xmin / xmin),
								0.5 * ln2(mur2 / s / xmin / xmin) };
	double bornQQ = suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_);

	divTerm result = bornQQ * 1.0 / M_PI * Ceps * (divTerm) {
		eps_2,
			eps_1,
			eps0, 0.0, 0.0
	};
	return result.eps0;
}