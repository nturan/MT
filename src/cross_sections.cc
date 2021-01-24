#include "cross_sections.h"


std::map<std::string, Parameters*> parameter_sets;
std::map<std::string, std::vector<Histogram*>*> histogram_sets;

std::map<std::string, std::pair<double, double>> lo_variables, 
nlo_2_variables, nlo_3_variables, nlo_z_variables, nlo_j_variables;

std::map<std::string, Integrand> lo_integrands, lo_nlo_2_integrands,
								 nlo_2_integrands, nlo_3_integrands;

std::map<std::string, std::vector<Integral*>*> integrals; 

PartonicContribution lo::partonic::Born = {
	{"gg", &lo::partonic::gg},
	{"qqb", &lo::partonic::qqb},
	{"qbq", &lo::partonic::qbq},
	{"qg", &NoContribution},
	{"gq", &NoContribution},
	{"qbg", &NoContribution},
	{"gqb", &NoContribution}
};
PartonicContribution nlo::partonic::Hard = {
	{"gg",  &nlo::partonic::hard::gg},
	{"qqb", &nlo::partonic::hard::qqb},
	{"qbq", &nlo::partonic::hard::qbq},
	{"qg", &nlo::partonic::hard::qg},
	{"gq", &nlo::partonic::hard::gq},
	{"qbg", &nlo::partonic::hard::qbg},
	{"gqb", &nlo::partonic::hard::gqb}
};
PartonicContribution nlo::partonic::Soft = {
	{"gg",  &nlo::partonic::soft::gg},
	{"qqb", &nlo::partonic::soft::qqb},
	{"qbq", &nlo::partonic::soft::qbq},
	{"qg",  &NoContribution},
	{"gq",  &NoContribution},
	{"qbg", &NoContribution},
	{"gqb", &NoContribution}
};
PartonicContribution nlo::partonic::Virt = {
	{"gg",  &nlo::partonic::virt::gg},
	{"qqb", &nlo::partonic::virt::qqb},
	{"qbq", &nlo::partonic::virt::qbq},
	{"qg",  &NoContribution},
	{"gq",  &NoContribution},
	{"qbg", &NoContribution},
	{"gqb", &NoContribution}
};
PartonicContribution nlo::partonic::Coll_left_z = {
	{"gg",  &nlo::partonic::coll_left_z::gg},
	{"qqb", &nlo::partonic::coll_left_z::qqb},
	{"qbq", &nlo::partonic::coll_left_z::qqb},
	{"qg",  &nlo::partonic::coll_left_z::qg},
	{"gq",  &nlo::partonic::coll_right_z::qg},
	{"qbg", &nlo::partonic::coll_left_z::qg},
	{"gqb", &nlo::partonic::coll_right_z::qg}
};
PartonicContribution nlo::partonic::Coll_right_z = {
	{"gg",  &nlo::partonic::coll_left_z::gg},
	{"qqb", &nlo::partonic::coll_left_z::qqb},
	{"qbq", &nlo::partonic::coll_left_z::qqb},
	{"qg",  &nlo::partonic::coll_right_z::qg},
	{"gq",  &nlo::partonic::coll_left_z::qg},
	{"qbg", &nlo::partonic::coll_right_z::qg},
	{"gqb", &nlo::partonic::coll_left_z::qg}
};
PartonicContribution nlo::partonic::Coll_1 = {
	{"gg",  &nlo::partonic::coll_1::gg},
	{"qqb", &nlo::partonic::coll_1::qqb},
	{"qbq", &nlo::partonic::coll_1::qqb},
	{"qg",  &NoContribution},
	{"gq",  &NoContribution},
	{"qbg", &NoContribution},
	{"gqb", &NoContribution}
};
PartonicContribution nlo::partonic::Coll_0 = {
	{"gg",  &nlo::partonic::coll_0::gg},
	{"qqb", &nlo::partonic::coll_0::qqb},
	{"qbq", &nlo::partonic::coll_0::qqb},
	{"qg",  &NoContribution},
	{"gq",  &NoContribution},
	{"qbg", &NoContribution},
	{"gqb", &NoContribution}
};

void InitializeIntegrands(std::vector<std::string> histogram_strings, double ecms, double m) {
	
	lo_variables["k1p"]     = std::make_pair(0.0, std::sqrt( ecms*ecms/4.0 - m*m ));
	lo_variables["phi1"]   = std::make_pair(-M_PI, M_PI);
	lo_variables["theta1"] = std::make_pair(0.0, M_PI);
	lo_variables["theta2"] = std::make_pair(0.0, M_PI);

	nlo_2_variables["k1p"]     = std::make_pair(0.0, std::sqrt(ecms * ecms / 4.0 - m * m));
	nlo_2_variables["phi1"]   = std::make_pair(-M_PI, M_PI);
	nlo_2_variables["theta1"] = std::make_pair(0.0, M_PI);
	nlo_2_variables["theta2"] = std::make_pair(0.0, M_PI);
	nlo_2_variables["z"]      = std::make_pair(0.0, 1.0);

	nlo_3_variables["k1p"]     = std::make_pair(0.0, std::sqrt(ecms * ecms / 4.0 - m * m));
	nlo_3_variables["phi1"]   = std::make_pair(-M_PI, M_PI);
	nlo_3_variables["theta1"] = std::make_pair(0.0, M_PI);
	nlo_3_variables["theta2"] = std::make_pair(0.0, M_PI);
	nlo_3_variables["k3p"]     = std::make_pair(0.0, ecms);
	nlo_3_variables["phi3"]   = std::make_pair(-M_PI, M_PI);
	nlo_3_variables["theta3"] = std::make_pair(0.0, M_PI);


	for (auto it = parameter_sets.begin(); it != parameter_sets.end(); ++it) {

		std::vector<Histogram*>* hists = new std::vector<Histogram*>();
		for (auto ij = histogram_strings.begin(); ij != histogram_strings.end(); ++ij) {
			hists->push_back(new Histogram(*ij, it->second));
		}
		using namespace std::placeholders;
		histogram_sets.insert({ it->first, hists });
		lo_integrands[it->first] = std::bind(&lo::Hadronic, _1, _2, it->second, hists);
		lo_nlo_2_integrands[it->first] = std::bind(&nlo::Hadronic2WithBorn, _1, _2, it->second, hists);
		nlo_2_integrands[it->first] = std::bind(&nlo::Hadronic2, _1, _2, it->second, hists);
		nlo_3_integrands[it->first] = std::bind(&nlo::Hadronic3, _1, _2, it->second, hists);
	}

	integrals["lo"]  = new std::vector<Integral*>{ new Integral(lo_variables,    lo_integrands) };
	integrals["nlo"] = new std::vector<Integral*>{ new Integral(nlo_2_variables, nlo_2_integrands),
												   new Integral(nlo_3_variables, nlo_3_integrands) };
	integrals["lo+nlo"] = new std::vector<Integral*>{ new Integral(nlo_2_variables, lo_nlo_2_integrands),
											          new Integral(nlo_3_variables, nlo_3_integrands) };
}

void ExecuteIntegralsAndPrintResults(std::vector<std::string> perturbation_order, int iterations, int calls) {
	for (auto it = perturbation_order.begin(); it != perturbation_order.end(); ++it) {
		std::map<std::string, std::tuple<double, double, double>> results;
		for (auto ij = integrals.find(*it)->second->begin(); ij != integrals.find(*it)->second->end(); ++ij) {
			std::cout << "#INITIALISING GRID" << std::endl;
			(*ij)->ExecuteVegas(1, 30, 10000, 1);
			std::cout << "#VEGAS_START" << std::endl;
			(*ij)->ExecuteVegas(2, iterations, calls, 1);
			std::cout << "#VEGAS_END" << std::endl;

			for (auto ik = (*ij)->results_.begin(); ik != (*ij)->results_.end(); ++ik) {
				double val_new, err_new, chi_new;
				double val_old, err_old, chi_old;
				std::tie(val_new, err_new, chi_new) = ik->second;
				std::tie(val_old, err_old, chi_old) = results[ik->first];
				val_old = val_old + val_new;
				err_old = std::sqrt(err_old * err_old + err_new * err_new);
				chi_old = (chi_old != 0 ? 0.5 : 1.0) * (chi_old + chi_new);
				results[ik->first] = std::make_tuple(val_old, err_old, chi_old);
			}
		}
		for (auto ij = results.begin(); ij != results.end(); ++ij) {
			double val, err, chi;
			std::tie(val, err, chi) = ij->second;
			std::cout << "#RESULTS for " << *it << " with " << ij->first << ": "
				<< val << " +/- " << err
				<< " with chi2 = " << chi << std::endl;
			std::cout << "Histograms for parameter set: " << ij->first << std::endl;
			for (auto ik = histogram_sets.find(ij->first)->second->begin();
				ik != histogram_sets.find(ij->first)->second->end(); ++ik) {
				(*ik)->Print();
				(*ik)->Clear();
			}
		}
	}

}


double lo::Hadronic( std::map<std::string, double> var, double& wgt, Parameters *p, std::vector<Histogram*>* histograms) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var, p);
	if (PS.dGamma_ == 0) {
		return 0.0;
	}
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	for (auto it = p->channels_.begin(); it != p->channels_.end(); ++it) {
		std::string channel = *it;
		res += p->Fs[channel](x1, x2) * lo::partonic::Born[channel]( PS, p );
	}
	for (auto it = histograms->begin(); it != histograms->end(); ++it) {
		(*it)->Fill(PS, res * wgt);
	}
	return res;
}

double nlo::Hadronic2WithBorn(std::map<std::string, double> var, double& wgt, Parameters* p, std::vector<Histogram*>* histograms) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var, p);
	if (PS.dGamma_ == 0) {
		return 0.0;
	}
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double z = PS.z_;
	for (auto ij = p->channels_.begin(); ij != p->channels_.end(); ++ij) {
		res += p->Fs[*ij](x1, x2) * (
			  lo::partonic::Born[*ij](PS, p)
			+ nlo::partonic::Soft[*ij](PS, p)
			+ nlo::partonic::Virt[*ij](PS, p)
			+ nlo::partonic::Coll_0[*ij](PS, p)
			+ nlo::partonic::Coll_1[*ij](PS, p));
		res += p->Fs[*ij](x1 / z, x2) * nlo::partonic::Coll_left_z[*ij](PS, p)/z;
		res += p->Fs[*ij](x1, x2 / z) * nlo::partonic::Coll_right_z[*ij](PS, p)/z;
	}
	for (auto it = histograms->begin(); it != histograms->end(); ++it) {
		(*it)->Fill(PS, res * wgt);
	}
	return res;
}

double nlo::HadronicConstZ(std::map<std::string, double> var, double& wgt, Parameters* p, std::vector<Histogram*>* histograms) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var, p);
	if (PS.dGamma_ == 0) {
		return 0.0;
	}
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	for (auto ij = p->channels_.begin(); ij != p->channels_.end(); ++ij) {
		res += p->Fs[*ij](x1, x2) * ( nlo::partonic::Soft[*ij](PS, p)
								    + nlo::partonic::Virt[*ij](PS, p)
			                        + nlo::partonic::Coll_0[*ij](PS, p));
	}
	for (auto it = histograms->begin(); it != histograms->end(); ++it) {
		(*it)->Fill(PS, res * wgt);
	}
	return res;
}

double nlo::HadronicZ(IntegrationVariablesMap var, double& wgt, Parameters* p, std::vector<Histogram*>* histograms) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var, p);

	if (PS.dGamma_ == 0) {
		return 0.0;
	}
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double z = PS.z_;
	for (auto ij = p->channels_.begin(); ij != p->channels_.end(); ++ij) {
		res += p->Fs[*ij](x1, x2) * nlo::partonic::Coll_1[*ij](PS, p);
		res += p->Fs[*ij](x1 / z, x2) * nlo::partonic::Coll_left_z[*ij](PS, p) / z;
		res += p->Fs[*ij](x1, x2 / z) * nlo::partonic::Coll_right_z[*ij](PS, p) / z;
	}
	for (auto it = histograms->begin(); it != histograms->end(); ++it) {
		(*it)->Fill(PS, res * wgt);
	}
	return res;
}

double nlo::Hadronic2(std::map<std::string, double> var, double& wgt, Parameters* p, std::vector<Histogram*>* histograms) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var, p);
	if (PS.dGamma_ == 0) {
		return 0.0;
	}
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double z = PS.z_;
	for (auto ij = p->channels_.begin(); ij != p->channels_.end(); ++ij) {
		res += p->Fs[*ij](x1, x2) * (nlo::partonic::Soft[*ij](PS, p)
			+ nlo::partonic::Virt[*ij](PS, p)
			+ nlo::partonic::Coll_0[*ij](PS, p)
			+ nlo::partonic::Coll_1[*ij](PS, p));
		res += p->Fs[*ij](x1 / z, x2) * nlo::partonic::Coll_left_z[*ij](PS, p)/z;
		res += p->Fs[*ij](x1, x2 / z) * nlo::partonic::Coll_right_z[*ij](PS, p)/z;
	}
	for (auto it = histograms->begin(); it != histograms->end(); ++it) {
		(*it)->Fill(PS, res * wgt);
	}
	return res;
}

double nlo::Hadronic3(std::map<std::string, double> var, double& wgt, Parameters* p, std::vector<Histogram*>* histograms) {
	double res = 0.0;
	PhaseSpaceGenerator PS(var, p);
	if (PS.dGamma_ == 0) {
		return 0.0;
	}
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	for (auto ij = p->channels_.begin(); ij != p->channels_.end(); ++ij) {
		res += p->Fs[*ij](x1, x2) * nlo::partonic::Hard[*ij](PS, p);
	}
	for (auto it = histograms->begin(); it != histograms->end(); ++it) {
		(*it)->Fill(PS, res * wgt);
	}
	return res;
}
double lo::partonic::gg(PhaseSpaceGenerator PS, Parameters* p) {
	double s	  = PS.s_;
	double dGamma = PS.dGamma_;
	double coup   = p->GetSquaredGs();
	double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * sgg_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
	
	return result;
}

double lo::partonic::qqb(PhaseSpaceGenerator PS, Parameters* p) {
	double s      = PS.s_;
	double dGamma = PS.dGamma_;
	double coup   = p->GetSquaredGs();
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double lo::partonic::qbq(PhaseSpaceGenerator PS, Parameters * p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * suub_ttb_(PS.p2_, PS.p1_, PS.k1_, PS.k2_) * kPbarn;
	//std::cout << result << std::endl;
	return result;
}

double NoContribution(PhaseSpaceGenerator, Parameters * p) { return 0.0; }

double nlo::partonic::soft::gg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double dGamma = PS.dGamma_;
//	double coup = p.GetSquaredGs();
	double m = p->GetTopQuarkMass();
	double mur2 = p->GetSquaredRenormalizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
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

	Mtgg = 4.0 * M_PI * M_PI * std::pow(alpha_s, 2) * (Nc * Nc - 1.0)
		/ Nc / std::pow(1.0 - b2 * y2, 2) * Agg;

	Mgg = 2.0 * (Nc * Nc * (1.0 + b2 * y2) - 2.0) * Mtgg;

	double Phigg = 1.0 / 64.0;
	divTerm softGG = Phigg * alpha_s / M_PI * Ceps *
		(SA * Mgg + Mtgg * (SB + SC(y) + SC(-y)));

	return 1.0 / 2.0 / s * dGamma * softGG.eps0 * kPbarn;
}

double nlo::partonic::soft::qqb(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	double m = p->GetTopQuarkMass();
	double mur2 = p->GetSquaredRenormalizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
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

	divTerm result = bornQQ * alpha_s / M_PI * Ceps * (divTerm) {
		eps_2,
			eps_1,
			eps0, 0.0, 0.0};
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * result.eps0 * kPbarn;
}

double nlo::partonic::soft::qbq(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double x1 = PS.x1_;
	double x2 = PS.x2_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	double m = p->GetTopQuarkMass();
	double mur2 = p->GetSquaredRenormalizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
	double m2 = m * m;
	double beta = (x2 - x1) / (x1 + x2);
	double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
	double b = std::sqrt(1.0 - 4.0 * m2 / s);
	double y =-1.0 / b * (2.0 * PS.k1_[3] / std::sqrt(s) / gamma + beta);
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

	divTerm result = bornQQ * alpha_s / M_PI * Ceps * (divTerm) {
		eps_2,
			eps_1,
			eps0, 0.0, 0.0
	};
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * result.eps0 * kPbarn;
}

double nlo::partonic::virt::gg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double alpha_s = p->GetAlphaS();
	divTerm cGamma = { 0.0, 0.0, 1.0, 0.0, 0.0 };
	std::complex<double> amp[3];
	double mur2 = p->GetSquaredRenormalizationScale();
	double p1[4] = { PS.p1_[0], PS.p1_[1], PS.p1_[2], PS.p1_[3] };
	double p2[4] = { PS.p2_[0], PS.p2_[1], PS.p2_[2], PS.p2_[3] };
	double k1[4] = { PS.k1_[0], PS.k1_[1], PS.k1_[2], PS.k1_[3] };
	double k2[4] = { PS.k2_[0], PS.k2_[1], PS.k2_[2], PS.k2_[3] };
	bsyggttsq_(p1, p2, k1, k2, &mur2, amp);
	divTerm result = cGamma * 4.0 * M_PI * std::pow(alpha_s, 3)
		* (divTerm) { amp[0].real(), amp[1].real(), amp[2].real(), 0.0, 0.0 };
	return 1.0 / 2.0 / s * dGamma * result.eps0 * kPbarn;
}

double nlo::partonic::virt::qqb(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetAlphaS();
	divTerm cGamma = { 0.0, 0.0, std::pow(coup, 3), 0.0, 0.0 };
	std::complex<double> amp[3];
	double mur2 = p->GetSquaredRenormalizationScale();
	bsyqqttsq_(PS.p2_, PS.p1_, PS.k1_, PS.k2_, &mur2, amp);
	divTerm result = cGamma * 4.0 * M_PI
		* (divTerm) { amp[0].real(), amp[1].real(), amp[2].real(), 0.0, 0.0 };
	return 1.0 / 2.0 / s * dGamma * result.eps0 * kPbarn;
}

double nlo::partonic::virt::qbq(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetAlphaS();
	divTerm cGamma = { 0.0, 0.0, std::pow(coup, 3), 0.0, 0.0 };
	complex<double> amp[3];
	double mur2 = p->GetSquaredRenormalizationScale();
	bsyqqttsq_(PS.p1_, PS.p2_, PS.k1_, PS.k2_, &mur2, amp);
	divTerm result = cGamma * 4.0 * M_PI
		* (divTerm) { amp[0].real(), amp[1].real(), amp[2].real(), 0.0, 0.0 };
	return 1.0 / 2.0 / s * dGamma * result.eps0 * kPbarn;
}

double nlo::partonic::coll_left_z::gg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double z = PS.z_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	double mur2 = p->GetSquaredRenormalizationScale();
	double muf2 = p->GetSquaredFactorizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();

	double Feps = 2.0 * mur2 / s / xmin * z;
	double Ceps = mur2 / muf2;
	double z2 = z * z;
	double result = (-2 * Nc * std::pow(1 - z + z2, 2) * (ln(Ceps) - ln(Feps) + 2 * ln(1 - z))) / ((-1 + z) * z);
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * alpha_s / 2.0 / M_PI * result * sgg_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double nlo::partonic::coll_left_z::qqb(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double z = PS.z_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();

	double mur2 = p->GetSquaredRenormalizationScale();
	double muf2 = p->GetSquaredFactorizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
	double Ceps = mur2 / muf2;
	double Feps = 2.0 * mur2 / s / xmin * z;
	double z2 = z * z;
	double result = -(((1 - 2 * z + z2 + (1 + z2) * ln(Ceps) - (1 + z2) * ln(Feps) +
		2 * ln(1 - z) + 2 * z2 * ln(1 - z))) / (-1 + z));
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * alpha_s / 2.0 / M_PI * CF * result * suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double nlo::partonic::coll_left_z::qg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double z = PS.z_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();

	double mur2 = p->GetSquaredRenormalizationScale();
	double muf2 = p->GetSquaredFactorizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
	double gg = sgg_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_);
	divTerm Feps = { 0.0, 0.0, 1.0, ln(2.0 * mur2 / s / xmin / std::pow(1 - z,2) * z), 0.0 };
	divTerm result = -alpha_s / 2.0 / M_PI * Feps.qEps() *
		gg * CF * (divTerm) { 0.0, 0.0, (1.0 + std::pow(1.0 - z, 2)) / z, -z, 0.0 };


	divTerm FepsC = { 0.0, 0.0, 1.0, ln(mur2 / muf2), 0.0 };
	double Pgq = CF * (1.0 + std::pow(1 - z, 2)) / z;
	divTerm resultC = alpha_s / 2.0 / M_PI * FepsC.qEps() * gg * Pgq;

	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * (result + resultC).eps0 * kPbarn;
}

double nlo::partonic::coll_right_z::qg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double z = PS.z_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();

	double mur2 = p->GetSquaredRenormalizationScale();
	double muf2 = p->GetSquaredFactorizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
	double qq = suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_);

	divTerm Feps = { 0.0, 0.0, 1.0, ln(2.0 * mur2 / s / xmin / std::pow(1 - z,2) * z), 0.0 };
	divTerm result = -alpha_s / 2.0 / M_PI * Feps.qEps() *
		qq * 0.5 * (divTerm) { 0.0, 0.0, z* z + std::pow(1.0 - z, 2), z* z + std::pow(1.0 - z, 2) - 1.0, 0.0 };


	divTerm FepsC = { 0.0, 0.0, 1.0, ln(mur2 / muf2), 0.0 };
	double Pqg = 0.5 * (z * z + std::pow(1 - z, 2));
	divTerm resultC = alpha_s / 2.0 / M_PI * FepsC.qEps() * qq * Pqg;

	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * (result + resultC).eps0 * kPbarn;
}

double nlo::partonic::coll_1::gg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double z = PS.z_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();

	double mur2 = p->GetSquaredRenormalizationScale();
	double muf2 = p->GetSquaredFactorizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
	double Feps = 2.0 * mur2 / s / xmin;
	double Ceps = mur2 / muf2;
	double m2 = std::pow(p->GetTopQuarkMass(), 2);
	double x = 4.0 * m2 / s;
	double result = z > x ? 2 * Nc * (ln(Ceps) - ln(Feps) + 2 * ln(1 - z)) / (-1 + z) : 0.0;
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * alpha_s / 2.0 / M_PI * result * 2 * sgg_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double nlo::partonic::coll_1::qqb(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double z = PS.z_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();

	double mur2 = p->GetSquaredRenormalizationScale();
	double muf2 = p->GetSquaredFactorizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
	double m2 = std::pow(p->GetTopQuarkMass(), 2);
	double Ceps = mur2 / muf2;
	double Feps = 2.0 * mur2 / s / xmin;
	double x = 4.0 * m2 / s;
	double result = z > x ? (4 * (ln(Ceps) - ln(Feps) + 2 * ln(1 - z))) / (-1 + z) : 0.0;
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * alpha_s / 2.0 / M_PI * CF * result * suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double nlo::partonic::coll_0::gg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	double mur2 = p->GetSquaredRenormalizationScale();
	double muf2 = p->GetSquaredFactorizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
	double m2 = std::pow(p->GetTopQuarkMass(), 2);
	double Feps = 2.0 * mur2 / s / xmin;
	double Ceps = mur2 / muf2;

	double x = 4.0 * m2 / s;
	divTerm result = { 0.0, (2 + 11 * Nc - 2 * NF + 12 * Nc * ln(xmin)) / 6.,
	(ln(Ceps) * (2 + 11 * Nc - 2 * NF + 12 * Nc * ln(1 - x)) -
		 12 * Nc * (ln(Feps) * (ln(1 - x) - ln(xmin)) - ln2(1 - x) + ln2(xmin))) / 6., 0.0, 0.0 };
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * alpha_s / 2.0 / M_PI * result.eps0 * 2 * sgg_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double nlo::partonic::coll_0::qqb(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	double mur2 = p->GetSquaredRenormalizationScale();
	double muf2 = p->GetSquaredFactorizationScale();
	double xmin = p->GetCutParameter();
	double alpha_s = p->GetAlphaS();
	double m2 = std::pow(p->GetTopQuarkMass(), 2);
	double Ceps = mur2 / muf2;
	double Feps = 2.0 * mur2 / s / xmin;
	double x = 4.0 * m2 / s;
	divTerm result = { 0.0,  3.0 + 4.0 * ln(xmin),
	ln(Ceps) * (3 + 4 * ln(1 - x)) + 4 * (ln(Feps) * (-ln(1 - x) + ln(xmin)) + ln2(1 - x) - ln2(xmin)), 0.0, 0.0 };
	return 1.0 / 2.0 / s * dGamma * std::pow(coup, 2) * alpha_s/2.0 / M_PI * CF * result.eps0 * suub_ttb_(PS.p1_, PS.p2_, PS.k1_, PS.k2_) * kPbarn;
}

double nlo::partonic::hard::gg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	if (PS.not_soft_ && PS.not_collinear_) {
		double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 3) * sgg_ttbg_(PS.p1_, PS.p2_, PS.k1_, PS.k2_, PS.k3_) * kPbarn;
		if (std::isnan(result)) {
			return 0.0;
		}
		else {
			return result;
		}
	}
	else {
		return 0.0;
	}
}

double nlo::partonic::hard::qqb(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	if (PS.not_soft_ && PS.not_collinear_) {
		double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 3) * suub_ttbg_(PS.p1_, PS.p2_, PS.k1_, PS.k2_, PS.k3_) * kPbarn;
		if (std::isnan(result)) {
			return 0.0;
		}
		else {
			return result;
		}
	}
	else {
		return 0.0;
	}
}

double nlo::partonic::hard::qbq(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	if (PS.not_soft_ && PS.not_collinear_) {
		double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 3) * subu_ttbg_(PS.p1_, PS.p2_, PS.k1_, PS.k2_, PS.k3_) * kPbarn;
		if (std::isnan(result)) {
			return 0.0;
		}
		else {
			return result;
		}
	}
	else {
		return 0.0;
	}
}

double nlo::partonic::hard::qg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	if (PS.not_collinear_) {
		double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 3) * sug_ttbu_(PS.p1_, PS.p2_, PS.k1_, PS.k2_, PS.k3_) * kPbarn;
		if (std::isnan(result)) {
			return 0.0;
		}
		else {
			return result;
		}
	}
	else {
		return 0.0;
	}
}

double nlo::partonic::hard::qbg(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	if (PS.not_collinear_) {
		double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 3) * subg_ttbub_(PS.p1_, PS.p2_, PS.k1_, PS.k2_, PS.k3_) * kPbarn;
		if (std::isnan(result)) {
			return 0.0;
		}
		else {
			return result;
		}
	}
	else {
		return 0.0;
	}
}

double nlo::partonic::hard::gq(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	if (PS.not_collinear_) {
		double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 3) * sgu_ttbu_(PS.p1_, PS.p2_, PS.k1_, PS.k2_, PS.k3_) * kPbarn;
		if (std::isnan(result)) {
			return 0.0;
		}
		else {
			return result;
		}
	}
	else {
		return 0.0;
	}
}

double nlo::partonic::hard::gqb(PhaseSpaceGenerator PS, Parameters* p) {
	double s = PS.s_;
	double dGamma = PS.dGamma_;
	double coup = p->GetSquaredGs();
	if (PS.not_collinear_) {
		double result = 1.0 / 2.0 / s * dGamma * std::pow(coup, 3) * sgub_ttbub_(PS.p1_, PS.p2_, PS.k1_, PS.k2_, PS.k3_) * kPbarn;
		if (std::isnan(result)) {
			return 0.0;
		}
		else {
			return result;
		}
	}
	else {
		return 0.0;
	}
}