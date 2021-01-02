#include "phase_space_generator.h"


PhaseSpaceGenerator::PhaseSpaceGenerator(std::map<std::string, double> var, Parameters& p)
{
	m_ = p.GetTopQuarkMass();
	ecms_ = p.GetColliderEnergy();
	if (var.size() == 4) {
		generatePS2(var);
	}
	else if (var.size() == 5) {
		generatePS2(var);
		z_ = var["z"];
	}
	else {
		std::cout << "Incorrect amount of variables were given"
			      << " to the phase space generator. Aborting..." << std::endl;
		exit(1);
	}
}



void PhaseSpaceGenerator::generatePS2(std::map<std::string, double> var)
{
	double E1 = var["E1"], phi1 = var["phi1"], theta1 = var["theta1"], theta2 = var["theta2"];

	double eta1 = -std::log(std::tan(theta1 / 2.0));
	double eta2 = -std::log(std::tan(theta2 / 2.0));

	double dEta1 = 1.0 / std::sin(theta1);
	double dEta2 = 1.0 / std::sin(theta2);

	double m = m_;
	double m2 = std::pow(m, 2);
	double Ecms = ecms_;
	double Shad = std::pow(Ecms, 2);

	double k1p = std::sqrt((E1 * E1 - m2)) / std::cosh(eta1);   /* transversal momenta      */

	double E2 = std::sqrt(m2 + std::pow(k1p * std::cosh(eta2), 2));

	if (E2 > Ecms / 2.0) { 
		dGamma_ = 0.0;
		return;
	}    /* bigger energies are kinematically impossible */
	if (E2 < m) {
		dGamma_ = 0.0;
		return;
	}

	double beta = std::sqrt(1.0 - m2 / E2 / E2);

	double sinPhi = std::sin(phi1);
	double cosPhi = std::cos(phi1);

	/* the total jacobian                                                         */
	dGamma_ = 1.0 / 8.0 / M_PI / M_PI * beta * k1p / std::cosh(eta1) / Shad * dEta1 * dEta2;



	x1_ = (E1 + E2 + k1p * sinh(eta1) + k1p * sinh(eta2)) / Ecms;
	x2_ = (E1 + E2 - k1p * sinh(eta1) - k1p * sinh(eta2)) / Ecms;
	/* parton momenta fractions x1 and x2 must be in the interval 0 < xi < 1      */
	if (x1_ > 1.0) { dGamma_ = 0.0; return; }
	if (x1_ < 0.0) { dGamma_ = 0.0; return; }
	if (x2_ > 1.0) { dGamma_ = 0.0; return; }
	if (x2_ < 0.0) { dGamma_ = 0.0; return; }




	s_ = x1_ * x2_ * Shad;
	/* not sure if it is necessary at this point. since E1 and E2 cannot          */
	/* go under m.                                                                */
	if (s_ < 4.0 * m2) { dGamma_ = 0.0;  return; }

	p1_[0] = x1_ * Ecms / 2.0; p2_[0] = x2_ * Ecms / 2.0;
	p1_[1] = 0.0;              p2_[1] = 0.0;
	p1_[2] = 0.0;              p2_[2] = 0.0;
	p1_[3] = x1_ * Ecms / 2.0; p2_[3] = -x2_ * Ecms / 2.0;

	/*     top quark                    anti top quark                            */
	k1_[0] = E1;                k2_[0] = E2;
	k1_[1] = k1p * cosPhi;      k2_[1] = -k1p * cosPhi;
	k1_[2] = k1p * sinPhi;      k2_[2] = -k1p * sinPhi;
	k1_[3] = k1p * sinh(eta1);  k2_[3] = k1p * sinh(eta2);
	/******************************************************************************/
}