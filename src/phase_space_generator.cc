#include "phase_space_generator.h"


PhaseSpaceGenerator::PhaseSpaceGenerator(
	std::map<std::string, double> var, 
	Parameters* p)
{
	m_ = p->GetTopQuarkMass();
	ecms_ = p->GetColliderEnergy();
	xmin_ = p->GetCutParameter();
	if (var.size() == 4) {
		generatePS2(var);
	}
	else if (var.size() == 5) {
		generatePS2(var);
		z_ = var["z"];
	}
	else if (var.size() == 7) {
		if (var.count("x2") > 0) {
			generatePS3_CMS(var);
		}
		else {
			generatePS3(var);
		}
	}
	else {
		std::cout << "Incorrect amount of variables were given"
			      << " to the phase space generator. Aborting..." << std::endl;
		exit(1);
	}
}



void PhaseSpaceGenerator::generatePS2(std::map<std::string, double> var)
{
	double k1p = var["k1p"], 
		phi1 = var["phi1"], theta1 = var["theta1"], theta2 = var["theta2"];

	double eta1 = -std::log(std::tan(theta1 / 2.0));
	double eta2 = -std::log(std::tan(theta2 / 2.0));

	double dEta1 = 1.0 / std::sin(theta1);
	double dEta2 = 1.0 / std::sin(theta2);

	double m = m_;
	double m2 = std::pow(m, 2);
	double Ecms = ecms_;
	double Shad = std::pow(Ecms, 2);

//	double k1p = std::sqrt((E1 * E1 - m2)) / std::cosh(eta1);   
	double E1 = std::sqrt(std::pow(k1p * std::cosh(eta1), 2) + m2);

	double E2 = std::sqrt(m2 + std::pow(k1p * std::cosh(eta2), 2));

	if (E2 > Ecms / 2.0) { 
		dGamma_ = 0.0;
		return;
	}    /* bigger energies are kinematically impossible */
	if (E1 < m) {
		dGamma_ = 0.0;
		return;
	}
	if (E2 < m) {
		dGamma_ = 0.0;
		return;
	}


	double sinPhi = std::sin(phi1);
	double cosPhi = std::cos(phi1);

	/* the total jacobian                                                         */
	dGamma_ = 1.0 / 8.0 / M_PI / M_PI / E1 / E2
		* k1p * k1p * k1p * std::cosh(eta1) * std::cosh(eta2) / Shad
		* dEta1 * dEta2;



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

void PhaseSpaceGenerator::generatePS3(std::map<std::string, double> var) {

	double k1p = var["k1p"], phi1 = var["phi1"], theta1 = var["theta1"], 
		theta2 = var["theta2"];
	double k3p = var["k3p"], phi3 = var["phi3"], theta3 = var["theta3"];
	double eta1 = -std::log(std::tan(theta1 / 2.0));
	double eta2 = -std::log(std::tan(theta2 / 2.0));
	double eta3 = -std::log(std::tan(theta3 / 2.0));

	double dEta1 = 1.0 / std::sin(theta1);
	double dEta2 = 1.0 / std::sin(theta2);
	double dEta3 = 1.0 / std::sin(theta3);

	double m = m_;
	double m2 = std::pow(m, 2);
	double Ecms = ecms_;
	double Shad = std::pow(Ecms, 2);

	/* transversal momenta                                                        */
	//double k1p = std::sqrt((E1 * E1 - m2)) / std::cosh(eta1);
	double E1 = std::sqrt(std::pow(k1p * std::cosh(eta1), 2) + m2);
	//double k3p = E3 / std::cosh(eta3);
	double E3 = k3p * std::cosh(eta3);
	double k2p = 
		std::sqrt(k1p * k1p + k3p * k3p + 2.0 * k1p * k3p * std::cos(phi1 - phi3));

	/* azimuth angle for k2 is determined by k1 and k3. atan2 function delivers   */
	/* arctan in every quadrant.                                                  */
	double phi2 = std::atan2(k1p * std::sin(phi1) + k3p * std::sin(phi3),
		k1p * std::cos(phi1) + k3p * std::cos(phi3));

	double E2 = std::sqrt(std::pow(k2p * std::cosh(eta2), 2) + m2);
	if (E1 < m) {
		dGamma_ = 0.0;
		return;
	}
	if (E2 < m) {
		dGamma_ = 0.0;
		return;
	}

	double k1[4] = 
	{ E1, k1p * std::cos(phi1), k1p * std::sin(phi1), k1p * std::sinh(eta1) };
	double k2[4] = 
	{ E2,-k2p * std::cos(phi2),-k2p * std::sin(phi2), k2p * std::sinh(eta2) };
	double k3[4] = 
	{ E3, k3p * std::cos(phi3), k3p * std::sin(phi3), k3p * std::sinh(eta3) };


	/* total jacobian                                                             */
//	double jac = 1.0 / std::pow(2 * M_PI, 5) / 4.0 / Shad * k1p / std::cosh(eta1)
//		* beta * k3p / std::cosh(eta3) * dEta1 * dEta2 * dEta3;

	double jac = 
		1.0 / std::pow(2.0 * M_PI, 5) * 2.0 / Shad / 2.0 / E1 / 2.0 / E2 / 2.0 / E3
		* k1p * k1p * k2p * k3p * k3p
		* std::cosh(eta1) * std::cosh(eta2) * std::cosh(eta3)
		* dEta1 * dEta2 * dEta3;

	x1_ = (k1[0] + k2[0] + k3[0] + k1[3] + k2[3] + k3[3]) / Ecms;
	x2_ = (k1[0] + k2[0] + k3[0] - k1[3] - k2[3] - k3[3]) / Ecms;

	/* parton momenta fractions x1 and x2 must be in the interval 0 < xi < 1      */
	if (x1_ > 1.0) { dGamma_ = 0.0; return; }
	if (x1_ < 0.0) { dGamma_ = 0.0; return; }
	if (x2_ > 1.0) { dGamma_ = 0.0; return; }
	if (x2_ < 0.0) { dGamma_ = 0.0; return; }

	double p1[4] = { x1_ * Ecms / 2.0, 0.0, 0.0, x1_ * Ecms / 2.0 };
	double p2[4] = { x2_ * Ecms / 2.0, 0.0, 0.0,-x2_ * Ecms / 2.0 };



	s_ = x1_ * x2_ * Shad;
	if (s_ < 4.0 * m2) { dGamma_ = 0.0;  return; }

	p1_[0] = p1[0];  p1_[1] = p1[1]; p1_[2] = p1[2]; p1_[3] = p1[3];
	p2_[0] = p2[0];  p2_[1] = p2[1]; p2_[2] = p2[2]; p2_[3] = p2[3];

	k1_[0] = k1[0];  k1_[1] = k1[1]; k1_[2] = k1[2]; k1_[3] = k1[3];
	k2_[0] = k2[0];  k2_[1] = k2[1]; k2_[2] = k2[2]; k2_[3] = k2[3];
	k3_[0] = k3[0];  k3_[1] = k3[1]; k3_[2] = k3[2]; k3_[3] = k3[3];
	/* beta factor of the lorentz boost to CMS     */
	double b = (x1_ - x2_) / (x1_ + x2_);  
	double p1s[4], k3s[4];
	boostZ(p1_, b, p1s);                          /* boost in z direction only */
	boostZ(k3_, b, k3s);

	/* calculating slicing conditions                                             */
	double y = cos3p(p1s, k3s);
	double xg = 2.0 * k3s[0] / std::sqrt(s_);
	bool soft = xg < xmin_;
	bool coll = y > 1.0 - xmin_ || y < -(1.0 - xmin_);
	not_soft_ = !soft;
	not_collinear_ = !coll;
	dGamma_ = jac;
	/* not sure if it is necessary at this point. since E1 and E2 cannot          */
	/* go under m.                                                                */
	/******************************************************************************/
}


void PhaseSpaceGenerator::generatePS3_E3(std::map<std::string, double> var) {

	double k1p = var["k1p"], phi1 = var["phi1"], theta1 = var["theta1"],
		theta2 = var["theta2"];
	double E3 = std::exp(var["E3"]), phi3 = var["phi3"], theta3 = var["theta3"];
	double eta1 = -std::log(std::tan(theta1 / 2.0));
	double eta2 = -std::log(std::tan(theta2 / 2.0));
	double eta3 = -std::log(std::tan(theta3 / 2.0));

	double dEta1 = 1.0 / std::sin(theta1);
	double dEta2 = 1.0 / std::sin(theta2);
	double dEta3 = 1.0 / std::sin(theta3);

	double m = m_;
	double m2 = std::pow(m, 2);
	double Ecms = ecms_;
	double Shad = std::pow(Ecms, 2);

	/* transversal momenta                                                        */
	//double k1p = std::sqrt((E1 * E1 - m2)) / std::cosh(eta1);
	double E1 = std::sqrt(std::pow(k1p * std::cosh(eta1), 2) + m2);
	double k3p = E3 / std::cosh(eta3);
//	double E3 = k3p * std::cosh(eta3);
	double k2p = 
		std::sqrt(k1p * k1p + k3p * k3p + 2.0 * k1p * k3p * std::cos(phi1 - phi3));

	/* azimuth angle for k2 is determined by k1 and k3. atan2 function delivers   */
	/* arctan in every quadrant.                                                  */
	double phi2 = std::atan2(k1p * std::sin(phi1) + k3p * std::sin(phi3),
		k1p * std::cos(phi1) + k3p * std::cos(phi3));

	double E2 = std::sqrt(std::pow(k2p * std::cosh(eta2), 2) + m2);
	if (E1 < m) {
		dGamma_ = 0.0;
		return;
	}
	if (E2 < m) {
		dGamma_ = 0.0;
		return;
	}

	double k1[4] = 
	{ E1, k1p * std::cos(phi1), k1p * std::sin(phi1), k1p * std::sinh(eta1) };
	double k2[4] = 
	{ E2,-k2p * std::cos(phi2),-k2p * std::sin(phi2), k2p * std::sinh(eta2) };
	double k3[4] = 
	{ E3, k3p * std::cos(phi3), k3p * std::sin(phi3), k3p * std::sinh(eta3) };

	/* total jacobian                                                             */
//	double jac = 1.0 / std::pow(2 * M_PI, 5) / 4.0 / Shad * k1p / std::cosh(eta1)
//		* beta * k3p / std::cosh(eta3) * dEta1 * dEta2 * dEta3;

	double jac = 
		1.0 / std::pow(2.0 * M_PI, 5) * 2.0 / Shad / 2.0 / E1 / 2.0 / E2 / 2.0 / E3
		* k1p * k1p * k2p * k3p * k3p * E3
		* std::cosh(eta1) * std::cosh(eta2)
		* dEta1 * dEta2 * dEta3;

	x1_ = (k1[0] + k2[0] + k3[0] + k1[3] + k2[3] + k3[3]) / Ecms;
	x2_ = (k1[0] + k2[0] + k3[0] - k1[3] - k2[3] - k3[3]) / Ecms;

	/* parton momenta fractions x1 and x2 must be in the interval 0 < xi < 1      */
	if (x1_ > 1.0) { dGamma_ = 0.0; return; }
	if (x1_ < 0.0) { dGamma_ = 0.0; return; }
	if (x2_ > 1.0) { dGamma_ = 0.0; return; }
	if (x2_ < 0.0) { dGamma_ = 0.0; return; }

	double p1[4] = { x1_ * Ecms / 2.0, 0.0, 0.0, x1_ * Ecms / 2.0 };
	double p2[4] = { x2_ * Ecms / 2.0, 0.0, 0.0,-x2_ * Ecms / 2.0 };



	s_ = x1_ * x2_ * Shad;
	if (s_ < 4.0 * m2) { dGamma_ = 0.0;  return; }

	p1_[0] = p1[0];  p1_[1] = p1[1]; p1_[2] = p1[2]; p1_[3] = p1[3];
	p2_[0] = p2[0];  p2_[1] = p2[1]; p2_[2] = p2[2]; p2_[3] = p2[3];

	k1_[0] = k1[0];  k1_[1] = k1[1]; k1_[2] = k1[2]; k1_[3] = k1[3];
	k2_[0] = k2[0];  k2_[1] = k2[1]; k2_[2] = k2[2]; k2_[3] = k2[3];
	k3_[0] = k3[0];  k3_[1] = k3[1]; k3_[2] = k3[2]; k3_[3] = k3[3];

	/* beta factor of the lorentz boost to CMS     */
	double b = (x1_ - x2_) / (x1_ + x2_);  
	double p1s[4], k3s[4];
	boostZ(p1_, b, p1s);                          /* boost in z direction only */
	boostZ(k3_, b, k3s);

	/* calculating slicing conditions                                             */
	double y = cos3p(p1s, k3s);
	double xg = 2.0 * k3s[0] / std::sqrt(s_);
	bool soft = xg < xmin_;
	bool coll = y > 1.0 - xmin_ || y < -(1.0 - xmin_);
	not_soft_ = !soft;
	not_collinear_ = !coll;
	dGamma_ = jac;
	/* not sure if it is necessary at this point. since E1 and E2 cannot          */
	/* go under m.                                                                */
	/******************************************************************************/
}

void PhaseSpaceGenerator::generatePS3_CMS(std::map<std::string, double> var) {
	double s = std::exp(var["lns"]), x2 = var["x2"], xE3 = var["xE3"];

	double sqrt_s_2 = std::sqrt(s) / 2.0;
	double dE3 = std::log(sqrt_s_2 - m_ * m_ / sqrt_s_2) + 35.0;
	double E3 = std::exp(dE3 * xE3 - 35.0);
	double costheta3 = var["costheta3"], phi3 = var["phi3"];
	double costheta12 = var["costheta12"], phi12 = var["phi12"];
	double dx2 = 1 - s / ecms_ / ecms_;
	x2 = x2*dx2 + s / ecms_ / ecms_;
	double x1 = s / x2 / ecms_ / ecms_;
	double s12 = s - 4 * E3 * sqrt_s_2;
	//if (x1 > 1.0 || x1 < 0) { dGamma_ = 0.0; return; }
	/*
	std::cout << "x1: " << x1 << "\n"
		<< "x2: " << x2 << "\n"
		<< "s12: " << s12 << "\n"
		<< "costheta3: " << costheta3 << "\n"
		<< "costheta12: " << costheta12 << "\n"
		<< "phi3: " << phi3 << "\n"
		<< "phi12: " << phi12 << std::endl; */


	double sintheta3  = std::sqrt( (1.0 - costheta3 ) * (1.0 + costheta3 ) );
	double sintheta12 = std::sqrt( (1.0 - costheta12) * (1.0 + costheta12) );

	/* I need to control s12 in any other way */
	//if (s < 4.0 * m_ * m_) { dGamma_ = 0.0; return; }
	//if (s12 < 4.0 * m_ * m_) { dGamma_ = 0.0; return; }
	//if (s12 > s) { dGamma_ = 0.0; return; }
	double beta3 = 1.0 - s12 / s;
	double beta12 = std::sqrt( 1.0 - 4.0*m_*m_/s12 );

	double e3[3] = { sintheta3 * std::cos(phi3), 
					 sintheta3 * std::sin(phi3), 
					 costheta3 };

	double k12[4] = {sqrt_s_2 * (1.0 + s12 / s),
				     sqrt_s_2 * beta3 * e3[0],
	                 sqrt_s_2 * beta3 * e3[1],
	                 sqrt_s_2 * beta3 * e3[2] };

	k3_[0] =   sqrt_s_2 * (1.0 - s12 / s);
	k3_[1] = - sqrt_s_2 * beta3 * e3[0];
	k3_[2] = - sqrt_s_2 * beta3 * e3[1];
	k3_[3] = - sqrt_s_2 * beta3 * e3[2];

	double sqrt_s12_2 = std::sqrt(s12) / 2.0;
	double e12[3] = { sintheta12 * std::cos(phi12),
		              sintheta12 * std::sin(phi12),
		              costheta12};

	double k1[4] = { sqrt_s12_2,
		             sqrt_s12_2 * beta12 * e12[0],
		             sqrt_s12_2 * beta12 * e12[1],
		             sqrt_s12_2 * beta12 * e12[2]};

	double k2[4] = { sqrt_s12_2,
	               - sqrt_s12_2 * beta12 * e12[0],
	               - sqrt_s12_2 * beta12 * e12[1],
	               - sqrt_s12_2 * beta12 * e12[2]};

	double gamma   = k12[0] / std::sqrt(s12);
	double beta[3] = { -k12[1] / k12[0], -k12[2] / k12[0],- k12[3] / k12[0] };
	double boost_vector[4] = {gamma, beta[0], beta[1], beta[2]};

	BoostInGeneratlDirection(k1, boost_vector, k1_);
	BoostInGeneratlDirection(k2, boost_vector, k2_);

	p1_[0] = sqrt_s_2; p1_[1] = 0.0; p1_[2] = 0.0; p1_[3] = sqrt_s_2;
	p2_[0] = sqrt_s_2; p2_[1] = 0.0; p2_[2] = 0.0; p2_[3] = -sqrt_s_2;
	//k1_[0] = k1[0]; k1_[1] = k1[1]; k1_[2] = k1[2]; k1_[3] = k1[3];
	//k2_[0] = k2[0]; k2_[1] = k2[1]; k2_[2] = k2[2]; k2_[3] = k2[3];
	/*
	std::cout << "beta12: " << beta12 << std::endl;
	std::cout << "beta3: " << beta3 << std::endl;
	std::cout << "k1: " << k1_[0] << " " << k1_[1] << " " << k1_[2] << " " << k1_[3] 
		<< " with " << calcInvariantMass(k1_) << "GeV mass "<< std::endl;
	std::cout << "k2: " << k2_[0] << " " << k2_[1] << " " << k2_[2] << " " << k2_[3] 
		<< " with " << calcInvariantMass(k2_) << "GeV mass " << std::endl; 
	std::cout << "k3: " << k3_[0] << " " << k3_[1] << " " << k3_[2] << " " << k3_[3] 
		<< " with " << calcInvariantMass(k3_) << "GeV mass " << std::endl; 
	
	std::cout << "p1 + p2 - k1 - k2 - k3: " 
		<< p1_[0] + p2_[0] - k1_[0] - k2_[0] - k3_[0] << " "
		<< p1_[1] + p2_[1] - k1_[1] - k2_[1] - k3_[1] << " "
		<< p1_[2] + p2_[2] - k1_[2] - k2_[2] - k3_[2] << " "
		<< p1_[3] + p2_[3] - k1_[3] - k2_[3] - k3_[3] << std::endl;

	*/
	
	//double y = cos3p(p1_, k3_);
	s_ = s;
	x1_ = x1;
	x2_ = x2;
	double xg = k3_[0] / sqrt_s_2;
	bool soft = xg < xmin_;
	//bool coll = y > 1.0 - xmin_ || y < -(1.0 - xmin_);
	not_soft_ = !soft;
	not_collinear_ = true;
	double jac = 1.0 / 64.0 / std::pow(2 * M_PI, 5) * beta12 * beta3;
	dGamma_ = jac*s*x1/sqrt_s_2*dx2*dE3*E3;
}