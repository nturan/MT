#include "qcd_parameters.h"

Parameters* Parameters::GetInstance() {
	if (has_instance) {
		return parameters;
	}
	else {
		parameters = new Parameters();
		return parameters;
	}
}

Parameters::Parameters(){ has_instance = true; }
Parameters::~Parameters(){ has_instance = false; }

void Parameters::SetColliderEnergy(double ecms){ ecms_ = ecms; }
void Parameters::SetRenormalizationScale(double mur){ mur_ = mur; }
void Parameters::SetFactorizationScale(double muf){	muf_ = muf; }
void Parameters::SetTopQuarkMass(double m){	m_ = m; }
void Parameters::SetAlphaS(double coup) { alpha_s_ = coup; }

double Parameters::GetColliderEnergy() { return ecms_; }
double Parameters::GetRenormalizationScale() { return mur_; }
double Parameters::GetFactorizationScale() { return muf_; }
double Parameters::GetTopQuarkMass() { return m_; }
double Parameters::GetAlphaS() { return alpha_s_; }

