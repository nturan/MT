#include "parameters.h"

Parameters::Parameters( std::string pdf_name ) {
	InitializePartonDistributionFunctionSets(pdf_name);
	std::cout << "Dont forget to:\n";
	std::cout << " o SetColliderEnergy()\n";
	std::cout << " o SetRenormalizationScale()\n";
	std::cout << " o SetFactorizationScale()\n";
	std::cout << " o SetTopQuarkMass()\n";
	std::cout << std::endl;
}

void Parameters::SetColliderEnergy(double ecms){ ecms_ = ecms; }
void Parameters::SetRenormalizationScale(double mur) { mur2_ = mur * mur; }

void Parameters::SetFactorizationScale(double muf){	
	muf2_     = muf*muf; 
	alpha_s_ = lhapdf->alphasQ2(muf2_);
	gs2_ = alpha_s_ * 4.0 * M_PI;

	using namespace std::placeholders;
	Fs["gg"] = std::bind(&Parameters::gg, this, _1, _2, muf2_);
	Fs["qqb"] = std::bind(&Parameters::qqb, this, _1, _2, muf2_);
	Fs["qbq"] = std::bind(&Parameters::qbq, this, _1, _2, muf2_);
	Fs["gq"] = std::bind(&Parameters::gq, this, _1, _2, muf2_);
	Fs["qg"] = std::bind(&Parameters::qg, this, _1, _2, muf2_);
	Fs["gqb"] = std::bind(&Parameters::gqb, this, _1, _2, muf2_);
	Fs["qbg"] = std::bind(&Parameters::qbg, this, _1, _2, muf2_);
}
void Parameters::SetTopQuarkMass(double m){	m_ = m; }
void Parameters::SetCutParameter(double xmin) { xmin_ = xmin; }

double Parameters::GetColliderEnergy() { return ecms_; }
double Parameters::GetSquaredRenormalizationScale() { return mur2_; }
double Parameters::GetSquaredFactorizationScale() { return muf2_; }
double Parameters::GetTopQuarkMass() { return m_; }
double Parameters::GetAlphaS() { return alpha_s_; }
double Parameters::GetSquaredGs() { return gs2_; }
double Parameters::GetCutParameter() { return xmin_; }

void Parameters::InitializePartonDistributionFunctionSets(std::string pdf_name)
{
	std::cout << "Initializing pdf sets...\n";
	try {
		lhapdfset = LHAPDF::getPDFSet(pdf_name);
	}
	catch (const LHAPDF::ReadError&) {
		std::cerr << "Cannot load PDF set from LHAPDF." << std::endl;
		std::cerr << "Please check the PDF name " << pdf_name << " and"
			<< std::endl;
		std::cerr << "ensure its presence in your LHAPDF setup." << std::endl;
		exit(1);
	}
	lhapdf = lhapdfset.mkPDF(0);
}


double Parameters::gg(double x1, double x2, double muf2) {
	if (x1 > 1.0) { return 0.0; }
	if (x2 > 1.0) { return 0.0; }
	return lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(21, x2, muf2) / x1 / x2;
}

double Parameters::qqb(double x1, double x2, double muf2) {
	if (x1 > 1.0) { return 0.0; }
	if (x2 > 1.0) { return 0.0; }
	return (lhapdf->xfxQ2(+1, x1, muf2) * lhapdf->xfxQ2(-1, x2, muf2) +
			lhapdf->xfxQ2(+2, x1, muf2) * lhapdf->xfxQ2(-2, x2, muf2) +
			lhapdf->xfxQ2(+3, x1, muf2) * lhapdf->xfxQ2(-3, x2, muf2) +
			lhapdf->xfxQ2(+4, x1, muf2) * lhapdf->xfxQ2(-4, x2, muf2) +
			lhapdf->xfxQ2(+5, x1, muf2) * lhapdf->xfxQ2(-5, x2, muf2)) / x1 / x2;
}

double Parameters::qbq(double x1, double x2, double muf2) {
	return qqb(x2, x1, muf2);
}

double Parameters::gq(double x1, double x2, double muf2) {
	if (x1 > 1.0) { return 0.0; }
	if (x2 > 1.0) { return 0.0; }

	return (lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(+1, x2, muf2) +
			lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(+2, x2, muf2) +
			lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(+3, x2, muf2) +
			lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(+4, x2, muf2) +
			lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(+5, x2, muf2)) / x1 / x2;
}

double Parameters::qg(double x1, double x2, double muf2) {
	return gq(x2, x1, muf2);
}

double Parameters::gqb(double x1, double x2, double muf2) {
	if (x1 > 1.0) { return 0.0; }
	if (x2 > 1.0) { return 0.0; }
	return (lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(-1, x2, muf2) +
			lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(-2, x2, muf2) +
			lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(-3, x2, muf2) +
			lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(-4, x2, muf2) +
			lhapdf->xfxQ2(21, x1, muf2) * lhapdf->xfxQ2(-5, x2, muf2)) / x1 / x2;
}

double Parameters::qbg(double x1, double x2, double muf2) {
	return gqb(x2, x1, muf2);
}
