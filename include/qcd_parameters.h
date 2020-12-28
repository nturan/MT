#pragma once

#include <vector>
#include <string>

class Parameters {
private:
	static bool has_instance;
	static Parameters* parameters;
	Parameters();

	double ecms_, mur_, muf_, m_, alpha_s_;
public:
	static Parameters* GetInstance();
	void SetColliderEnergy(double ecms);
	void SetRenormalizationScale(double mur);
	void SetFactorizationScale(double muf);
	void SetTopQuarkMass(double m);
	void SetAlphaS(double coup);
	double GetColliderEnergy();
	double GetRenormalizationScale();
	double GetFactorizationScale();
	double GetTopQuarkMass();
	double GetAlphaS();
	std::vector<std::string> channels;
	~Parameters();
 };
