#pragma once
#include <map>
#include <string>
#include "parameters.h"
#include "helper.h"

class PhaseSpaceGenerator {
private:
	void generatePS2(std::map<std::string, double> var);
	void generatePS3(std::map<std::string, double> var);
public:
	PhaseSpaceGenerator(std::map<std::string, double> var, Parameters &p);
	double x1_, x2_, z_, s_, p1_[4], p2_[4], k1_[4], k2_[4], k3_[4], dGamma_;
	double m_;
	double ecms_;
	double xmin_;
	bool not_collinear_, not_soft_;
};