#pragma once
#include <map>
#include <string>

class PhaseSpaceGenerator {

public:
	PhaseSpaceGenerator(std::map<std::string, double> var);
	double x1_, x2_, s_, p1_[4], p2_[4], k1_[4], k2_[4], k3_[4], dGamma_;
};