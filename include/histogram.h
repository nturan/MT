#pragma once

#include "phase_space_generator.h"
#include <sstream>
#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <map>

typedef std::function<double(double[4])> ParameterFunction;
struct Bin {
	int n;
	double x;
	double v;
	double e2;
};

class Histogram {
private:
	std::string name_;
	std::pair<double, double> limits_;
	unsigned int number_of_bins_;
	double bin_width_;
	ParameterFunction parameter_;
	std::vector<Bin> data_;
	std::vector<std::string> particles_;
	Bin underflow_, overflow_;
	double ecms_, mur_, muf_, m_;
	std::string pdf_name_;
	std::map<std::string, ParameterFunction> parameters_{
		{"PT",  &calcTransversalMomentum},
		{"E", &calcEnergy},
		{"Y", &calcRapidity},
		{"ETA", &calcPseudoRapidity},
		{"PHI", &calcAzimuth},
		{"M", &calcInvariantMass}
	};
public:
	Histogram(std::string name, std::pair<double, double> limits, unsigned int number_of_bins, 
		ParameterFunction parameter, std::vector<std::string> particles, Parameters& p);
	Histogram(std::string initialization_string, Parameters& p);
	~Histogram();
	void Fill(PhaseSpaceGenerator PS, double weight);
	void Print();

};