#pragma once
/*
	This file was written by Turan Nuraliyev.
	It includes implementation of the Histogram class. 
	The Histogram class allows to create fill and print 
	weighted histograms during vegas integration.
	
	nuraliyev.turan@gmail.com
*/

#include "phase_space_generator.h"
#include <sstream>
#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <map>

typedef std::function<double(double[4])> ParameterFunction;
struct Bin { int n;	double x; double v; double e2; };

class Histogram {
public:
	/**
	 * @brief Standard constructor for Histogram
	 * @param name a string name for the histogram to print later
	 * @param limits start and end points as std::pair
	 * @param number_of_bins number of bins to create as unsigned int 
	 * @param parameter parameter function as std::function taking double 
	 * 4 arrays and returning double
	 * @param particles names of particles as vector of strings: top, atop, jet
	 * @param p theory parameters as pointer to an instance of the Parameters class
	*/
	Histogram(std::string name, 
		std::pair<double, double> limits, 
		unsigned int number_of_bins,
		ParameterFunction parameter, 
		std::vector<std::string> particles, 
		Parameters* p);
	/**
	 * @brief compact constructor for Histogram
	 * @param initialization_string histogram creation string as in madanalysis,
	 * e.g. "PT(top) 50 0 500"
	 * @param p theory parameters as pointer to an instance of the Parameters class
	*/
	Histogram(std::string initialization_string, Parameters* p);
	/**
	 * @brief Function to fill Histogram for given phase space point and weight
	 * @param PS phase space point as instance of the PhaseSpaceGenerator
	 * @param weight event weight as double
	*/
	void Fill(PhaseSpaceGenerator PS, double weight);
	void Print();
	void Clear();
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


};