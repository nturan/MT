#pragma once

#include "cross_sections.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

typedef std::tuple<double, double, double, double, double> Event;
double LoMaxIntegrandFinderIntegral(std::map<std::string, double> v, double& wgt,
	Parameters* p, std::vector<Histogram*>* histograms, double* max_integrand);

class EventGenerator {
private:
	void GenerateLOEvents(unsigned int number, Parameters* p);
	void GenerateNLOEvents(unsigned int number, Parameters* p);
	std::vector<Event> events;
public:
	EventGenerator(unsigned int number, std::string perturbation_order, Parameters* p);
	EventGenerator(std::string file_name);
	std::vector<Event> GetEvents();
	void Print();

};