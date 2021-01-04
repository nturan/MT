#include "histogram.h"

Histogram::Histogram(std::string name, std::pair<double, double> limits, unsigned int number_of_bins,
	ParameterFunction parameter, std::vector<std::string> particles) {
	name_ = name;
	limits_ = limits;
	number_of_bins_ = number_of_bins;
	bin_width_ = (limits_.second - limits_.first) / number_of_bins_;
	for (unsigned int i = 0; i < number_of_bins_; ++i) {
		Bin new_bin;
		new_bin.x = limits_.first + bin_width_ * i + 0.5 * bin_width_;
		new_bin.n = 0;
		new_bin.v = 0.0;
		new_bin.e2 = 0.0;
		data_.push_back(new_bin);
	}
	underflow_.x = 0.0;
	underflow_.n = 0;
	underflow_.v = 0.0;
	underflow_.e2 = 0.0;

	overflow_.x = 0.0;
	overflow_.n = 0;
	overflow_.v = 0.0;
	overflow_.e2 = 0.0;

	parameter_ = parameter;
	particles_ = particles;
}

void Histogram::Fill(PhaseSpaceGenerator PS, double weight){
	std::map<std::string, double*> momenta;
	momenta["top"]  = PS.k1_;
	momenta["atop"] = PS.k2_;
	momenta["jet"]  = PS.k3_;
	double p[4] = { 0.0 };
	for (auto i = particles_.begin(); i != particles_.end(); ++i) {
		p[0] += momenta[*i][0];
		p[1] += momenta[*i][1];
		p[2] += momenta[*i][2];
		p[3] += momenta[*i][3];
	}
	double x = parameter_(p);

	double v = weight / bin_width_;
	if (x < limits_.first) {
		underflow_.n  += 1;
		underflow_.v  += v;
		underflow_.e2 += v * v;
	}
	else if (x < limits_.second) {
		unsigned int i = static_cast<int>(trunc((x - limits_.first) / bin_width_));
		data_.at(i).n += 1;
		data_.at(i).v  += v;
		data_.at(i).e2 += v * v;
	}
	else {
		overflow_.n += 1;
		overflow_.v += v;
		overflow_.e2 += v * v;
	}
}

void Histogram::Print(){
	std::setprecision(17);
	std::cout << "#START: HISTOGRAM" << std::endl;
	std::cout << "NAME: " << name_ << std::endl;
	std::cout << "FROM: " << limits_.first << std::endl;
	std::cout << "TO: " << limits_.second << std::endl;
	std::cout << "UNDERFLOW: " << underflow_.n << " " << underflow_.v << " " << std::sqrt(underflow_.e2) << std::endl;
	std::cout << "DATA_START:" << std::endl;
	for (auto it = data_.begin(); it != data_.end(); ++it) {
		std::cout << it->x << " " << it->n << " " << it->v << " " << std::sqrt(it->e2) << std::endl;
	}
	std::cout << "DATA_END " << std::endl;
	std::cout << "OVERFLOW: " << overflow_.n << " " << overflow_.v << " " << std::sqrt(overflow_.e2) << std::endl;
	std::cout << "#END: HISTOGRAM" << std::endl;
}
