#pragma once

#include <vector>
#include <string>

namespace qcd {
	class Parameters {
	public:
		double ecms_, mur2_, muf2_, m_;
		std::vector<std::string> channels;
  };
}