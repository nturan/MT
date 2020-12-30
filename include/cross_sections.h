#pragma once

#include <map>
#include "amps.h"
#include "parameters.h"
#include "phase_space_generator.h"
#include "helper.h"

typedef double PartonicCrossSection(PhaseSpaceGenerator,Parameters&);

typedef double HadronicCrossSection(std::map<std::string, double> v, double* wgt, Parameters& p);

typedef std::map<std::string, std::function<double(PhaseSpaceGenerator, Parameters&)>> PartonicContribution;


PartonicCrossSection NoContribution;
namespace sigma::leading_order {
	HadronicCrossSection Hadronic;
	namespace partonic {
		extern PartonicContribution Born;
		PartonicCrossSection gg, qqb, qbq; 
	}
}

namespace sigma::next_to_leading_order {
	HadronicCrossSection Hadronic2, Hadronic3;
	namespace partonic {
		extern PartonicContribution Soft, Virt, Coll_left_z, Coll_right_z, Coll_1, Coll_0;
		namespace soft {
			PartonicCrossSection gg, qqb, qbq;
		}
		namespace virt {
			PartonicCrossSection gg, qqb, qbq;
		}
		namespace coll_left_z {
			PartonicCrossSection gg, qqb, qbq, qg;
		}
		namespace coll_right_z {
			PartonicCrossSection qg;
		}
		namespace coll_1 {
			PartonicCrossSection gg, qqb, qbq;
		}
		namespace coll_0 {
			PartonicCrossSection gg, qqb, qbq;
		}
		
	}
}

