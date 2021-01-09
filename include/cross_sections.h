#pragma once

#include <map>
#include "parameters.h"
#include "phase_space_generator.h"
#include "helper.h"
#include "histogram.h"

typedef double PartonicCrossSection(PhaseSpaceGenerator,Parameters*);

typedef double HadronicCrossSection(std::map<std::string, double> v, double& wgt, 
	Parameters* p, std::vector<Histogram*>* histograms);

typedef std::map<std::string, std::function<double(PhaseSpaceGenerator, Parameters*)>> PartonicContribution;


PartonicCrossSection NoContribution;
namespace lo {
	HadronicCrossSection Hadronic;
	namespace partonic {
		PartonicCrossSection gg, qqb, qbq; 

		extern PartonicContribution Born;
	}
}

namespace nlo {
	HadronicCrossSection Hadronic2, Hadronic3;
	namespace partonic {
		extern PartonicContribution Hard, Soft, Virt, Coll_left_z, Coll_right_z, Coll_1, Coll_0;
		namespace hard {
			PartonicCrossSection gg, qqb, qbq, qg, gq, qbg, gqb;
		}
		namespace soft {
			PartonicCrossSection gg, qqb, qbq;
		}
		namespace virt {
			PartonicCrossSection gg, qqb, qbq;
		}
		namespace coll_left_z {
			PartonicCrossSection gg, qqb, qg;
		}
		namespace coll_right_z {
			PartonicCrossSection qg;
		}
		namespace coll_1 {
			PartonicCrossSection gg, qqb;
		}
		namespace coll_0 {
			PartonicCrossSection gg, qqb;
		}
		
	}
}

