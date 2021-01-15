#pragma once

#include <map>
#include "parameters.h"
#include "phase_space_generator.h"
#include "helper.h"
#include "histogram.h"
#include "n_dimensional_integral.h"

typedef double PartonicCrossSection(PhaseSpaceGenerator,Parameters*);

typedef double HadronicCrossSection(std::map<std::string, double> v, double& wgt, 
	Parameters* p, std::vector<Histogram*>* histograms);

typedef std::map<std::string, std::function<double(PhaseSpaceGenerator, Parameters*)>> PartonicContribution;


PartonicCrossSection NoContribution;

extern std::map<std::string, std::pair<double, double>> lo_variables, 
nlo_2_variables, nlo_3_variables, nlo_z_variables, nlo_j_variables;

extern std::map<std::string, Parameters*> parameter_sets;
extern std::map<std::string, std::vector<Histogram*>*> histogram_sets;
extern std::map<std::string, std::vector<Integral*>*> integrals;

extern std::map<std::string, Integrand> lo_integrands, lo_nlo_2_integrands,
										nlo_2_integrands, nlo_3_integrands;

extern void InitializeIntegrands( std::vector<std::string> histogram_strings, double ecms, double m );
extern void ExecuteIntegralsAndPrintResults( std::vector<std::string> perturbation_order, int iterations, int calls );

namespace lo {
	HadronicCrossSection Hadronic;
	namespace partonic {
		PartonicCrossSection gg, qqb, qbq; 

		extern PartonicContribution Born;
	}
}

namespace nlo {
	HadronicCrossSection Hadronic2WithBorn, Hadronic2, Hadronic3, HadronicConstZ, HadronicZ;
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

