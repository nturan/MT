#pragma once

#include <vector>
#include <string>

namespace qcd {
	class Parameters {
	private:
		static bool has_instance;
		static Parameters* parameters;
		Parameters();

		double ecms_, mur2_, muf2_, m_;
	public:
		static Parameters* GetInstance();
		void SetColliderEnergy(double ecms);
		void SetRenormalizationScale(double mur);
		void SetFactorizationScale(double muf);
		void SetTopQuarkMass(double m);
		double GetColliderEnergy();
		double GetSquaredRenormalizationScale();
		double GetSquaredFactorizationScale();
		double GetTopQuarkMass();
		std::vector<std::string> channels;
		~Parameters();
  };
}