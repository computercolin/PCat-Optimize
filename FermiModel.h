#ifndef _FermiModel_
#define _FermiModel_

#include <vector>
#include "MyModel.h"

class FermiModel : public MyModel
{
	protected:
		// PSF
		double s_min, s_max; // same prior bounds for all energy bins and PSF classes
		std::vector<double> lim;
		std::vector<double> score;
		std::vector<double> stail;
		std::vector<double> gtail;
		std::vector<double> fcore;

		void add_source_flux(int ibin, int ipsf, double lc, double bc, double M);
		double pixelLogLikelihood(double data, double lambda) const;
	public:
		FermiModel();

		void fromPrior();
		double perturb();

		void print(std::ostream& out) const;
		// Return string with column information
		std::string description() const;
};

#endif

