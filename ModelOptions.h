#ifndef _ModelOptions_
#define _ModelOptions_

#include "MyRJObject.h"
#include "MyDistribution.h"

class ModelOptions
{
	private:
		// if fixed: nmax is number of sources
		// else: nmax is max number of sources
		int nmax;
		bool fixed;

		// lower flux limit (ph cm^-2 s^-1)
		double fluxlo;
		// lower bound on upper flux limit (ph cm^-2 s^-1)
		double fluxhi_min;
		// flux where norm = dN/dlogF defined
		double fluxnorm;
		// bounds on norm
		double norm_min, norm_max;
		// which bin to apply lower flux limit to
		int midbin;
		// bounds on PSF radii (in radians)
		double smin, smax;
		// bound on PSF evaluation
		double slim;
		// bounds on background (ph cm^-2 s^-1 sr^-1)
		double bg_min, bg_max;
		// bounds on templates
		std::vector<double> tem_min;
		std::vector<double> tem_max;

	public:
		ModelOptions();
		void load(const char* modeloptions_file);
		MyRJObject<MyDistribution> objects();
		double get_smin() { return smin; }
		double get_smax() { return smax; }
		double get_slim() { return slim; }
		double get_bg_min() { return bg_min; }
		double get_bg_max() { return bg_max; }
		std::vector<double>& get_tem_min() { return tem_min; }
		std::vector<double>& get_tem_max() { return tem_max; }


	// Singleton
	private:
		static ModelOptions instance;
	public:
		static ModelOptions& get_instance() { return instance; }
};

#endif
