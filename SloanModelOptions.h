#ifndef _SloanModelOptions_
#define _SloanModelOptions_

#include "MyRJObject.h"
#include "MyDistribution.h"

class SloanModelOptions
{
	private:
		// if fixed: nmax is number of sources
		// else: nmax is max number of sources
		int nmax;
		bool fixed;

		// lower flux limit (counts)
		double fluxlo;
		// lower bound on upper flux limit (counts)
		double fluxhi_min;
		// flux where norm = dN/dlogF defined
		double fluxnorm;
		// bounds on norm
		double norm_min, norm_max;
		// which bin to apply lower flux limit to
		int midbin;
		// bounds on background (counts / pixel)
		double bg_min, bg_max;

	public:
		SloanModelOptions();
		void load(const char* modeloptions_file);
		MyRJObject<MyDistribution> objects();
		double get_bg_min() { return bg_min; }
		double get_bg_max() { return bg_max; }

	// Singleton
	private:
		static SloanModelOptions instance;
	public:
		static SloanModelOptions& get_instance() { return instance; }
};

#endif
