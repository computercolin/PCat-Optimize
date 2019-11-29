#ifndef _MyModel_
#define _MyModel_

#include <Model.h>
#include <vector>
#include "MyRJObject.h"
#include "MyDistribution.h"
#include "MyModelGlobals.h"

using namespace std;

class MyModel : public DNest3::Model
{
	protected:
        MyModelGlobals* globals;
		MyRJObject<MyDistribution> objects;

		int nbin, npsf, npix;
		double pixel_area;
		// The model source flux image
		std::vector<double> image;
		// not sure what visibility the next two methods should have
		void update_lambdas();
		void calculate_image(); // this should be protected so FermiModel PSF changes can recalculate images from scratch
		// The model total image
		std::vector<double> lambda;

		virtual void add_source_flux(int ibin, int ipsf, double xc, double yc, double M) = 0;
		virtual double pixelLogLikelihood(double data, double lambda) const = 0;
		// How many steps since image was computed from scratch
		int staleness;

		// isotropic background
		double bg_min, bg_max; // same prior bounds for all energy bins
		std::vector<double> bg;

                // emission templates
                int ntem;
                std::vector<double> tem; // coefficients
	public:
		MyModel(MyRJObject<MyDistribution> objects, int nbin, int npsf, int npix, double pixel_area,
			vector<double> data, vector<double> exposure, double bg_min, double bg_max,
			int ntem, vector<double> tem_min, vector<double> tem_max, vector<double> etemplate);

		// Generate the point from the prior
		void fromPrior();

		// Metropolis-Hastings proposals
		double perturb();

		// Likelihood function
		double logLikelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;

		// MyModel(const MyModel& o) = delete;  // Copy constructor needed for datastructure initalization.
		// MyModel& operator=(const MyModel& o) = delete;
};

#endif
