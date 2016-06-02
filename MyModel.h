#ifndef _MyModel_
#define _MyModel_

#include <Model.h>
#include <vector>
#include "MyRJObject.h"
#include "ModelOptions.h"
#include "MyDistribution.h"

class MyModel:public DNest3::Model
{
	private:
		MyRJObject<MyDistribution> objects;

		// The model source flux image
		std::vector<long double> image;
		void update_lambdas();
		void calculate_image();
		// The model total image
		std::vector<long double> lambda;

		// How many steps since image was computed from scratch
		int staleness;

		// isotropic background
		double bg_min, bg_max; // same prior bounds for all energy bins
		std::vector<double> bg;

                // emission templates
                std::vector<double> tem_min;
		std::vector<double> tem_max;
                std::vector<double> tem;

		// PSF
		double s_min, s_max; // same prior bounds for all energy bins and PSF classes
		std::vector<double> score;
		//std::vector<double> gcore;
		std::vector<double> stail;
		std::vector<double> gtail;
		std::vector<double> fcore;
	public:
		MyModel();

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
};

#endif

