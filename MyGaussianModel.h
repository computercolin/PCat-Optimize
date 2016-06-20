#include "MyModel.h"
#include <vector>

class MyGaussianModel:public MyModel
{
	private:
		double gain_min, gain_max;
		std::vector<double> gain;

	public:
		MyGaussianModel();

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

