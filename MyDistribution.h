#ifndef _MyDistribution_
#define _MyDistribution_

#include <vector>
#include <Distribution.h>

// Hyperparameters setting interim prior for galaxy properties
class MyDistribution:public Distribution
{
	private:
		// Limits
		double l_min, l_max, b_min, b_max;
		double fluxlo;
		double fluxhi_min;
		double flux_norm, norm_min, norm_max;

		// Lower limit and 1/slope for Pareto interim prior
		// for masses
		double fluxhi, norm;
		double gamma;

		int midbin;
		double sdev_color_scale;
		std::vector< double > mean_colors;
		std::vector< double > sdev_colors;

		double perturb_parameters();

	public:
		MyDistribution(double l_min, double l_max,
					double b_min, double b_max,
					double fluxlo,
					double fluxhi_min,
					double flux_norm,
					double norm_min, double norm_max,
					int midbin);

		void fromPrior();

		double log_pn(const int n) const;

		double log_pdf(const std::vector<double>& vec) const;
		double angular_area() const; // for splits and merges
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
};

#endif

