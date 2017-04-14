#ifndef _MyDistribution_
#define _MyDistribution_

#include <vector>
#include <Distribution.h>

// inherit from Distribution to get hyperparameter change and drag moves for free
class MyDistribution:public Distribution
{
	private:
		// Limits
		double l_min, l_max, b_min, b_max;
		bool radec; // sets whether cos(b) factor is used
		double fluxlo;
		double fluxhi_min;
		double flux_norm, norm_min, norm_max;

		double penalty;
		double slope;

		// Lower limit and 1/slope for Pareto interim prior
		// for masses
		double fluxhi, norm;
		double gamma;

		int nbin;
		int midbin;
		double sdev_color_scale;
		std::vector< double > mean_colors;
		std::vector< double > sdev_colors;

		double perturb_parameters();

	public:
		MyDistribution(double l_min, double l_max,
					double b_min, double b_max, bool radec,
					double fluxlo,
					double fluxhi_min,
					double flux_norm,
					double norm_min, double norm_max,
					double penalty, double slope,
					int nbin, int midbin);

		void fromPrior();

		double log_pn(const int n) const;

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
		std::string description() const;
};

#endif

