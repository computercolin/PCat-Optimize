#include "MyDistribution.h"
#include <RandomNumberGenerator.h>
#include <Utils.h>
#include <cmath>
#include "Data.h"
#include <boost/math/special_functions/erf.hpp>
#include <iostream>

using namespace DNest3;

MyDistribution::MyDistribution(double l_min, double l_max,
					double b_min, double b_max,
					double fluxlo_min, double fluxlo_max,
					double fluxhi_min,
					double flux_norm,
					double norm_min, double norm_max)
:l_min(l_min)
,l_max(l_max)
,b_min(b_min)
,b_max(b_max)
,fluxlo_min(fluxlo_min)
,fluxlo_max(fluxlo_max)
,fluxhi_min(fluxhi_min)
,flux_norm(flux_norm)
,norm_min(norm_min)
,norm_max(norm_max)
,sdev_color_scale(0.64)
,mean_colors(Data::get_instance().get_nbin() - 1)
,sdev_colors(Data::get_instance().get_nbin() - 1)
{

}

void MyDistribution::fromPrior()
{
	//fluxlo = exp(log(fluxlo_min) + log(fluxlo_max/fluxlo_min)*randomU());
	fluxlo = fluxlo_min; // fixed fmin
	fluxhi = fluxhi_min / randomU();
	norm = exp(log(norm_min) + log(norm_max/norm_min)*randomU());
	gamma = randomU();
	for (int i=0; i<Data::get_instance().get_nbin()-1; i++){
		mean_colors[i] = tan(M_PI * (randomU() - 0.5)); // full Cauchy
		sdev_colors[i] = sdev_color_scale * tan(0.5 * M_PI * randomU()); //half Cauchy
	}
}

double MyDistribution::perturb_parameters()
{
	double logH = 0.;

	int which = randInt(4) + 1;

	if(which == 0)
	{
		/*fluxlo = log(fluxlo);
		fluxlo += log(fluxlo_max/fluxlo_min)*randh();
		fluxlo = mod(fluxlo - log(fluxlo_min),
			log(fluxlo_max/fluxlo_min)) + log(fluxlo_min);
		fluxlo = exp(fluxlo);*/
	}
        else if(which == 1)
        {
                fluxhi = fluxhi_min / fluxhi;
                fluxhi += randh();
                fluxhi = mod(fluxhi, 1.);
                fluxhi = fluxhi_min / fluxhi;
        }
	else if(which == 2)
	{
		norm = log(norm);
		norm += log(norm_max/norm_min)*randh();
		norm = mod(norm - log(norm_min),
			log(norm_max/norm_min)) + log(norm_min);
		norm = exp(norm);
	}
	else if(which == 3)
	{
		gamma += randh();
		gamma = mod(gamma, 1.);
	}
	else if(which == 4){
		int i = randInt(Data::get_instance().get_nbin()-1);
		mean_colors[i] = atan(mean_colors[i]) / M_PI + 0.5;
		mean_colors[i] += randh();
		mean_colors[i] = mod(mean_colors[i], 1.);
		mean_colors[i] = tan(M_PI * (mean_colors[i] - 0.5));
		sdev_colors[i] = 2 * atan(sdev_colors[i] / sdev_color_scale) / M_PI;
		sdev_colors[i] += randh();
		sdev_colors[i] = mod(sdev_colors[i], 1.);
		sdev_colors[i] = sdev_color_scale * tan(0.5 * M_PI * sdev_colors[i]);
	}


	return logH;
}

// kind of hacky... returns P(N | this distribution)
double MyDistribution::log_pn(const int n) const
{
	// alpha = a - 1
	double alpha = 1./gamma - 1;
	double navg = norm * pow(flux_norm / fluxlo, alpha) * (1 - pow(fluxlo / fluxhi, alpha)) / alpha;
	return n * log(navg) - navg - lgamma(n + 1);
}

// x, y, flux
double MyDistribution::log_pdf(const std::vector<double>& vec) const
{
	// alpha = a - 1
	double alpha = 1./gamma - 1.;
	if(vec[0] < l_min || vec[0] > l_max || vec[1] < b_min || vec[1] > b_max){
		return -1E300;
	}

	int midbin = Data::get_instance().get_nbin() / 2;

	double logp = 0;
	for (int i=0; i<Data::get_instance().get_nbin(); i++){
		int ii = (i < midbin) ? i : i - 1;
		if (i == midbin){
			if (vec[2+i] < fluxlo || vec[2+i] > fluxhi){
				return -1E300;
			}
			logp += log(alpha) + alpha*log(fluxlo) - (alpha+1)*log(vec[2+i]) - log(1 - pow(fluxlo / fluxhi, alpha));
		}
		else{
			if (vec[2+i] < 0){
				return -1E300;
			}
			double colour = vec[2+i] / vec[2+midbin];
			// constant term matters because we change the colour parameters
			logp -= log(vec[2+i]) + 0.5*log(2. * M_PI * sdev_colors[ii] * sdev_colors[ii]) + pow(log(colour) - mean_colors[ii], 2)/(2.*sdev_colors[ii]*sdev_colors[ii]);
		}
	}
 
	// no spatial distribution part since sources don't change in drag

	return logp;
}

double MyDistribution::angular_area() const
{
	return (sin(b_max) - sin(b_min)) * (l_max - l_min);
}

void MyDistribution::from_uniform(std::vector<double>& vec) const
{
	// alpha = a - 1
	double alpha = 1./gamma - 1.;

	vec[0] = l_min + (l_max - l_min)*vec[0];
	vec[1] = asin(sin(b_min) + (sin(b_max) - sin(b_min))*vec[1]);
	int midbin = Data::get_instance().get_nbin() / 2;
	vec[2+midbin] = fluxlo*pow(1. + (pow(fluxlo/fluxhi, alpha) - 1.) * vec[2+midbin], -1./alpha);
	for (int i=0; i<Data::get_instance().get_nbin(); i++){
		int ii = (i < midbin) ? i : i - 1;
		if (i != midbin) {
			double colour;
			//if (vec[2+i] > 1e-10){ //kludgey?
				colour = exp(mean_colors[ii] + sqrt(2*sdev_colors[ii]*sdev_colors[ii]) * boost::math::erf_inv(2*vec[2+i] - 1));
			//}
			//else{
			//	colour = 0;
			//}
			vec[2+i] = colour * vec[2+midbin];
		}
	}
}

void MyDistribution::to_uniform(std::vector<double>& vec) const
{
	// alpha = a - 1
	double alpha = 1./gamma - 1.;

	vec[0] = (vec[0] - l_min)/(l_max - l_min);
	vec[1] = (sin(vec[1]) - sin(b_min))/(sin(b_max) - sin(b_min));
	int midbin = Data::get_instance().get_nbin() / 2;
	for (int i=0; i<Data::get_instance().get_nbin(); i++){
		int ii = (i < midbin) ? i : i - 1;
		if (i != midbin){
			double colour = vec[2+i] / vec[2+midbin];
			vec[2+i] = 0.5*(1 + erf((log(colour) - mean_colors[ii])/sqrt(2*sdev_colors[ii]*sdev_colors[ii])));
		}
	}
	vec[2+midbin] = (1. - pow(fluxlo/vec[2+midbin], alpha))/(1. - pow(fluxlo/fluxhi, alpha));
}

void MyDistribution::print(std::ostream& out) const
{
	//out<<fluxlo<<' '<<fluxhi<<' '<<norm<<' '<<gamma<<' ';
	for (int i=0; i<Data::get_instance().get_nbin() - 1; i++){
		out<<mean_colors[i]<<' '<<sdev_colors[i]<<' ';
	}
}

