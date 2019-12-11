#include "SloanModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "SloanData.h"
#include <cmath>
#include <vector>
#include "SloanModelOptions.h"

using namespace std;
using namespace DNest3;

SloanModel::SloanModel()
	: MyModel(
        SloanModelOptions::get_instance().objects(),
        SloanData::get_instance().get_nband(),
        1, //only allowing one PSF class for now
        SloanData::get_instance().get_height() * SloanData::get_instance().get_width(),
        1., //all quantities will be per pixel
        SloanData::get_instance().get_image(),
        // exposure is all ones
        vector<double>(SloanData::get_instance().get_nband()*SloanData::get_instance().get_height()*SloanData::get_instance().get_width(), SloanData::get_instance().get_exposure()),
        SloanModelOptions::get_instance().get_bg_min(),
        SloanModelOptions::get_instance().get_bg_max(),
        1, // templates:
        vector<double>(1, 0.99), // who
        vector<double>(1, 1.01), // needs
        vector<double>(SloanData::get_instance().get_nband()*SloanData::get_instance().get_height()*SloanData::get_instance().get_width(), 0)
	     ) // them?
    // ,sloanGlobals(&SloanModelGlobals::get_instance())
	,width(SloanData::get_instance().get_width())
	,height(SloanData::get_instance().get_height())
	,psf_size(SloanData::get_instance().get_psf_size())
	,psf_resampling(SloanData::get_instance().get_psf_resampling())
	,psfs(SloanData::get_instance().get_psfs())
	,bias(SloanData::get_instance().get_bias())
	,gain(SloanData::get_instance().get_gain())
    ,gain_inv(1 / SloanData::get_instance().get_gain())
{
}

void SloanModel::add_source_flux(int ibin, int ipsf, double xc, double yc, double M)
{
	// check that ipsf = 0?
	(void)ipsf;
	// assert that psf_size is odd?
	// determine pixels to evaluate over to stay within PSF box
	int bound = (psf_size - 1) / 2;
	int xmin = xc > bound - 1 ? ceil(xc - bound) : 0;
	int xmax = xc < width - bound - 1 ? floor(xc + bound) + 1 : width;
	int ymin = yc > bound - 1 ? ceil(yc - bound) : 0;
	int ymax = yc < height - bound - 1 ? floor(yc + bound) + 1 : height;

	for (int y=ymin; y<ymax; y++){
		for (int x=xmin; x<xmax; x++){
			int jimg = ibin*width*height + y*width + x;
			double psf_x = (x - xc + bound)*psf_resampling;
			double psf_y = (y - yc + bound)*psf_resampling;
			int f_psf_x = floor(psf_x);
			int c_psf_x = ceil(psf_x);
			int f_psf_y = floor(psf_y);
			int c_psf_y = ceil(psf_y);
			// maybe an actual interpolation routine would be faster?
			double psf = psfs[f_psf_y*psf_size*psf_resampling + f_psf_x] * (c_psf_x - psf_x) * (c_psf_y - psf_y)
					+ psfs[f_psf_y*psf_size*psf_resampling + c_psf_x] * (psf_x - f_psf_x) * (c_psf_y - psf_y)
					+ psfs[c_psf_y*psf_size*psf_resampling + f_psf_x] * (c_psf_x - psf_x) * (psf_y - f_psf_y)
					+ psfs[c_psf_y*psf_size*psf_resampling + c_psf_x] * (psf_x - f_psf_x) * (psf_y - f_psf_y);
			image[jimg] += M*psf;
		}
	}
}

double SloanModel::logLikelihood() const
{
    double logL = 0.;
    int maxInd = nbin*npsf*npix;
    for (int i=0; i<maxInd; i++){
        logL += std::pow(lambda[i] - globals->data[i], 2)
                 / (-2 * gain_inv * (globals->data[i] - bias));
    }

    return logL;
}

double SloanModel::pixelLogLikelihood(double data, double lambda) const
{
	// http://classic.sdss.org/dr7/algorithms/fluxcal.html
	// double variance = (data - bias)/gain; //ignoring dark noise and read noise
    // double variance = (data - bias) * gain_inv; //ignoring dark noise and read noise

	// if variance calculated on data, not lambda
	// then normalization term constant
	// return -(lambda - data)*(lambda - data)/(2 * variance);
    return std::pow(lambda - data, 2)/(-2 * gain_inv * (data - bias));

	//VARIANCE CALCULATED ON LAMBDA
	//double variance = (lambda-bias)/gain;
	//if (variance < 1.) { variance = 1.; }
	// the sqrt(2 pi) stays constant
	//return -0.5*log(variance)-(lambda - data)*(lambda - data)/(2 * variance);
}
