#include <cmath>
#include <vector>
#include <immintrin.h>
#include "SloanModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "SloanData.h"
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
        DNest3::vec_align32<double>(SloanData::get_instance().get_nband() * SloanData::get_instance().get_height() * SloanData::get_instance().get_width(), SloanData::get_instance().get_exposure()),
        SloanModelOptions::get_instance().get_bg_min(),
        SloanModelOptions::get_instance().get_bg_max(),
        1, // templates:
        DNest3::vec_align32<double>(1, 0.99), // who
        DNest3::vec_align32<double>(1, 1.01), // needs
        DNest3::vec_align32<double>(SloanData::get_instance().get_nband() * SloanData::get_instance().get_height() * SloanData::get_instance().get_width(), 0)
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
    double *dataPos = &globals->data[0];
    double *dataLast = &globals->data[nbin * npsf * npix -1];
    double *dataLastChunk = (double *) ((uintptr_t) (dataLast - 3) & ~0b11111);
    const double *lambdaPos = &this->lambda[0];
    __m256d gain_pd = _mm256_set1_pd(gain_inv);
    __m256d neg2_pd = _mm256_set1_pd(-2);
    __m256d bias_pd = _mm256_set1_pd(bias);
    __m128d logL_sd = _mm_setzero_pd();
    for (;dataPos <= dataLastChunk;) {
        __m256d _data = _mm256_load_pd(dataPos);
        __m256d numerator = _mm256_sub_pd(_mm256_load_pd(lambdaPos), _data);
        numerator = _mm256_mul_pd(numerator, numerator);
        __m256d denominator = _mm256_mul_pd(_mm256_mul_pd(neg2_pd, gain_pd), _mm256_sub_pd(_data, bias_pd));
        numerator = _mm256_div_pd(numerator, denominator);

        // Sum all 4 packed doubles, accumulate on logL_sd.
        numerator = _mm256_hadd_pd(numerator, numerator);
        __m128d tmp = _mm256_extractf128_pd(numerator, 1);
        tmp = _mm_add_sd(_mm256_castpd256_pd128(numerator), tmp);
        logL_sd = _mm_add_sd(logL_sd, tmp);

        dataPos += 4; lambdaPos += 4;
    }

    double logL = -1;
    _mm_store_sd(&logL, logL_sd);
    for (; dataPos <= dataLast; dataPos++, lambdaPos++) {
        logL_sd += std::pow(*lambdaPos - *dataPos, 2) / (-2 * gain_inv * (*dataPos - bias));
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
