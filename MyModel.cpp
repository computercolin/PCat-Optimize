#include <boost/align/aligned_allocator.hpp>
#include <cmath>
#include <immintrin.h>
#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"

using namespace std;
using namespace DNest3;

MyModel::MyModel(MyRJObject<MyDistribution> objects, int nbin, int npsf, int npix, double pixel_area,
                 vec_align32<double> data, vec_align32<double> exposure, double bg_min, double bg_max,
                 int ntem, vec_align32<double> tem_min, vec_align32<double> tem_max, vec_align32<double> etemplate)
:globals(&MyModelGlobals::get_instance())
,objects(objects)
,nbin(nbin)
,npsf(npsf)
,npix(npix)
,pixel_area(pixel_area)
,image(nbin*npsf*npix)
,lambda(nbin*npsf*npix)
,staleness(100) // start stale
,bg_min(bg_min)
,bg_max(bg_max)
,bg(nbin*npsf)
,ntem(ntem)
,tem(ntem)
{
    globals->data = data;
    globals->exposure = exposure;
    globals->tem_min = tem_min;
    globals->tem_max = tem_max;
    globals->etemplate = etemplate;
}


void MyModel::fromPrior()
{
	objects.fromPrior();
	for (int i=0; i<nbin; i++){
		for (int j=0; j<npsf; j++){
			int jpsf = i*npsf + j;
			bg[jpsf] = exp(log(bg_min) + log(bg_max/bg_min)*randomU());
		}
	}

	for (int i=0; i<ntem; i++){
		tem[i] = exp(log(globals->tem_min[i]) + log(globals->tem_max[i]/globals->tem_min[i])*randomU());
	}
	calculate_image();
}

void MyModel::update_lambdas() {
    int bin_x_psf_max = nbin * npsf;
    int nbin_npsf_npix_prod = nbin * npsf * npix;
    double *image0 = &image[0];
    double *lambda0 = &lambda[0];
    double *exposure0= &globals->exposure[0];
    double *tem0 = &tem[0];
    double *etemplate0 = &globals->etemplate[0];
    __m256d pixel_area_pd = _mm256_set1_pd(pixel_area);
    for (int bin_x_psf = 0; bin_x_psf < bin_x_psf_max; bin_x_psf++) {
        int binPsfNpix_offset = bin_x_psf * npix;
        double bgval = bg[bin_x_psf];
        __m256d bgval_pd = _mm256_set1_pd(bgval);
        int bin_x_psf_x_pix = binPsfNpix_offset;
        int bin_x_psf_x_pix_last = binPsfNpix_offset + npix -1;
        int bin_x_psf_x_pix_avx_max = (bin_x_psf_x_pix_last -3) & ~0b11;
        for (int i_tem = 0; i_tem < ntem; i_tem++) {
            int i_jetemplate = i_tem * nbin_npsf_npix_prod + bin_x_psf_x_pix;
            for (; bin_x_psf_x_pix <= bin_x_psf_x_pix_avx_max; bin_x_psf_x_pix+=4) {
                __m256d lam_pd = _mm256_load_pd(image0 + bin_x_psf_x_pix);
                lam_pd = _mm256_add_pd(lam_pd, bgval_pd);
                __m256d tmp = _mm256_mul_pd(_mm256_set1_pd(*(tem0 + i_tem)),
                                            _mm256_load_pd(etemplate0 + i_jetemplate));
                lam_pd = _mm256_add_pd(lam_pd, tmp);
                lam_pd = _mm256_mul_pd(lam_pd, _mm256_load_pd(exposure0 + bin_x_psf_x_pix));
                lam_pd = _mm256_mul_pd(lam_pd, pixel_area_pd);
                _mm256_store_pd(lambda0 + bin_x_psf_x_pix, lam_pd);
            }

            // Loop peel. Not 4-doubles wide, process individually without avx.
            for (; bin_x_psf_x_pix <= bin_x_psf_x_pix_last; bin_x_psf_x_pix++) {
                lambda[bin_x_psf_x_pix] = image[bin_x_psf_x_pix] + bgval;
                lambda[bin_x_psf_x_pix] += *(tem0 + i_tem) * *(etemplate0 + i_jetemplate);
                lambda[bin_x_psf_x_pix] *= *(exposure0 + bin_x_psf_x_pix) * pixel_area;
            }
        }
    }

//    for (int bin_x_psf = 0; bin_x_psf < bin_x_psf_max; bin_x_psf++) {
//        int binpsfnpix_offset = bin_x_psf * npix;
//        double bgval = bg[bin_x_psf];
//        double *image_pos= &image[binpsfnpix_offset];
//        double *exposure_pos = &globals->exposure[binpsfnpix_offset];
//        double *lambda_pos = &lambda[binpsfnpix_offset];
//        double *lambda_last = lambda_pos + npix -1;
//        double *lambda_last_chunk = (double *) ((uintptr_t) (lambda_last - 3) & ~0b11111);
//        int bin_x_psf_x_pix = binpsfnpix_offset;
//        for (; lambda_pos <= lambda_last_chunk; bin_x_psf_x_pix+=4, lambda_pos+=4, image_pos+=4, exposure_pos+=4) {
//            __m256d lam_pd = _mm256_add_pd(_mm256_load_pd(lambda_pos), _mm256_set1_pd(bgval));
//            for (int i_tem = 0; i_tem < ntem; i_tem++) {
//                int i_jetemplate = i_tem * nbin_npsf_npix_prod + bin_x_psf_x_pix;
//                lam_pd = _mm256_mul_pd(lam_pd, _mm256_load_pd(&globals->etemplate[i_jetemplate]));
//            }
//            lam_pd = _mm256_mul_pd(lam_pd, _mm256_load_pd(exposure_pos));
//            lam_pd = _mm256_mul_pd(lam_pd, _mm256_set1_pd(pixel_area));
//            _mm256_store_pd(lambda_pos, lam_pd);
//        }
//
//        // Loop peel, final data not 256d packable.
//        for (; lambda_pos <= lambda_last; bin_x_psf_x_pix++, lambda_pos++, image_pos++, exposure_pos++) {
//            lambda[bin_x_psf_x_pix] = image[bin_x_psf_x_pix] + bgval;
//            for (int i_tem = 0; i_tem < ntem; i_tem++) {
//                int jetemplate = i_tem * nbin_npsf_npix_prod + bin_x_psf_x_pix;
//                lambda[bin_x_psf_x_pix] += tem[i_tem] * globals->etemplate[jetemplate];
//            }
//            lambda[bin_x_psf_x_pix] *= globals->exposure[bin_x_psf_x_pix] * pixel_area;
//        }
//    }
}

void MyModel::calculate_image()
{
	// Components
	const vector< vector<double> >& components = objects.get_components();

	// Diff
	const vector< vector<double> >& added = objects.get_added();
	const vector< vector<double> >& removed = objects.get_removed();
	bool update = (added.size() + removed.size() < components.size()) &&
				(staleness < 100);
	int num = (update)?(added.size() + removed.size()):(components.size());

	if(update)
		staleness++;
	else
	{
		staleness = 0;
		// Zero the image
		image.assign(nbin*npsf*npix, 0);
	}

	double xc, yc, M;

	const vector<double>* component;
	double coeff = 1.; // switches to -1 for removing components

	for(int k=0; k<num; k++)
	{
		if(update)
		{
			if(k < (int)added.size())
			{
				component = &(added[k]);
				coeff = 1.;
			}
			else
			{
				component = &(removed[k - (int)added.size()]);
				coeff = -1.;
			}
		}
		else
			component = &(components[k]);

		xc = (*component)[0]; yc = (*component)[1];

		for (int i=0; i<nbin; i++){
			M = (*component)[2+i]; //assumes components are x,y,fluxes
			for (int ii=0; ii<npsf; ii++){
				add_source_flux(i, ii, xc, yc, coeff*M);
			}
		}
	}

	update_lambdas();
}

double MyModel::perturb()
{
	double logH = 0.;

	double r = randomU();
	if(r <= 0.97)
	{
		logH += objects.perturb();
		calculate_image();
	}
	else if (r <= 0.98)
	{
		int i = randInt(nbin*npsf);
		bg[i] = log(bg[i]);
		bg[i] += log(bg_max/bg_min)*randh();
		bg[i] = mod(bg[i] - log(bg_min), log(bg_max/bg_min)) + log(bg_min);
		bg[i] = exp(bg[i]);
		update_lambdas();
	}
	else if (r <= 0.99)
	{
		int i = randInt(ntem);
		tem[i] = log(tem[i]);
		tem[i] += log(globals->tem_max[i]/globals->tem_min[i])*randh();
		tem[i] = mod(tem[i] - log(globals->tem_min[i]), log(globals->tem_max[i]/globals->tem_min[i])) + log(globals->tem_min[i]);
		tem[i] = exp(tem[i]);
		update_lambdas(); // point source image has not changed; circumvent RJObject added/removed component arrays
	}
	else
	{
		// perturb PSF
		staleness = 100; // force complete recalculation (and circumvent RJObject added/removed component arrays)
		calculate_image();
	}

	return logH;
}

double MyModel::logLikelihood() const
{
	double logL = 0.;
        for (int i=0; i<nbin*npsf*npix; i++){
		logL += pixelLogLikelihood(globals->data[i], lambda[i]);
	}

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	out<<setprecision(5);
	for(size_t i=0; i<image.size(); i++)
		out<<lambda[i]<<' ';
	out<<setprecision(10);
	objects.print(out); out<<' ';
	int i_lim = nbin*npsf;
	for(int i=0; i<i_lim; i++){
		out<<bg[i]<<' ';
	}
	for(int i=0; i<ntem; i++){
		out<<tem[i]<<' ';
	}
	out<<staleness;
}

string MyModel::description() const
{
	return string("pixel lambdas |"+objects.get_dist().description()+" | N, point sources | isotropic backgrounds | template coefficients | staleness");
}
