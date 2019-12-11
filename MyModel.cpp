#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace DNest3;

MyModel::MyModel(MyRJObject<MyDistribution> objects, int nbin, int npsf, int npix, double pixel_area,
	vector<double> data, vector<double> exposure, double bg_min, double bg_max,
	int ntem, vector<double> tem_min, vector<double> tem_max, vector<double> etemplate)
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

void MyModel::update_lambdas(){
        for (int i=0; i<nbin; i++){
        	for (int ii=0; ii<npsf; ii++){
			int jpsf = i*npsf + ii;
        	        for (int iii=0; iii<npix; iii++){
        	                int jimg = jpsf*npix + iii;
                                lambda[jimg] = image[jimg] + bg[jpsf];
                                for (int iv=0; iv<ntem; iv++){
          	                      int jetemplate = iv*nbin*npsf*npix + jimg;
                                      lambda[jimg] += tem[iv]*globals->etemplate[jetemplate];
                                }
                                lambda[jimg] *= globals->exposure[jimg] * pixel_area;
			}
		}
	}
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
