#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>

using namespace std;
using namespace DNest3;

MyModel::MyModel()
:objects(ModelOptions::get_instance().objects())
,image(Data::get_instance().get_nbin() * Data::get_instance().get_npsf() * Data::get_instance().get_npix())
,lambda(Data::get_instance().get_nbin() * Data::get_instance().get_npsf() * Data::get_instance().get_npix())
,staleness(100) // start stale
,bg_min(ModelOptions::get_instance().get_bg_min())
,bg_max(ModelOptions::get_instance().get_bg_max())
,bg(Data::get_instance().get_nbin()*Data::get_instance().get_npsf())
,tem_min(ModelOptions::get_instance().get_tem_min())
,tem_max(ModelOptions::get_instance().get_tem_max())
,tem(Data::get_instance().get_ntem())
,s_min(0.001) // fixed PSF scale, 0.05 degrees
,s_max(0.050) // 3 degrees
,score(Data::get_instance().get_nbin() * Data::get_instance().get_npsf())
//,gcore(Data::get_instance().get_nbin() * Data::get_instance().get_npsf())
,stail(Data::get_instance().get_nbin() * Data::get_instance().get_npsf())
,gtail(Data::get_instance().get_nbin() * Data::get_instance().get_npsf())
,fcore(Data::get_instance().get_nbin() * Data::get_instance().get_npsf())
{
}

void MyModel::fromPrior()
{
	objects.fromPrior();
	for (int i=0; i<Data::get_instance().get_nbin(); i++){
		for (int j=0; j<Data::get_instance().get_npsf(); j++){
			int jpsf = i*Data::get_instance().get_npsf() + j;
			bg[jpsf] = exp(log(bg_min) + log(bg_max/bg_min)*randomU());
			score[jpsf] = exp(log(s_min) + log(s_max/s_min)*randomU());
			//gcore[jpsf] = 1. / randomU();
			stail[jpsf] = 0;
			while (stail[jpsf] < score[jpsf]){
				stail[jpsf] = exp(log(s_min) + log(s_max/s_min)*randomU());
			}
			gtail[jpsf] = 1. / randomU();
			fcore[jpsf] = randomU();
		}
	}

	for (int i=0; i<Data::get_instance().get_ntem(); i++){
		tem[i] = exp(log(tem_min[i]) + log(tem_max[i]/tem_min[i])*randomU());
	}
	calculate_image();
}

void MyModel::update_lambdas(){
        const vector< double >& exposure = Data::get_instance().get_exposure();
        const vector< double >& etemplate = Data::get_instance().get_etemplate();
        // HARDCODED nside
        const double pixel_area = 4 * M_PI / (12. * 256. * 256.);
        
        for (int i=0; i<Data::get_instance().get_nbin(); i++){
        	for (int ii=0; ii<Data::get_instance().get_npsf(); ii++){
			int jpsf = i*Data::get_instance().get_npsf() + ii;
        	        for (int iii=0; iii<Data::get_instance().get_npix(); iii++){
        	                int jimg = i*Data::get_instance().get_npsf()*Data::get_instance().get_npix() + ii*Data::get_instance().get_npix() + iii;
                                lambda[jimg] = image[jimg] + bg[jpsf];
                                for (int iv=0; iv<Data::get_instance().get_ntem(); iv++){
          	                      int jetemplate = iv*Data::get_instance().get_nbin()*Data::get_instance().get_npix() + i*Data::get_instance().get_npix() + iii;
                                      lambda[jimg] += tem[iv]*etemplate[jetemplate];
                                }
                                lambda[jimg] *= exposure[jimg] * pixel_area;
			}
		}
	}
}

void MyModel::calculate_image()
{
	// Get coordinate stuff from data
	const vector< double >& l = Data::get_instance().get_l_rays();
	const vector< double >& b = Data::get_instance().get_b_rays();
	const vector< vector < double > >& unit_vectors = Data::get_instance().get_unit_vectors();

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
		image.assign(Data::get_instance().get_nbin()*Data::get_instance().get_npsf()*Data::get_instance().get_npix(), 0);
	}

	double lc, bc, M;
	double ll, bb, r;

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

		lc = (*component)[0]; bc = (*component)[1];
		// unit vector
		double xc = cos(bc) * cos(lc);
		double yc = cos(bc) * sin(lc);
		double zc = sin(bc);

		const double lim = 0.125; // 7 degrees = 0.122 rad
		const double coslim = cos(lim);
		for (int i=0; i<Data::get_instance().get_nbin(); i++){
			M = (*component)[2+i]; //assumes flux bins with independent intensities
			for (int ii=0; ii<Data::get_instance().get_npsf(); ii++){
				int jpsf = i*Data::get_instance().get_npsf() + ii;
				double sc = score[jpsf];
				double isc = 1./sc;
				double ist = 1./stail[jpsf];
				//double igc = 1./gcore[jpsf];
				double igt = 1./gtail[jpsf];
				for (int iii=0; iii<Data::get_instance().get_npix(); iii++){
					int jimg = i*Data::get_instance().get_npsf()*Data::get_instance().get_npix() + ii*Data::get_instance().get_npix() + iii;
					double cosdt = unit_vectors[iii][0] * xc + unit_vectors[iii][1] * yc + unit_vectors[iii][2] * zc;
					if (cosdt > coslim){
						image[jimg] += coeff*M*fcore[jpsf]/(2.*M_PI)*isc*isc*exp(-(1-cosdt)*isc*isc);
						image[jimg] += coeff*M*(1.-fcore[jpsf])/(2.*M_PI)*ist*ist*(1.-igt)*pow(1.+(1-cosdt)*igt*ist*ist,-gtail[jpsf]);
					}

					/*ll = fmod(l[iii] - lc + M_PI, 2*M_PI) - M_PI; // returns longitude difference [-pi, +pi]
					bb = b[iii] - bc;
					double maxb = max(fabs(b[iii]), fabs(bc));
					// only calculate r if within boxed limit
					if (fabs(bb) < lim && fabs(ll)*cos(maxb) < lim)
					{
						// Haversine formula
						double lH = sin(0.5 * ll);
						double bH = sin(0.5 * bb);
						double cc = cos(b[iii]) * cos(bc);
						r = 2. * asin(sqrt(bH*bH + cc*lH*lH));
						if (r*r < lim*lim)
						{
							image[jimg] += coeff*M*fcore[jpsf]/(2.*M_PI)*isc*isc*exp(-0.5*r*r*isc*isc);
							image[jimg] += coeff*M*(1.-fcore[jpsf])/(2.*M_PI)*ist*ist*(1.-igt)*pow(1.+0.5*igt*r*r*ist*ist,-gtail[jpsf]);
						}
					}*/
				}
			}
		}
	}

	update_lambdas();
}

double MyModel::perturb()
{
	double logH = 0.;
	// see comment below
	int oldN = objects.get_components().size();
	logH -= objects.get_dist().log_pn(oldN);

	double r = randomU();
	if(r <= 0.97)
	{
		logH += objects.perturb();
		calculate_image();
	}
	else if (r <= 0.98)
	{
		//for (int i=0; i<Data::get_instance().get_nbin()*Data::get_instance().get_npsf(); i++){
			int i = randInt(Data::get_instance().get_nbin()*Data::get_instance().get_npsf());
			bg[i] = log(bg[i]);
			bg[i] += log(bg_max/bg_min)*randh();
			bg[i] = mod(bg[i] - log(bg_min), log(bg_max/bg_min)) + log(bg_min);
			bg[i] = exp(bg[i]);
		//}
		update_lambdas();
	}
	else if (r <= 0.99)
	{
		//for (int i=0; i<Data::get_instance().get_ntem(); i++){
			int i = randInt(Data::get_instance().get_ntem());
			tem[i] = log(tem[i]);
			tem[i] += log(tem_max[i]/tem_min[i])*randh();
			tem[i] = mod(tem[i] - log(tem_min[i]), log(tem_max[i]/tem_min[i])) + log(tem_min[i]);
			tem[i] = exp(tem[i]);
		//}
		update_lambdas(); // point source image has not changed; circumvent RJObject added/removed component arrays
	}
	else
	{
		int jpsf = randInt(Data::get_instance().get_nbin()*Data::get_instance().get_npsf());
		//for (int jpsf=0; jpsf<Data::get_instance().get_nbin()*Data::get_instance().get_npsf(); jpsf++){
			score[jpsf] = log(score[jpsf]);
			score[jpsf] += log(s_max/s_min)*randh();
			score[jpsf] = mod(score[jpsf] - log(s_min), log(s_max/s_min)) + log(s_min);
			score[jpsf] = exp(score[jpsf]);
			//score[jpsf] = mod(score[jpsf] + randh(), 1.);
			//gcore[jpsf] = 1./gcore[jpsf];
			//gcore[jpsf] = mod(gcore[jpsf] + randh(), 1.);
			//gcore[jpsf] = 1./gcore[jpsf];

                        stail[jpsf] = log(stail[jpsf]);
                        stail[jpsf] += log(s_max/s_min)*randh();
                        stail[jpsf] = mod(stail[jpsf] - log(s_min), log(s_max/s_min)) + log(s_min);
                        stail[jpsf] = exp(stail[jpsf]);
                        gtail[jpsf] = 1./gtail[jpsf];
                        gtail[jpsf] = mod(gtail[jpsf] + randh(), 1.);
                        gtail[jpsf] = 1./gtail[jpsf];

			fcore[jpsf] = mod(fcore[jpsf] + randh(), 1.);

			if (stail[jpsf] < score[jpsf]){
				return -1E300;
			}
		//}
		staleness = 100; // force complete recalculation (and circumvent RJObject added/removed component arrays)
		calculate_image();
	}

	// this is HACKY, but it allows us to avoid editing RJObject
	// if number of components or hyperparameters change
	// need term from ratio of P(N | norm, fmin, fmax, a)
	// if only components have moved, this term is unnecessary
	// but should be zero
	int newN = objects.get_components().size();
	logH += objects.get_dist().log_pn(newN);

	return logH;
}

double MyModel::logLikelihood() const
{
	const vector< double >& data = Data::get_instance().get_image();

	double logL = 0.;
        for (int i=0; i<Data::get_instance().get_nbin()*Data::get_instance().get_npsf()*Data::get_instance().get_npix(); i++){
		logL += data[i]*log(lambda[i]) - lambda[i];
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
	for(int i=0; i<Data::get_instance().get_nbin()*Data::get_instance().get_npsf(); i++){
		out<<bg[i]<<' ';
	}
	for(int i=0; i<Data::get_instance().get_ntem(); i++){
		out<<tem[i]<<' ';
	}
	for(int i=0; i<Data::get_instance().get_nbin()*Data::get_instance().get_npsf(); i++){
		//out<<score[i]<<' '<<gcore[i]<<' '<<stail[i]<<' '<<gtail[i]<<' '<<fcore[i]<<' ';
		out<<score[i]<<' '<<stail[i]<<' '<<gtail[i]<<' '<<fcore[i]<<' ';
	}
	out<<staleness<<' ';
}

string MyModel::description() const
{
	return string("pixel lambdas | flux distribution fmin, fmax, norm, gamma | point sources | isotropic backgrounds | template coefficients | PSF score, gcore, stail, gtail | staleness");
}

