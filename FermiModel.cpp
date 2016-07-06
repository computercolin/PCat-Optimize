#include "FermiModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "FermiData.h"
#include <cmath>
#include <vector>
#include "FermiModelOptions.h"

using namespace std;
using namespace DNest3;

FermiModel::FermiModel()
	: MyModel(
	FermiModelOptions::get_instance().objects(),
	FermiData::get_instance().get_nbin(),
	FermiData::get_instance().get_npsf(),
	FermiData::get_instance().get_npix(),
	FermiData::get_instance().get_pixel_area(),
	FermiData::get_instance().get_image(),
	FermiData::get_instance().get_exposure(),
	FermiModelOptions::get_instance().get_bg_min(),
	FermiModelOptions::get_instance().get_bg_max(),
	FermiData::get_instance().get_ntem(),
	FermiModelOptions::get_instance().get_tem_min(),
	FermiModelOptions::get_instance().get_tem_max(),
	FermiData::get_instance().get_etemplate())
	,s_min(FermiModelOptions::get_instance().get_smin())
	,s_max(FermiModelOptions::get_instance().get_smax())
	,lim(FermiModelOptions::get_instance().get_slim())
	,score(nbin*npsf)
	,stail(nbin*npsf)
	,gtail(nbin*npsf)
	,fcore(nbin*npsf)
{
}

void FermiModel::fromPrior()
{
	for (int i=0; i<nbin; i++){
		for (int j=0; j<npsf; j++){
			int jpsf = i*npsf + j;
			score[jpsf] = exp(log(s_min) + log(s_max/s_min)*randomU());
			stail[jpsf] = 0;
			while (stail[jpsf] < score[jpsf]){
				stail[jpsf] = exp(log(s_min) + log(s_max/s_min)*randomU());
			}
			gtail[jpsf] = 1. / randomU();
			fcore[jpsf] = randomU();
		}
	}

	MyModel::fromPrior();
}

void FermiModel::add_source_flux(int ibin, int ipsf, double lc, double bc, double M)
{
	// Get coordinate stuff from data
	// maybe this fetch is expensive?
	const vector< vector < double > >& unit_vectors = FermiData::get_instance().get_unit_vectors();

	double xc = cos(bc) * cos(lc);
	double yc = cos(bc) * sin(lc);
	double zc = sin(bc);

	double coslim = cos(lim[ibin]);
	int jpsf = ibin*npsf + ipsf;
	double sc = score[jpsf];
	double isc = 1./sc;
	double ist = 1./stail[jpsf];
	double igt = 1./gtail[jpsf];
	for (int iii=0; iii<npix; iii++){
		int jimg = jpsf*npix + iii;
		double cosdt = unit_vectors[iii][0] * xc + unit_vectors[iii][1] * yc + unit_vectors[iii][2] * zc;
		if (cosdt > coslim){
			image[jimg] += M*fcore[jpsf]/(2.*M_PI)*isc*isc*exp(-(1-cosdt)*isc*isc);
			image[jimg] += M*(1.-fcore[jpsf])/(2.*M_PI)*ist*ist*(1.-igt)*pow(1.+(1-cosdt)*igt*ist*ist,-gtail[jpsf]);
		}
	}
}

double FermiModel::perturb()
{
	double r = randomU();
	if(r <= 0.99)
	{
		return MyModel::perturb();
	}
	else
	{
		int jpsf = randInt(nbin*npsf);
		score[jpsf] = log(score[jpsf]);
		score[jpsf] += log(s_max/s_min)*randh();
		score[jpsf] = mod(score[jpsf] - log(s_min), log(s_max/s_min)) + log(s_min);
		score[jpsf] = exp(score[jpsf]);
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
		staleness = 100; // force complete recalculation (and circumvent RJObject added/removed component arrays)
		calculate_image();

		return 0;
	}
}

double FermiModel::pixelLogLikelihood(double data, double lambda) const
{
	return data*log(lambda) - lambda;
}

void FermiModel::print(std::ostream& out) const
{
	MyModel::print(out);

	for(int i=0; i<nbin*npsf; i++){
		out<<' '<<score[i]<<' '<<stail[i]<<' '<<gtail[i]<<' '<<fcore[i];
	}
}

string FermiModel::description() const
{
	return string(MyModel::description()+" | score, stail, gtail, fcore for each bin x PSF");
}

