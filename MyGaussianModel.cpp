#include "MyModel.h"
#include "MyGaussianModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>

using namespace std;
using namespace DNest3;

MyGaussianModel::MyGaussianModel()
:MyModel()
,gain_min(0.01)
,gain_max(10.)
,gain(Data::get_instance().get_nbin() * Data::get_instance().get_npsf())
{
}

void MyGaussianModel::fromPrior()
{
	for (int i=0; i<Data::get_instance().get_nbin()*Data::get_instance().get_npsf(); i++){
		gain[i] = exp(log(gain_min) + log(gain_max/gain_min)*randomU());
	}
	MyModel::fromPrior();
}

double MyGaussianModel::perturb()
{
	double r = randomU();

	if (r <= 0.99){
		return MyModel::perturb();
	}
	else{
		int i = randInt(Data::get_instance().get_nbin()*Data::get_instance().get_npsf());
		gain[i] = log(gain[i]);
		gain[i] += log(gain_max/gain_min)*randh();
		gain[i] = mod(gain[i] - log(gain_min), log(gain_max/gain_min)) + log(gain_min);
		gain[i] = exp(gain[i]);
		return 0; // proposal respects prior
	}
}

double MyGaussianModel::logLikelihood() const
{
	const vector< double >& data = Data::get_instance().get_image();

	double logL = 0.;
        for (int i=0; i<Data::get_instance().get_nbin()*Data::get_instance().get_npsf(); i++){
		for (int j=0; j<Data::get_instance().get_npix(); j++){
			int k = i*Data::get_instance().get_npix() + j;
			logL -= 0.5*log(gain[i]) + 0.5*log(this->lambda[k]) + (data[k] - this->lambda[k])*(data[k] - this->lambda[k])/(2*gain[i]*this->lambda[k]);
		}
	}

	return logL;
}

void MyGaussianModel::print(std::ostream& out) const
{
	MyModel::print(out);
	for (int i=0; i<Data::get_instance().get_nbin()*Data::get_instance().get_npsf(); i++){
		out<<gain[i]<<' ';
	}
}

string MyGaussianModel::description() const
{
	return string("pixel lambdas | ndim, maxN, flux distribution fmax, norm, gamma | color mean, sdev | N, point sources | isotropic backgrounds | template coefficients | PSF score, stail, gtail, fcore | staleness | gains");
}

