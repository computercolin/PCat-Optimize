#include "MyRJObject.h"
#include "ModelOptions.h"
#include "Data.h"
#include "MyDistribution.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

ModelOptions ModelOptions::instance;

ModelOptions::ModelOptions()
{

}

void ModelOptions::load(const char* modeloptions_file)
{
	fstream fin(modeloptions_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<modeloptions_file<<"."<<endl;
	fin>>nmax>>fixed;
	fin>>fluxlo>>fluxhi_min;
	fin>>fluxnorm>>norm_min>>norm_max>>midbin;
	fin>>smin>>smax>>slim;
	fin>>bg_min>>bg_max;
	int ntem = Data::get_instance().get_ntem();
	tem_min.assign(ntem, 0);
	tem_max.assign(ntem, 0);
	for (int i=0; i<Data::get_instance().get_ntem(); i++){
		fin>>tem_min[i]>>tem_max[i];
	}
	fin.close();
	cout<<"# Using "<<Data::get_instance().get_npix()<<" pixels"<<endl;
	cout<<"# Using "<<Data::get_instance().get_nbin()<<" energy bins"<<endl;
	cout<<"# Using "<<Data::get_instance().get_npsf()<<" PSF classes"<<endl;
	cout<<"# Using maximum number of sources "<<nmax<<" (fixed="<<fixed<<")"<<endl;
	cout<<"# Using flux lower limit between "<<fluxlo<<" ph cm^-2 s^-1"<<endl;
	cout<<"# Using flux upper limit above "<<fluxhi_min<<" ph cm^-2 s^-1"<<endl;
	cout<<"# Using normalization at "<<fluxnorm<<" ph cm^-2 s^-1 between "<<norm_min<<" and "<<norm_max<<endl;
	cout<<"# Applying flux limits to energy bin "<<midbin<<endl;
	cout<<"# Using PSF radii between "<<smin<<" and "<<smax<<", limiting PSF evaluation to "<<slim<<endl;
	cout<<"# Using background between "<<bg_min<<" and "<<bg_max<<" ph cm^-2 s^-1 sr^-1"<<endl;
	for (int i=0; i<Data::get_instance().get_ntem(); i++){
		cout<<"# Using template "<<i<<" between "<<tem_min[i]<<" and "<<tem_max[i]<<endl;
	}
}

MyRJObject<MyDistribution> ModelOptions::objects()
{
	Data data = Data::get_instance();

	// assumes objects have independent flux in each energy bin
	return MyRJObject<MyDistribution>(2+Data::get_instance().get_nbin(), nmax, fixed, MyDistribution(
			data.get_lmin(), data.get_lmax(),
			data.get_bmin(), data.get_bmax(),
                        fluxlo, fluxhi_min, fluxnorm, norm_min, norm_max, midbin));
}
