#include "MyRJObject.h"
#include "FermiModelOptions.h"
#include "FermiData.h"
#include "MyDistribution.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

FermiModelOptions FermiModelOptions::instance;

FermiModelOptions::FermiModelOptions()
{

}

void FermiModelOptions::load(ptree pt)
{
	nmax = pt.get<int>("modeloptions.nmax");
	fixed = pt.get<bool>("modeloptions.nfixed");
	fluxlo = pt.get<double>("modeloptions.fmin");
	fluxhi_min = pt.get<double>("modeloptions.fmax_lo");
	fluxnorm = pt.get<double>("modeloptions.fnorm");
	norm_min = pt.get<double>("modeloptions.norm_lo");
	norm_max = pt.get<double>("modeloptions.norm_hi");
	midbin = pt.get<int>("modeloptions.midbin");
	smin = pt.get<double>("modeloptions.sigma_lo");
	smax = pt.get<double>("modeloptions.sigma_hi");
	BOOST_FOREACH(ptree::value_type &v, pt.get_child("modeloptions.eval_lims"))
		slim.push_back(v.second.get_value<double>());
	bg_min = pt.get<double>("modeloptions.background_lo");
	bg_max = pt.get<double>("modeloptions.background_hi");
	BOOST_FOREACH(ptree::value_type &v, pt.get_child("modeloptions.template_los"))
		tem_min.push_back(v.second.get_value<double>());
	BOOST_FOREACH(ptree::value_type &v, pt.get_child("modeloptions.template_his"))
		tem_max.push_back(v.second.get_value<double>());

	cout<<"# Using "<<FermiData::get_instance().get_npix()<<" pixels"<<endl;
	cout<<"# Using "<<FermiData::get_instance().get_nbin()<<" energy bins"<<endl;
	cout<<"# Using "<<FermiData::get_instance().get_npsf()<<" PSF classes"<<endl;
	cout<<"# Using maximum number of sources "<<nmax<<" (fixed="<<fixed<<")"<<endl;
	cout<<"# Using flux lower limit between "<<fluxlo<<" ph cm^-2 s^-1"<<endl;
	cout<<"# Using flux upper limit above "<<fluxhi_min<<" ph cm^-2 s^-1"<<endl;
	cout<<"# Using normalization at "<<fluxnorm<<" ph cm^-2 s^-1 between "<<norm_min<<" and "<<norm_max<<endl;
	cout<<"# Applying flux limits to energy bin "<<midbin<<endl;
	cout<<"# Using PSF radii between "<<smin<<" and "<<smax<<", limiting PSF evaluation to";
	for (int i=0; i<FermiData::get_instance().get_nbin(); i++){
		cout<<" "<<slim[i];
	}
	cout<<endl;
	cout<<"# Using background between "<<bg_min<<" and "<<bg_max<<" ph cm^-2 s^-1 sr^-1"<<endl;
	for (int i=0; i<FermiData::get_instance().get_ntem(); i++){
		cout<<"# Using template "<<i<<" between "<<tem_min[i]<<" and "<<tem_max[i]<<endl;
	}
}

MyRJObject<MyDistribution> FermiModelOptions::objects()
{
	FermiData data = FermiData::get_instance();

	// assumes objects have independent flux in each energy bin
	return MyRJObject<MyDistribution>(2+FermiData::get_instance().get_nbin(), nmax, fixed, MyDistribution(
			data.get_lmin(), data.get_lmax(),
			data.get_bmin(), data.get_bmax(), true,
                        fluxlo, fluxhi_min, fluxnorm, norm_min, norm_max, 0., 0.,
			FermiData::get_instance().get_nbin(), midbin));
}
