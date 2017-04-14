#include "MyRJObject.h"
#include "SloanModelOptions.h"
#include "SloanData.h"
#include "MyDistribution.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

SloanModelOptions SloanModelOptions::instance;

SloanModelOptions::SloanModelOptions()
{

}

void SloanModelOptions::load(ptree pt)
{
	nmax = pt.get<int>("modeloptions.nmax");
	fixed = pt.get<bool>("modeloptions.nfixed");
	fluxlo = pt.get<double>("modeloptions.fmin");
	fluxhi_min = pt.get<double>("modeloptions.fmax_lo");
	fluxnorm = pt.get<double>("modeloptions.fnorm");
	norm_min = pt.get<double>("modeloptions.norm_lo");
	norm_max = pt.get<double>("modeloptions.norm_hi");
	midbin = pt.get<int>("modeloptions.midbin");
	bg_min = pt.get<double>("modeloptions.background_lo");
	bg_max = pt.get<double>("modeloptions.background_hi");
	penalty = pt.get<double>("modeloptions.penalty", 0);
	slope = pt.get<double>("modeloptions.slope", 0);

	cout<<"# Using "<<SloanData::get_instance().get_width()<<"x"<<SloanData::get_instance().get_height()<<" images"<<endl;
	if (fixed) {
		cout<<"# Using exactly "<<nmax<<" sources"<<endl;
	}
	else{
		cout<<"# Using maximum number of sources "<<nmax<<endl;
	}
	cout<<"# Using flux lower limit between "<<fluxlo<<" counts"<<endl;
	cout<<"# Using flux upper limit above "<<fluxhi_min<<" counts"<<endl;
	cout<<"# Using normalization at "<<fluxnorm<<" counts between "<<norm_min<<" and "<<norm_max<<endl;
	cout<<"# Applying flux limits to band "<<midbin<<endl;
	cout<<"# Using background between "<<bg_min<<" and "<<bg_max<<" ph cm^-2 s^-1 sr^-1"<<endl;
	if (penalty) {
		cout<<"# Using a penalty for each source of "<<penalty<<endl;
	}
	if (slope) {
		cout<<"# Fixing flux distribution slope to "<<slope<<endl;
	}
}

MyRJObject<MyDistribution> SloanModelOptions::objects()
{
	SloanData data = SloanData::get_instance();

	// assumes objects have independent flux in each energy bin
	return MyRJObject<MyDistribution>(2+data.get_nband(), nmax, fixed, MyDistribution(
			data.get_xmin(), data.get_xmax(),
			data.get_ymin(), data.get_ymax(), false,
                        fluxlo, fluxhi_min, fluxnorm, norm_min, norm_max,
			penalty, slope,
			data.get_nband(), midbin));
}
