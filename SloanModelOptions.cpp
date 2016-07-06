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

void SloanModelOptions::load(const char* modeloptions_file)
{
	fstream fin(modeloptions_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<modeloptions_file<<"."<<endl;
	fin>>nmax>>fixed;
	fin>>fluxlo>>fluxhi_min;
	fin>>fluxnorm>>norm_min>>norm_max>>midbin;
	fin>>bg_min>>bg_max;
	cout<<"# Using "<<SloanData::get_instance().get_width()<<"x"<<SloanData::get_instance().get_height()<<" images"<<endl;
	cout<<"# Using "<<SloanData::get_instance().get_height()<<" bands"<<endl;
	cout<<"# Using maximum number of sources "<<nmax<<" (fixed="<<fixed<<")"<<endl;
	cout<<"# Using flux lower limit between "<<fluxlo<<" counts"<<endl;
	cout<<"# Using flux upper limit above "<<fluxhi_min<<" counts"<<endl;
	cout<<"# Using normalization at "<<fluxnorm<<" counts between "<<norm_min<<" and "<<norm_max<<endl;
	cout<<"# Applying flux limits to band "<<midbin<<endl;
	cout<<"# Using background between "<<bg_min<<" and "<<bg_max<<" ph cm^-2 s^-1 sr^-1"<<endl;
}

MyRJObject<MyDistribution> SloanModelOptions::objects()
{
	SloanData data = SloanData::get_instance();

	// assumes objects have independent flux in each energy bin
	return MyRJObject<MyDistribution>(2+data.get_nband(), nmax, fixed, MyDistribution(
			data.get_xmin(), data.get_xmax(),
			data.get_ymin(), data.get_ymax(), false,
                        fluxlo, fluxhi_min, fluxnorm, norm_min, norm_max,
			data.get_nband(), midbin));
}
