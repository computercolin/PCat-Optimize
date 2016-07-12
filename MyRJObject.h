#ifndef _MyRJObject_
#define _MyRJObject_

#include <RJObject.h>
#include <iostream>

template<class Distribution>
class MyRJObject : public RJObject<Distribution> {
	public:
		MyRJObject(int num_dimensions, int max_num_components, bool fixed, const Distribution& dist);
		double perturb();
		void print(std::ostream& out) const;
};

template<class Distribution>
MyRJObject<Distribution>::MyRJObject(int num_dimensions, int max_num_components, bool fixed, const Distribution& dist)
	: RJObject<Distribution>(num_dimensions, max_num_components, fixed, dist){

}

template<class Distribution>
double MyRJObject<Distribution>::perturb()
{
	if(this->max_num_components == 0)
		return 0.;

	this->added.resize(0);
	this->removed.resize(0);

	double logH = 0.;

	int which = (this->fixed)?(1 + DNest3::randInt(2)):(DNest3::randInt(3)); //leaving out merges/splits for now

	if(which == 0)
	{
		logH -= this->get_dist().log_pn(this->num_components);
		// Do some birth or death
		logH += this->perturb_num_components(
				pow(10., - 6.*DNest3::randomU()));
		logH += this->get_dist().log_pn(this->num_components);
	}
	else if(which == 1)
	{
		// Change the hyperparameters
		logH -= this->get_dist().log_pn(this->num_components);
		if(DNest3::randomU() <= 0.97)
		{
			logH += this->dist.perturb1(this->components, this->u_components);
		}
		else
		{
			this->removed = this->components;
			logH += this->dist.perturb2(this->components, this->u_components);
			this->added = this->components;
		}
		logH += this->get_dist().log_pn(this->num_components);
	}
	else if(which == 2)
	{
		logH += this->perturb_components(pow(10., - 6.*DNest3::randomU()));
	}

	return logH;
}

template<class Distribution>
void MyRJObject<Distribution>::print(std::ostream& out) const
{
	this->dist.print(out); out<<' ';
	out<<this->num_components<<' ';

	// Write out components
	for(int j=0; j<this->num_dimensions; j++)
	{
		for(int i=0; i<this->num_components; i++)
			out<<this->components[i][j]<<' ';

		// Pad with zeros (turned-off components)
		for(int i=this->num_components; i<this->max_num_components; i++)
			out<<0.<<' ';
	}
}

#endif
