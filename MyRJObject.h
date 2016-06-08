#ifndef _MyRJObject_
#define _MyRJObject_

#include <RJObject.h>
#include <iostream>

template<class Distribution>
class MyRJObject : public RJObject<Distribution> {
	public:
		MyRJObject(int num_dimensions, int max_num_components, bool fixed, const Distribution& dist);
		double perturb();

	private:
		double perturb_mergesplit();
		std::vector< std::vector<int> > find_pairs(double maxkickl, double maxkickb);
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
		// Do some birth or death
		logH += this->perturb_num_components(
				pow(10., - 6.*DNest3::randomU()));
	}
	else if(which == 1)
	{
		// Change the hyperparameters
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
	}
	else if(which == 2)
	{
		logH += this->perturb_components(pow(10., - 6.*DNest3::randomU()));
	}
	else if(which == 3)
	{
		logH += perturb_mergesplit();
	}

	return logH;
}

//FIXME 
template<class Distribution>
double MyRJObject<Distribution>::perturb_mergesplit(){
	double logH = 0;
	double maxkickb = 0.035 / 0.23;
	double maxkickl = 0.035; // about 2 deg
	int midbin = (this->num_dimensions - 2) / 2;
	//split
	if (DNest3::randomU() <= 0.5){
		// can't split if no sources or at max number
		if (this->num_components == 0 || this->num_components == this->max_num_components){
			return -1E300;
		}

		int i = DNest3::randInt(this->num_components);

		double kick_l = (DNest3::randomU()*2 - 1) * maxkickl;
		double kick_b = (DNest3::randomU()*2 - 1) * maxkickb;
		double fractn = DNest3::randomU();

		//create new component
		std::vector<double> newcomponent(this->num_dimensions, 0.);
		newcomponent[0] = this->components[i][0] + (1-fractn) * kick_l;
		newcomponent[1] = this->components[i][1] + (1-fractn) * kick_b;
		for (int j=2; j<this->num_dimensions; j++){
			if (j != 2 + midbin){
				newcomponent[j] = DNest3::randomU() * this->components[i][j];
			}
			else {
				newcomponent[j] = fractn * this->components[i][j];
			}
		}

		// perturb component i
		std::vector<double> modcomponent(this->num_dimensions, 0.);
		modcomponent[0] = this->components[i][0] - fractn * kick_l;
		modcomponent[1] = this->components[i][1] - fractn * kick_b;
		for (int j=2; j<this->num_dimensions; j++){
			// true for the midbin too
			modcomponent[j] = this->components[i][j] - newcomponent[j];
		}

		double modprior = this->dist.log_pdf(modcomponent);
		double newprior = this->dist.log_pdf(newcomponent);
		if (modprior > -1E300 && newprior > -1E300){
			double oldprior = this->dist.log_pdf(this->components[i]);
			logH += newprior + modprior - oldprior;
			logH += log(cos(modcomponent[1])) + log(cos(newcomponent[1])) - log(cos(this->components[i][1]));
			logH += log(2*maxkickl) + log(2*maxkickb) - log(this->dist.angular_area());
			logH += log(this->num_components) + log(this->num_components+1);
			for (int j=2; j<this->num_dimensions; j++){
				logH += log(this->components[i][j]);
			}
			
			this->num_components++;
			this->removed.push_back(this->components[i]);

			this->components[i].assign(modcomponent.begin(), modcomponent.end());
			this->added.push_back(modcomponent);
			this->dist.to_uniform(modcomponent);
			this->u_components[i].assign(modcomponent.begin(), modcomponent.end());

			this->components.push_back(newcomponent);
			this->added.push_back(newcomponent);
			this->dist.to_uniform(newcomponent);
			this->u_components.push_back(newcomponent);

			std::vector< std::vector<int> > pairs = find_pairs(maxkickl, maxkickb);
			logH -= log(2 * pairs.size()); // multiply by two to get ordered pairs
			return logH;
		}
		else{
			return -1E300;
		}
	}
	//merge
	else{
		// can't merge with fewer than two sources
		if (this->num_components < 2){
			return -1E300;
		}
		// count interacting pairs and choose one, with j > i
		std::vector< std::vector<int> > pairs = find_pairs(maxkickl, maxkickb);
		//if no interacting pairs, can't merge
		if (pairs.empty()){
			return -1E300;
		}
		int num_pairs = pairs.size();
		int k = DNest3::randInt(num_pairs);
		int i = pairs[k][0];
		int j = pairs[k][1];

		// make up merged component
		std::vector<double> newcomponent(this->num_dimensions, 0.);
		newcomponent[0] = (this->components[i][0]*this->components[i][2+midbin] + this->components[j][0]*this->components[j][2+midbin]) /
			(this->components[i][2+midbin] + this->components[j][2+midbin]);
		newcomponent[1] = (this->components[i][1]*this->components[i][2+midbin] + this->components[j][1]*this->components[j][2+midbin]) /
			(this->components[i][2+midbin] + this->components[j][2+midbin]);
		for (int l=2; l<this->num_dimensions; l++){
			newcomponent[l] = this->components[i][l] + this->components[j][l];
		}

		double newprior = this->dist.log_pdf(newcomponent);
		if (newprior > -1E300){
			this->num_components--;
			double oldprior = this->dist.log_pdf(this->components[i]) + this->dist.log_pdf(this->components[j]);
			logH += newprior - oldprior;
			logH -= log(cos(this->components[i][1])) + log(cos(this->components[j][1])) - log(cos(newcomponent[1]));
			logH -= log(2*maxkickl) + log(2*maxkickb) - log(this->dist.angular_area());
			logH -= log(this->num_components) + log(this->num_components+1);
			logH += log(2*num_pairs); //multiply by 2 to get ordered pairs
			for (int l=2; l<this->num_dimensions; l++){
				logH -= log(newcomponent[l]);
			}

			this->removed.push_back(this->components[j]);
			this->components.erase(this->components.begin()+j);
			this->u_components.erase(this->u_components.begin()+j);
			// since j>i, erasing j does not move i

			this->removed.push_back(this->components[i]);
			this->added.push_back(newcomponent);
			this->components[i].assign(newcomponent.begin(), newcomponent.end());
			this->dist.to_uniform(newcomponent);
			this->u_components[i].assign(newcomponent.begin(), newcomponent.end());

			return logH;
		}
		else{
			return -1E300;
		}
	}
}

// returns unordered pairs
template<class Distribution>
std::vector< std::vector<int> > MyRJObject<Distribution>::find_pairs(double maxkickl, double maxkickb){
	std::vector< std::vector<int> > pairs;

	for (int i=0; i<this->num_components; i++){
		for (int j=i+1; j<this->num_components; j++){
			double dist_l = this->components[j][0] - this->components[i][0];
			double dist_b = this->components[j][1] - this->components[i][1];
			if ((dist_l > -2*maxkickl) && (dist_l < 2*maxkickl) && (dist_b > -2*maxkickb) && (dist_b < 2*maxkickb)){
				std::vector<int> newpair(2);
				newpair[0] = i;
				newpair[1] = j;
				pairs.push_back(newpair);
			}
		}
	}

	return pairs;
}

#endif
