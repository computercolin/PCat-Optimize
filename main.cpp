#include <iostream>
#include <CommandLineOptions.h>
#include <MTSampler.h>
#include "ModelOptions.h"
#include "MyOptions.h"
#include "MyModel.h"
#include "Data.h"

using namespace std;
using namespace DNest3;

int main(int argc, char** argv)
{
	//TODO clean up main using DNest3 Start.h
	CommandLineOptions options(argc, argv);

	std::cout<<"# Using "<<options.get_numThreads()<<" thread"<<
		((options.get_numThreads() == 1)?("."):("s."))<<std::endl;

	std::cout<<"# Target compression factor between levels = ";
	std::cout<<options.get_compression()<<std::endl;

	// Seed random number generator
	std::cout<<"# Seeding random number generator with "<<
		options.get_seed_long()<<"."<<std::endl;
	RandomNumberGenerator::initialise_instance();
	RandomNumberGenerator::get_instance().set_seed(options.get_seed_long());

	// Load sampler options from file
	// Added option to specify output file name
	MyOptions samplerOptions(options.get_configFile().c_str(), options.set_gzip());

	// Load data; must be done before sampler created
        Data::get_instance().load(("Data/"+options.get_dataFile()+"_cts.txt").c_str(),
                                ("Data/"+options.get_dataFile()+"_exp.txt").c_str(),
				("Data/"+options.get_dataFile()+"_pix.txt").c_str(),
				("Data/"+options.get_dataFile()+"_tem.txt").c_str());

	// Load model options from file
	ModelOptions::get_instance().load(("run-"+options.get_configFile()+"/OPTIONS-MODEL").c_str());

	// Create sampler
	MTSampler<MyModel> sampler(options.get_numThreads(),
					options.get_compression_double(),
					samplerOptions);

	// Load levels file if requested
	if(options.get_levelsFile().compare("") != 0)
		sampler.loadLevels(options.get_levelsFile().c_str());

	sampler.run();

	return 0;
}

