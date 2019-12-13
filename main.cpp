#include <iostream>
#include <CommandLineOptions.h>
#include <MTSampler.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "FermiModelOptions.h"
#include "MyOptions.h"
#include "FermiModel.h"
#include "FermiData.h"
#include "SloanModelOptions.h"
#include "SloanData.h"
#include "SloanModel.h"
#include "buildinfo.h"

using namespace std;
using namespace DNest3;
using boost::property_tree::ptree;

constexpr char BUILDINFO_CMDLINE_FLAG[] = "--info";

int main(int argc, char** argv)
{
    if (argc == 2 && std::strcmp(argv[1], BUILDINFO_CMDLINE_FLAG) == 0) {
        // If only command line arg is --info, just print buildinfo.
        std::printf(BUILDINFO);
        std::printf("### SRC LAST COMMIT\n%s###", SRC_LAST_COMMIT);
        std::exit(0);
    }
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

	ptree pt;
	read_xml("run-"+options.get_configFile()+"/OPTIONS-MODEL", pt);
	
	// TODO look at Start file?
	// must load data, then model options, before creating sampler
	// then load levels file if requested
	string model = pt.get<string>("modeloptions.model");
	if (model == "fermi"){
        std::cout<<"ERROR: Fermimodel has been disabled for testing. Uncomment and fix vector allocators to continue."<<std::endl;
        std::exit(1);
//		FermiData::get_instance().load(("Data/"+options.get_dataFile()+"_cts.txt").c_str(),
//			("Data/"+options.get_dataFile()+"_exp.txt").c_str(),
//			("Data/"+options.get_dataFile()+"_pix.txt").c_str(),
//			("Data/"+options.get_dataFile()+"_tem.txt").c_str());
//		FermiModelOptions::get_instance().load(pt);
//		MTSampler<FermiModel> sampler(options.get_numThreads(), options.get_compression_double(), samplerOptions);
//		if(options.get_levelsFile().compare("") != 0)
//			sampler.loadLevels(options.get_levelsFile().c_str());
//		sampler.run();
	}
	if (model == "sloan"){
		SloanData::get_instance().load(("Data/"+options.get_dataFile()+"_cts.txt").c_str(),
			("Data/"+options.get_dataFile()+"_psf.txt").c_str(),
			("Data/"+options.get_dataFile()+"_pix.txt").c_str());
		SloanModelOptions::get_instance().load(pt);
		MTSampler<SloanModel> sampler(options.get_numThreads(), options.get_compression_double(), samplerOptions);
		if(options.get_levelsFile().compare("") != 0)
			sampler.loadLevels(options.get_levelsFile().c_str());
		sampler.run();
	}

	return 0;
}

