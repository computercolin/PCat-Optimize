#include "MyOptions.h"
#include <string>


MyOptions::MyOptions(const char* outname, bool useGzip)
:Options(("run-"+std::string(outname)+"/OPTIONS").c_str(), useGzip)
{
	// this doesn't work outside the braces
	// I guess we have to explicity state that we change the protected members
	sampleFile = "run-"+std::string(outname)+"/sample.txt";
	sampleInfoFile = "run-"+std::string(outname)+"/sample_info.txt";
	levelsFile = "run-"+std::string(outname)+"/levels.txt";
}
