#include "SloanData.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

SloanData SloanData::instance;

SloanData::SloanData()
{

}

void SloanData::load(const char* counts_file, const char* psf_file, const char* pixel_file)
{
	fstream fin(pixel_file, ios::in);
	// pixel file
        if (!fin)
                cerr<<"# ERROR: couldn't open file "<<pixel_file<<"."<<endl;
	fin>>width>>height>>nband; // read number of pixels, bins, and PSF classes
	fin>>xmin>>xmax>>ymin>>ymax; //read x,y bounds
	fin>>bias>>gain; //read bias and gain
	// FIXME assuming aligned images
        /*x_rays.assign(band*height*width, 0);
        y_rays.assign(band*height*width, 0);
        for (int i=0; i<band*height*width; i++){
                fin>>x_rays[i]>>y_rays[i];
	}*/
        fin.close();

	// PSF file
	fin.open(psf_file, ios::in);
	if (!fin)
		cerr<<"# ERROR: couldn't open file "<<psf_file<<"."<<endl;
	fin>>psf_size>>psf_resampling;
	psfs.assign(nband*psf_size*psf_size, 0);
	for (int i=0; i<nband*psf_size*psf_size; i++){
		fin>>psfs[i];
	}
	fin.close();

	// counts file
	fin.open(counts_file, ios::in);
	if (!fin)
		cerr<<"# ERROR: couldn't open file "<<counts_file<<"."<<endl;
	image.assign(nband*height*width, 0);
	for (int i=0; i<nband*height*width; i++)
		fin>>image[i];
	fin.close();
}
