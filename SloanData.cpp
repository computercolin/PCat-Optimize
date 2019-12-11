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
	fin>>bias>>gain>>exposure; //read bias and gain
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
	int maxsize = nband*psf_size*psf_size*psf_resampling*psf_resampling;
	psfs.assign(maxsize, 0);
	for (int i=0; i<maxsize; i++){
		fin>>psfs[i];
	}
	fin.close();

	// counts file
	fin.open(counts_file, ios::in);
	if (!fin)
		cerr<<"# ERROR: couldn't open file "<<counts_file<<"."<<endl;
	int imgMaxsize = nband*height*width;
	image.assign(imgMaxsize, 0);
	for (int i=0; i<imgMaxsize; i++)
		fin>>image[i];
	fin.close();
}
