#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

Data Data::instance;

Data::Data()
{

}

void Data::load(const char* counts_file, const char* exposure_file, const char* pixel_file, const char* etemplate_file)
{
	fstream fin(pixel_file, ios::in);
	// pixel file
        if (!fin)
                cerr<<"# ERROR: couldn't open file "<<pixel_file<<"."<<endl;
	fin>>npix>>nbin>>npsf; // read number of pixels, bins, and PSF classes
	fin>>lmin>>lmax>>bmin>>bmax>>pixel_area; // read l, b bounds, pixel area
        l_rays.assign(npix, 0);
        b_rays.assign(npix, 0);
        double junk;
        for (int i=0; i<npix; i++){
                fin>>l_rays[i]>>b_rays[i]>>junk; // read l, b, throw away pixel number
		vector< double > unit_vector(3);
		unit_vector[0] = cos(b_rays[i]) * cos(l_rays[i]);
		unit_vector[1] = cos(b_rays[i]) * sin(l_rays[i]);
		unit_vector[2] = sin(b_rays[i]);
		unit_vectors.push_back(unit_vector);
	}
        fin.close();

	// counts file
	fin.open(counts_file, ios::in);
	if (!fin)
		cerr<<"# ERROR: couldn't open file "<<counts_file<<"."<<endl;
	image.assign(nbin*npsf*npix, 0);
	for (int i=0; i<nbin*npsf*npix; i++)
		fin>>image[i];
	fin.close();

	// exposure file
	fin.open(exposure_file, ios::in);
	if (!fin)
		cerr<<"# ERROR: couldn't open file "<<exposure_file<<"."<<endl;
	exposure.assign(nbin*npsf*npix, 0);
	for (int i=0; i<nbin*npsf*npix; i++)
		fin>>exposure[i];
	fin.close();

	// template file
	fin.open(etemplate_file, ios::in);
	if (!fin)
		cerr<<"# ERROR: couldn't open file"<<etemplate_file<<"."<<endl;
	fin>>ntem;
	etemplate.assign(ntem*nbin*npsf*npix, 0);
	for (int i=0; i<ntem*nbin*npsf*npix; i++)
		fin>>etemplate[i];
	fin.close();
}
