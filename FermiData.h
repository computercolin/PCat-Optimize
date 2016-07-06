#ifndef _FermiData_
#define _FermiData_

#include <vector>

class FermiData
{
	private:
		// Number of pixels
		int npix;
		// Number of energy bins
		int nbin;
		// Number of PSF classes
		int npsf;
		// Number of templates
		int ntem;

		// bounds on longitude and latitude
		double lmin, lmax, bmin, bmax;
		// pixel area, in steradians
		double pixel_area;

		// Coordinates of pixel centers
		std::vector< double > l_rays;
		std::vector< double > b_rays;
		std::vector< std::vector< double > > unit_vectors;

		// The pixels
		std::vector< double > image;
		std::vector< double > exposure;
		std::vector< double > etemplate;

	public:
		FermiData();
		void load(const char* image_file, const char* exposure_file, const char* pixel_file, const char* etemplate_file);

		// Getters
		int get_npix() const { return npix; }
		int get_nbin() const { return nbin; }
		int get_npsf() const { return npsf; }
		int get_ntem() const { return ntem; }
		double get_lmin() const { return lmin; }
		double get_lmax() const { return lmax; }
		double get_bmin() const { return bmin; }
		double get_bmax() const { return bmax; }
		double get_pixel_area() const { return pixel_area; }
		const std::vector< double >& get_l_rays() const
			{ return l_rays; }
		const std::vector< double >& get_b_rays() const
			{ return b_rays; }
		const std::vector< std::vector< double > >& get_unit_vectors() const
			{ return unit_vectors; }
		const std::vector< double >& get_image() const
			{ return image; }
		const std::vector< double >& get_exposure() const
			{ return exposure; }
		const std::vector< double >& get_etemplate() const
			{ return etemplate; }

	// Singleton
	private:
		static FermiData instance;
	public:
		static FermiData& get_instance() { return instance; }
};

#endif

