#ifndef _SloanData_
#define _SloanData_

#include <vector>

class SloanData
{
	private:
		// Number of pixels
		int width, height;
		// Number of bands
		int nband;
		// bounds on longitude and latitude
		double xmin, xmax, ymin, ymax;
		double bias, gain; // Sloan 'gain' = inverse gain

		// Coordinates of pixel centers
		std::vector< double > x_rays;
		std::vector< double > y_rays;

		// PSFs
		int psf_size;
		std::vector< double > psfs;

		// The pixels
		std::vector< double > image;

	public:
		SloanData();
		void load(const char* image_file, const char* psf_file, const char* pixel_file);

		// Getters
		int get_width() const { return width; }
		int get_height() const { return height; }
		int get_nband() const { return nband; }
		double get_xmin() const { return xmin; }
		double get_xmax() const { return ymax; }
		double get_ymin() const { return xmin; }
		double get_ymax() const { return ymax; }
		int get_psf_size() const { return psf_size; }
		double get_bias() const { return bias; }
		double get_gain() const { return gain; }
		const std::vector< double >& get_x_rays() const
			{ return x_rays; }
		const std::vector< double >& get_y_rays() const
			{ return y_rays; }
		const std::vector< double >& get_image() const
			{ return image; }
		const std::vector< double >& get_psfs() const
			{ return psfs; }

	// Singleton
	private:
		static SloanData instance;
	public:
		static SloanData& get_instance() { return instance; }
};

#endif

