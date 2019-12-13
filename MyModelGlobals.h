#ifndef PORTILLO_PCAT_DNesT_MYMODELGLOBALS_H
#define PORTILLO_PCAT_DNesT_MYMODELGLOBALS_H


#include <vector>
#include "Utils.h"


struct MyModelGlobals {
    public:
        // isotropic background
        double bg_min, bg_max; // same prior bounds for all energy bins

        DNest3::vec_align32<double> data;
        DNest3::vec_align32<double> exposure;
        DNest3::vec_align32<double> etemplate; //actual templates

        // emission templates
        DNest3::vec_align32<double> tem_min;
        DNest3::vec_align32<double> tem_max;

        static MyModelGlobals& get_instance() {
            static MyModelGlobals instance;
            return instance;
        }

        MyModelGlobals(MyModelGlobals const&) = delete;
        void operator=(MyModelGlobals const&) = delete;

    private:
         MyModelGlobals() {}

//        MyModelGlobals(MyRJObject<MyDistribution> objects, int nbin, int npsf, int npix, double pixel_area,
//                vector<double> data, vector<double> exposure, double bg_min, double bg_max,
//                int ntem, vector<double> tem_min, vector<double> tem_max, vector<double> etemplate);
};


#endif //PORTILLO_PCAT_DNesT_MYMODELGLOBALS_H
