#ifndef PORTILLO_PCAT_DNEST_MYMODELGLOBALS_H
#define PORTILLO_PCAT_DNEST_MYMODELGLOBALS_H


#include <vector>

struct MyModelGlobals {
    public:
        // isotropic background
        double bg_min, bg_max; // same prior bounds for all energy bins

        std::vector<double> data;
        std::vector<double> exposure;
        std::vector<double> etemplate; //actual templates

        // emission templates
        std::vector<double> tem_min;
        std::vector<double> tem_max;

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


#endif //PORTILLO_PCAT_DNEST_MYMODELGLOBALS_H
