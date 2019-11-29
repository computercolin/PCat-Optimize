#ifndef PORTILLO_PCAT_DNEST_SLOANMODELGLOBALS_H
#define PORTILLO_PCAT_DNEST_SLOANMODELGLOBALS_H

#include <mutex>
#include "MyModelGlobals.h"


/** NOTE: currently unused.
 *  SloanModel only has a few primitive values that are shared, constant, and could be shared globally.
 */



struct SloanModelGlobals {
    public:
        double width;
        double height;
        static SloanModelGlobals& get_instance() {
            static SloanModelGlobals instance;
            return instance;
        }
        SloanModelGlobals(SloanModelGlobals const&) = delete;
        void operator=(SloanModelGlobals const&) = delete;

    private:
        SloanModelGlobals() {}
};

#endif //PORTILLO_PCAT_DNEST_SLOANMODELGLOBALS_H
