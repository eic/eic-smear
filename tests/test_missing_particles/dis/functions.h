#ifndef EJANA_FUNCTIONS_H
#define EJANA_FUNCTIONS_H

#include <cmath>

namespace dis {

    /**
    Helper function to convert eta to theta (radians)

    Detector acceptances require theta, not eta
    */
    inline double eta_to_theta(double eta) {
        return 2. * std::atan(std::exp(-eta));
    }

}



#endif //EJANA_FUNCTIONS_H
