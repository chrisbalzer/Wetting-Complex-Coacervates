#include "../../SCFT.h"

void SCFT::mixDensityPicardRelax() {
    int max_index = M*(numComponents + 1);
    double temp;

    // Update component density
    for (int i = 0; i < max_index; i++) {
        if(i < max_index - M) {
            temp = (1-mixF)*X[i] + mixF*(X[i] + F[i]);
            if(temp > 1) temp = 1;
            if(temp < 0) temp = 0.9*X[i];
            X[i] = temp;
        }
        else {
            temp = (1-0.25*mixF)*X[i] + 0.25*mixF*(X[i] + F[i]);
            if(fabs(temp) < 1E-20) temp = 0;
            X[i] = temp;
        }
    }
}
