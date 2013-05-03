#ifndef _LOCAL_EXTERNAL_PARAMS_H_
#define _LOCAL_EXTERNAL_PARAMS_H_

#include "HOOMDMath.h"

struct LocalExternalParams 
       {
       unsigned int perp_index;
       unsigned int parallel_index;
       Scalar fraction_length;
       Scalar order_parameter;
       Scalar interface_width;
       int    periodicity;
      
       LocalExternalParams(unsigned int perpIndex, unsigned int parallelIndex, Scalar fractionLength, Scalar orderParameter, Scalar interfaceWidth, int Periodicity)
            {
            perp_index = perpIndex;
            parallel_index = parallelIndex;
            fraction_length = fractionLength;
            order_parameter = orderParameter;
            interface_width = interfaceWidth;
            periodicity = Periodicity;
            }

       LocalExternalParams()
            {}

       };

#endif

