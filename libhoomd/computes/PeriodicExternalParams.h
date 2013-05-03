#ifndef _PERIODIC_EXTERNAL_PARAMS_H_
#define _PERIODIC_EXTERNAL_PARAMS_H_

#include "HOOMDMath.h"

struct PeriodicExternalParams 
       {
       Scalar order_parameter;
       int3 lattice_vector_1;
       int3 lattice_vector_2;
       int3 lattice_vector_3;
       Scalar interface_width;
       int    periodicity;
      
       PeriodicExternalParams(Scalar orderParameter, int3 latticeVector1, int3 latticeVector2, int3 latticeVector3, Scalar interfaceWidth, int Periodicity)
            {
            order_parameter = orderParameter;
            lattice_vector_1 = latticeVector1;
            lattice_vector_2 = latticeVector2;
            lattice_vector_3 = latticeVector3;
            interface_width = interfaceWidth;
            periodicity = Periodicity;
            }

       PeriodicExternalParams()
            {}
       };
#endif

