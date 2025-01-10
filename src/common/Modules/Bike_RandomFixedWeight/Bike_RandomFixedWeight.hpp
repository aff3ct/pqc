#ifndef BIKE_RANDOMFIXEDWEIGHT_H
#define BIKE_RANDOMFIXEDWEIGHT_H

#include <streampu.hpp>

#include "../../Tools/tools.hpp"


namespace spu
{
  namespace module
  {

      class Bike_RandomFixedWeight  : public Stateful {      
  
      private:
    
	  int frame_size;
	  int weight;

      public:

	  Bike_RandomFixedWeight(int frame_size, int weight);
	  virtual ~Bike_RandomFixedWeight();

      protected:

	  virtual void random_fixed_weight(int* input, int* output, const int frame_id);

    }; 
  }
}

#endif // BIKE_RANDOM_FIXED_WEIGHT_H

