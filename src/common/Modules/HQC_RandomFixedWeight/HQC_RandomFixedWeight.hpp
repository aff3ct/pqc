#ifndef HQC_RANDOMFIXEDWEIGHT_H
#define HQC_RANDOMFIXEDWEIGHT_H

#include "streampu.hpp"

#include "../../Tools/tools.hpp"


namespace spu
{
  namespace module
  {

      class HQC_RandomFixedWeight  : public Module {      
  
      private:
    
	  int frame_size;
	  int weight;

      public:

	  HQC_RandomFixedWeight(int frame_size, int weight);
	  virtual ~HQC_RandomFixedWeight();

      protected:

	  virtual void random_fixed_weight(int* input, int* output, const int frame_id);

    }; 
  }
}

#endif // HQC_RANDOM_FIXED_WEIGHT_H

