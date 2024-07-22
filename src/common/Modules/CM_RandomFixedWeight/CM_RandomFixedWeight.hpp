#ifndef CM_RANDOMFIXEDWEIGHT_H
#define CM_RANDOMFIXEDWEIGHT_H

#include "streampu.hpp"

#include "../../Tools/tools.hpp"


namespace spu
{
  namespace module
  {

    class CM_RandomFixedWeight  : public Module {      
  
    private:
    
      int frame_size;
      int weight;
      int N;
      int tau;

    public:

      CM_RandomFixedWeight(int frame_size, int weight, int N, int tau);
      virtual ~CM_RandomFixedWeight();

    protected:

      virtual void random_fixed_weight(int* input, int* output, const int frame_id);

    }; 
  }
}

#endif // CM_RANDOM_FIXED_WEIGHT_H

