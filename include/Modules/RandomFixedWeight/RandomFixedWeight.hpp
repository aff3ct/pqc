#ifndef RANDOMFIXEDWEIGHT_H
#define RANDOMFIXEDWEIGHT_H

#include <streampu.hpp>

#include "../../Tools/tools.hpp"


namespace spu
{
  namespace module
  {

    class RandomFixedWeight  : public Stateful {      
  
    private:
    
      int frame_size;
      int weight;
      int N;
      int tau;

    public:

      RandomFixedWeight(int frame_size, int weight, int N, int tau);
      virtual ~RandomFixedWeight();

    protected:

      virtual void random_fixed_weight(int* input, int* output, const int frame_id);

    }; 
  }
}

#endif // RANDOM_FIXED_WEIGHT_H

