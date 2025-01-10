#ifndef BIKE_DECODER_H
#define BIKE_DECODER_H

#include <streampu.hpp>

#include "../../Tools/tools.hpp"
#include "../../Tools/codes.hpp"
#include "../../Tools/Bike/Bike_secret_key.hpp"
#include "../../Tools/Bike/Bike_public_key.hpp"



namespace spu
{
    namespace module
    {

	class Bike_Decoder  : public Stateful {      
	    
      
	private:
    
	    int frame_size;    	/* it is 2*r in Bike : output_size */
	    int input_size;	/* it is r in Bike */
    
	    int weight;		/* hamming weight of the secret key ? */
	    int NbIter;	      	/* nb of iterations in BGF algo. */
	    int tau;		/* tau for BGF algo. */
	    
	public:
      
	    Bike_Decoder(int weight, int NbIter, int tau, Bike_secret_key& SK);
	    virtual ~Bike_Decoder();

      
	protected:

	    virtual void bike_decoder(int* input, int* output,
				      const Bike_secret_key& SK, const int frame_id);
	}; 
    }
}

#endif // BIKE_DECODER_H

