#ifndef BIKE_ENCODER_H
#define BIKE_ENCODER_H

#include "streampu.hpp"

#include "../../Tools/tools.hpp"
#include "../../Tools/codes.hpp"
#include "../../Tools/Bike/Bike_secret_key.hpp"
#include "../../Tools/Bike/Bike_public_key.hpp"



namespace spu
{
    namespace module
    {

	class Bike_Encoder  : public Module {      

      
	private:
    
	    int frame_size;    	/* it is 2*r in Bike */
	    int output_size;	/* it is r in Bike */
      
	public:
      
	    Bike_Encoder(Bike_public_key& PK);
	    virtual ~Bike_Encoder();

      
	protected:

	    virtual void bike_encoder(int* input, int* output, const Bike_public_key& PK,				      const int frame_id);

	}; 
    }
}

#endif // BIKE_ENCODER_H

