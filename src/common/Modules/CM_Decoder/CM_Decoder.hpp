#ifndef CM_DECODER_H
#define CM_DECODER_H

#include "../../../../streampu/include/streampu.hpp"
#include "../../Tools/tools.hpp"
#include "../../Tools/codes.hpp"
#include "../../Tools/ClassicMcEliece/CM_secret_key.hpp"
#include "../../Tools/ClassicMcEliece/CM_public_key.hpp"



namespace spu
{
    namespace module
    {

	class CM_Decoder  : public Stateful {      

      
	private:
    
	    int frame_size;
	    int synd_size;
	    int weight;		// hamming weight of the error ?

	  
	public:
      
	    CM_Decoder(int frame_size, int synd_size, int weight, CM_secret_key& SK);
	    virtual ~CM_Decoder();

      
	protected:

	    virtual void cm_decoder(int* input, int* output, const CM_secret_key& SK,
				    const int frame_id);

	}; 
    }
}

#endif // CM_DECODER_H

