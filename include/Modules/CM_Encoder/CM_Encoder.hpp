#ifndef CM_ENCODER_H
#define CM_ENCODER_H

#include <streampu.hpp>

#include "../../Tools/tools.hpp"
#include "../../Tools/codes.hpp"
#include "../../Tools/ClassicMcEliece/CM_secret_key.hpp"
#include "../../Tools/ClassicMcEliece/CM_public_key.hpp"



namespace spu
{
    namespace module
    {

	class CM_Encoder  : public Stateful {      

      
	private:
    
	    int frame_size;
	    int out_size;
      
	public:
      
	    CM_Encoder(int frame_size, int out_size, CM_public_key& PK);
	    virtual ~CM_Encoder();

      
	protected:

	    virtual void cm_encoder(int* input, int* output, const CM_public_key& PK,
				    const int frame_id);

	}; 
    }
}

#endif // CM_ENCODER_H

