#ifndef HQC_DECODER_H
#define HQC_DECODER_H

#include "streampu.hpp"

#include "../../Tools/tools.hpp"
#include "../../Tools/codes.hpp"
#include "../../Tools/HQC/HQC_secret_key.hpp"
#include "../../Tools/HQC/HQC_public_key.hpp"



namespace spu
{
    namespace module
    {

	class HQC_Decoder  : public Stateful {      

      
	private:
	    
	    int input_size1;  	/* it is n in HQC */
	    int input_size2;	/* it is n1 * n2 in HQC */
	    int output_size;	/* it is k in HQC */

	    int r;		/* RM code is duplicated `r` times */
	    

	    
	public:
	    
	    HQC_Decoder(HQC_secret_key& SK, HQC_public_key& PK, int k, int len,
			int r);
	    virtual ~HQC_Decoder();

      
	protected:

	    virtual void hqc_decoder(int* input1, int* input2, int* output, const HQC_secret_key& SK,
				     const HQC_public_key& PK, const int frame_id);

	}; 
    }
}

#endif // HQC_DECODER_H

