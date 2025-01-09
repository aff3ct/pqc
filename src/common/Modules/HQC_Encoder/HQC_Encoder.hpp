#ifndef HQC_ENCODER_H
#define HQC_ENCODER_H

#include "streampu.hpp"

#include "../../Tools/tools.hpp"
#include "../../Tools/codes.hpp"
#include "../../Tools/HQC/HQC_secret_key.hpp"
#include "../../Tools/HQC/HQC_public_key.hpp"



namespace spu
{
    namespace module
    {

	class HQC_Encoder  : public Stateful {      

      
	private:
	    
	    int frame_size;    	/* it is k in HQC */
	    int output_size1;	/* it is n in HQC */
	    int output_size2;	/* it is n1 * n2 in HQC */

	    int n1;		/* length of the RS code */
	    int r;		/* RM code is duplicated `r` times */
	    int w;		/* weight of e, r1, r2 */
	    

	    
	public:
	    
	    HQC_Encoder(HQC_public_key& PK, int k, int len, int r, int w);
	    virtual ~HQC_Encoder();

      
	protected:

	    virtual void hqc_encoder(int* input, int* output1, int* output2, const HQC_public_key& PK,
				     const int frame_id);

	}; 
    }
}

#endif // HQC_ENCODER_H

