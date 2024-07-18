#include "CM_Decoder.hpp"
#include <iostream>

using namespace spu;
using namespace spu::module;
using namespace std;

CM_Decoder:: CM_Decoder(int frame_size, int synd_size, int weight, CM_secret_key& SK) :
    Module(),
    frame_size(frame_size),
    synd_size(synd_size),
    weight(weight) {

    this->set_name("CM_Decoder");
    this->set_short_name("CM_Decoder");

    auto &t = create_task("cm_decoder");
    auto input   = create_socket_in<int>(t, "input", synd_size);
    auto output  = create_socket_out<int>(t, "output", frame_size);

    
    this->create_codelet(t, [input, output, &SK](Module &m, runtime::Task &t,
						 const size_t frame_id) -> int {
	static_cast<CM_Decoder&>(m).cm_decoder(static_cast<int*>(t[input].get_dataptr()),
					       static_cast<int*>(t[output].get_dataptr()),
					       static_cast<CM_secret_key&>(SK),
					       frame_id);
	return 0;
    }
	);
}


CM_Decoder:: ~CM_Decoder() {
}


void
CM_Decoder:: cm_decoder(int* input, int* output, const CM_secret_key& SK, const int frame_id) {
    fq_ctx_t* ctx = SK.get_ctx(); /* finite field F_2 */

    /* temp vectors for conversion between F_2 and int values */
    fq_struct* tmp_s = _fq_vec_init(this->synd_size, *ctx);
    fq_struct* tmp_e = _fq_vec_init(this->frame_size, *ctx);

    
    /* conversion int to F_2 */
    _int_vec_2_fq(tmp_s, input, this->synd_size, *ctx);

    
    /* decoding */
    int b = CM_syndrome_decoding(tmp_e, tmp_s, SK.get_alpha(), SK.g, this->frame_size, this->weight,
				 *ctx);
    
    // CM_encoding(tmp_s, tmp_e, PK.T, this->frame_size, *ctx_q);

    /* reverse conversion F_2 to int */
    _fq_vec_2_int(output, tmp_e, this->frame_size, *ctx);
    
    /* clear memory */
    _fq_vec_clear(tmp_e, this->frame_size, *ctx);
    _fq_vec_clear(tmp_s, this->synd_size, *ctx);
}
