#include "CM_Encoder.hpp"
#include <iostream>

using namespace spu;
using namespace spu::module;
using namespace std;

CM_Encoder:: CM_Encoder(int frame_size, int out_size, CM_public_key& PK) :
    Module(),
    frame_size(frame_size),
    out_size(out_size) {

    this->set_name("CM_Encoder");
    this->set_short_name("CM_Encoder");

    auto &t = create_task("cm_encoder");
    auto input   = create_socket_in<int>(t, "input", frame_size);
    auto output  = create_socket_out<int>(t, "output", out_size);


    
    this->create_codelet(t, [input, output, &PK](Module &m, runtime::Task &t,
						const size_t frame_id) -> int {
	static_cast<CM_Encoder&>(m).cm_encoder(static_cast<int*>(t[input].get_dataptr()),
					       static_cast<int*>(t[output].get_dataptr()),
					       static_cast<CM_public_key&>(PK),
					       frame_id);
	return 0;
    }			 	);

}


CM_Encoder:: ~CM_Encoder() {
}


void
CM_Encoder:: cm_encoder(int* input, int* output, const CM_public_key& PK, const int frame_id) {
    fq_ctx_t* ctx_q = PK.get_ctx_q(); /* finite field F_2 */

    /* temp vectors for conversion between F_2 and int values */
    fq_struct* tmp_e = _fq_vec_init(this->frame_size, *ctx_q);
    fq_struct* tmp_s = _fq_vec_init(this->out_size, *ctx_q);

    
    /* conversion int to F_2 */
    _int_vec_2_fq(tmp_e, input, this->frame_size, *ctx_q);

    /* encoding */
    CM_encoding(tmp_s, tmp_e, PK.T, this->frame_size, *ctx_q);

    /* reverse conversion F_2 to int */
    _fq_vec_2_int(output, tmp_s, this->out_size, *ctx_q);
    
    /* clear memory */
    _fq_vec_clear(tmp_e, this->frame_size, *ctx_q);
    _fq_vec_clear(tmp_s, this->out_size, *ctx_q);
}
