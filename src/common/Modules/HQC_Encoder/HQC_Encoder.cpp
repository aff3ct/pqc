#include "HQC_Encoder.hpp"
#include <iostream>

using namespace spu;
using namespace spu::module;
using namespace std;

HQC_Encoder:: HQC_Encoder(HQC_public_key& PK,  int k,  int len, int r, int w) :
    Module(),
    frame_size(k),
    output_size1(PK.get_n()),
    output_size2(len),
    n1(PK.get_n1()),
    r(r),
    w(w) {

    this->set_name("HQC_Encoder");
    this->set_short_name("HQC_Encoder");

    auto &t = create_task("hqc_encoder");
    auto input   = create_socket_in<int>(t, "input", frame_size);
    auto output1  = create_socket_out<int>(t, "output1", output_size1);
    auto output2  = create_socket_out<int>(t, "output2", output_size2);

    this->create_codelet(t, [input, output1, output2, &PK](Module &m, runtime::Task &t,
						 const size_t frame_id) -> int {
	static_cast<HQC_Encoder&>(m).hqc_encoder(static_cast<int*>(t[input].get_dataptr()),
						 static_cast<int*>(t[output1].get_dataptr()),
						 static_cast<int*>(t[output2].get_dataptr()),
						 static_cast<HQC_public_key&>(PK),
						 frame_id);
	return 0;
    }			 	);
}


HQC_Encoder:: ~HQC_Encoder() {
}


void
HQC_Encoder:: hqc_encoder(int* input, int* output1, int* output2, const HQC_public_key& PK,
			  const int frame_id) {

    fq_ctx_t* ctx_q = PK.get_ctx_q(); /* finite field F_2 */
    fq_ctx_t* ctx = PK.get_ctx(); /* finite field F_2^m */

    /* temporary poly for computation */
    fq_poly_t tmp_pol; fq_poly_init(tmp_pol, *ctx_q);
    fq_poly_zero(tmp_pol, *ctx_q);
    
    /* temporary vectors for conversion between F_2 and int values */
    fq_struct* tmp_vec = _fq_vec_init(this->frame_size, *ctx_q);
    _fq_vec_zero(tmp_vec, this->frame_size, *ctx_q);

    fq_struct* tmp_vec1 = _fq_vec_init(this->output_size1, *ctx_q);
    _fq_vec_zero(tmp_vec1, this->output_size1, *ctx_q);

    fq_struct* tmp_vec2 = _fq_vec_init(this->output_size2, *ctx_q);
    _fq_vec_zero(tmp_vec2, this->output_size2, *ctx_q);

    
    
    /* put input into a polynomial */
    _int_vec_2_fq(tmp_vec, input, this->frame_size, *ctx_q);
    fq_poly_set_coeffs(tmp_pol, tmp_vec, this->frame_size, *ctx_q);

    
    fq_poly_t res1; fq_poly_init(res1, *ctx_q);
    fq_poly_t res2; fq_poly_init(res2, *ctx_q); 


    printf("just before encoding \n");
    
    /* encoding */
    HQC_encoding(res1, res2, tmp_pol, PK.h, PK.s, PK.alpha, this->output_size1,
		 this->output_size2, this->n1, this->r, this->w, this->w, *ctx,
		 *ctx_q);

    printf("just after encoding \n");
	
    /* put it in vec format */
    _fq_vec_zero(tmp_vec1, this->output_size1, *ctx_q);
    fq_poly_get_coeffs(tmp_vec1, res1, this->output_size1, *ctx_q);
    fq_poly_get_coeffs(tmp_vec2, res2, this->output_size2, *ctx_q);

    printf("just after vec format \n");
    
    /* reverse conversion F_2 to int */
    _fq_vec_2_int(output1, tmp_vec1, this->output_size1, *ctx_q);
    _fq_vec_2_int(output2, tmp_vec2, this->output_size2, *ctx_q);

    printf("just after conversion \n");

    /* clear memory */
    /* _fq_vec_clear(tmp_e, this->frame_size, *ctx_q); */
    _fq_vec_clear(tmp_vec, this->frame_size, *ctx_q);
    _fq_vec_clear(tmp_vec1, this->output_size1, *ctx_q);
    _fq_vec_clear(tmp_vec2, this->output_size2, *ctx_q);
    fq_poly_clear(tmp_pol, *ctx_q);
}
