#include "HQC_Decoder.hpp"
#include <iostream>

using namespace spu;
using namespace spu::module;
using namespace std;

HQC_Decoder:: HQC_Decoder(HQC_secret_key& SK, HQC_public_key& PK, int k, int len,
			  int r) :
    Module(),
    input_size1(PK.get_n()),
    input_size2(len),
    output_size(k),
    r(r) {

    this->set_name("HQC_Decoder");
    this->set_short_name("HQC_Decoder");

    auto &t = create_task("hqc_decoder");
    auto input1   = create_socket_in<int>(t, "input1", input_size1);
    auto input2  = create_socket_in<int>(t, "input2", input_size2);
    auto output  = create_socket_out<int>(t, "output", output_size);

    this->create_codelet(t, [input1, input2, output, &PK, &SK](Module &m, runtime::Task &t,
								const size_t frame_id) -> int {
	static_cast<HQC_Decoder&>(m).hqc_decoder(static_cast<int*>(t[input1].get_dataptr()),
						 static_cast<int*>(t[input2].get_dataptr()),
						 static_cast<int*>(t[output].get_dataptr()),
						 static_cast<HQC_secret_key&>(SK),
						 static_cast<HQC_public_key&>(PK),
						 frame_id);
	return 0;
    }			 	);
}


HQC_Decoder:: ~HQC_Decoder() {
}


void
HQC_Decoder:: hqc_decoder(int* input1, int* input2, int* output, const HQC_secret_key& SK,
			  const HQC_public_key& PK, const int frame_id) {
  
    fq_ctx_t* ctx_q = PK.get_ctx_q(); /* finite field F_2 */
    fq_ctx_t* ctx = PK.get_ctx(); /* finite field F_2^m */

    /* temporary poly for computation */
    fq_poly_t res; fq_poly_init(res, *ctx_q);    fq_poly_zero(res, *ctx_q);
    fq_poly_t in1; fq_poly_init(in1, *ctx_q); fq_poly_zero(in1, *ctx_q);
    fq_poly_t in2; fq_poly_init(in2, *ctx_q);  fq_poly_zero(in2, *ctx_q);
    
    /* temporary vectors for conversion between F_2 and int values */
    fq_struct* tmp_vec = _fq_vec_init(this->output_size, *ctx_q);
    _fq_vec_zero(tmp_vec, this->output_size, *ctx_q);

    fq_struct* tmp_vec1 = _fq_vec_init(this->input_size1, *ctx_q);
    _fq_vec_zero(tmp_vec1, this->input_size1, *ctx_q);

    fq_struct* tmp_vec2 = _fq_vec_init(this->input_size2, *ctx_q);
    _fq_vec_zero(tmp_vec1, this->input_size2, *ctx_q);

    
    /* put input into a polynomial */
    _int_vec_2_fq(tmp_vec1, input1, this->input_size1, *ctx_q);
    fq_poly_set_coeffs(in1, tmp_vec1, this->input_size1, *ctx_q);

    _int_vec_2_fq(tmp_vec2, input2, this->input_size2, *ctx_q);
    fq_poly_set_coeffs(in2, tmp_vec2, this->input_size2, *ctx_q);

    printf("just before decoding \n");
    
    /* decoding */
    HQC_decoding(res, in1, in2, SK.y, PK.alpha, this->input_size1, PK.get_n1(),
		 this->output_size, this->r, *ctx, *ctx_q);

    printf("just after decoding \n");
    
    /* put it in vec format */
    fq_poly_get_coeffs(tmp_vec, res, this->output_size, *ctx_q);

    printf("just before conversion \n");
    /* reverse conversion F_2 to int */
    _fq_vec_2_int(output, tmp_vec, this->output_size, *ctx_q);

    /* clear memory */
    _fq_vec_clear(tmp_vec, this->output_size, *ctx_q);
    _fq_vec_clear(tmp_vec1, this->input_size1, *ctx_q);
    _fq_vec_clear(tmp_vec2, this->input_size2, *ctx_q);
    fq_poly_clear(res, *ctx_q);
    fq_poly_clear(in1, *ctx_q);
    fq_poly_clear(in2, *ctx_q);
}
