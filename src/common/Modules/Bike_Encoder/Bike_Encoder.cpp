#include "Bike_Encoder.hpp"
#include <iostream>

using namespace spu;
using namespace spu::module;
using namespace std;

Bike_Encoder:: Bike_Encoder(Bike_public_key& PK) :
    Module(),
    frame_size(2*PK.get_r()),
    output_size(PK.get_r()) {

    this->set_name("Bike_Encoder");
    this->set_short_name("Bike_Encoder");

    auto &t = create_task("bike_encoder");
    auto input   = create_socket_in<int>(t, "input", frame_size);
    auto output  = create_socket_out<int>(t, "output", output_size);

    this->create_codelet(t, [input, output, &PK](Module &m, runtime::Task &t,
						 const size_t frame_id) -> int {
	static_cast<Bike_Encoder&>(m).bike_encoder(static_cast<int*>(t[input].get_dataptr()),
						   static_cast<int*>(t[output].get_dataptr()),
						   static_cast<Bike_public_key&>(PK),
						   frame_id);
	return 0;
    }			 	);
}


Bike_Encoder:: ~Bike_Encoder() {
}


void
Bike_Encoder:: bike_encoder(int* input, int* output, const Bike_public_key& PK,
			    const int frame_id) {
    fq_ctx_t* ctx_q = PK.get_ctx_q(); /* finite field F_2 */

    fq_t tmp; fq_init(tmp, *ctx_q);
    
    /* temporary vectors for conversion between F_2 and int values */
    fq_struct* tmp_s = _fq_vec_init(this->output_size, *ctx_q);

    fq_poly_t s;
    fq_poly_init(s, *ctx_q); 
    
    /* encoding */
    Bike_encoding(s, input, PK.h, this->output_size, *ctx_q);

    /* put it in vec format */
    for (int i = 0; i < this->output_size; ++i) {
	fq_poly_get_coeff(&tmp_s[i], s, i, *ctx_q);
	/* fq_set(tmp_s[i], tmp, *ctx_q); */
    }
    
    /* reverse conversion F_2 to int */
    _fq_vec_2_int(output, tmp_s, this->output_size, *ctx_q);
    
    /* clear memory */
    /* _fq_vec_clear(tmp_e, this->frame_size, *ctx_q); */
    _fq_vec_clear(tmp_s, this->output_size, *ctx_q);
    fq_poly_clear(s, *ctx_q);
}
