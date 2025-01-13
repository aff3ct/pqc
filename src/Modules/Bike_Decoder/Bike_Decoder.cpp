#include <iostream>
#include "Modules/Bike_Decoder/Bike_Decoder.hpp"

using namespace spu;
using namespace spu::module;
using namespace std;

Bike_Decoder:: Bike_Decoder(int weight, int NbIter, int tau, Bike_secret_key& SK) :
    Stateful(),
    frame_size(2*SK.get_r()),
    input_size(SK.get_r()),
    weight(weight),
    NbIter(NbIter),
    tau(tau) {

    this->set_name("Bike_Decoder");
    this->set_short_name("Bike_Decoder");

    auto &t = create_task("bike_decoder");
    auto input   = create_socket_in<int>(t, "input", input_size);
    auto output  = create_socket_out<int>(t, "output", frame_size);

    this->create_codelet(t, [input, output, &SK](Module &m, runtime::Task &t,
						 const size_t frame_id) -> int {
	static_cast<Bike_Decoder&>(m).bike_decoder(static_cast<int*>(t[input].get_dataptr()),
						   static_cast<int*>(t[output].get_dataptr()),
						   static_cast<Bike_secret_key&>(SK),
						   frame_id);
	return 0;
    }			 	);
}


Bike_Decoder:: ~Bike_Decoder() {
}


void
Bike_Decoder:: bike_decoder(int* input, int* output, const Bike_secret_key& SK,
			    const int frame_id) {

    fq_ctx_t* ctx_q = SK.get_ctx_q(); /* finite field F_2 */
    int r = this->input_size;
    fq_poly_t P; fq_poly_init(P, *ctx_q); fq_poly_set_cyclic(P, r, *ctx_q);
    
    /* transform syndrome from vector to polynomial  */
    fq_struct* s =  _fq_vec_init(r, *ctx_q);
    _int_vec_2_fq(s, input, r, *ctx_q);

    /* compute s = s * h0 mod P */
    fq_poly_t sp; fq_poly_init(sp, *ctx_q); fq_poly_zero(sp, *ctx_q);
    fq_poly_set_coeffs(sp, s, r, *ctx_q);
    fq_poly_mulmod(sp, sp, SK.h0, P, *ctx_q);

    /* put sp in a vector (maybe should not ?) */
    fq_struct* ss = _fq_vec_init(r, *ctx_q);
    for (int k = 0; k < r; k++) {
	fq_poly_get_coeff(&ss[k], sp, k, *ctx_q);
    }

    /* decoding */
    fq_struct* res = _fq_vec_init(this->frame_size, *ctx_q);
    int b = Bike_decoding_v2(res, ss, SK.h0, SK.h1, r, (this->weight)/2,
			     this->NbIter, this->tau, *ctx_q);

    
    /* reverse conversion F_2 to int */
    _fq_vec_2_int(output, res, this->frame_size, *ctx_q);
    
    
    /* clear memory */
    fq_poly_clear(P, *ctx_q);
    fq_poly_clear(sp, *ctx_q);
    _fq_vec_clear(s, r, *ctx_q);
    _fq_vec_clear(ss, r, *ctx_q);
    _fq_vec_clear(res, r, *ctx_q);
}
