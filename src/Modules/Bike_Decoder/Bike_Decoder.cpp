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

    int r = this->input_size;

    GF2X P = gf2x_set_cyclic(r); /* X^r + 1 */ 
    vec_GF2 res; res.SetMaxLength(this->frame_size); res.SetLength(this->frame_size);

    GF2X s = int_vec_to_GF2X(input,r); /* Converts input into a GF2 polynomial */
    s = MulMod(s,SK.h0,P);        
    
    int b = Bike_decoding_v2(res, s, SK.h0, SK.h1,r,(this->weight)/2, this->NbIter, this->tau);
    vec_gf2_to_int(res,output);
    
}
