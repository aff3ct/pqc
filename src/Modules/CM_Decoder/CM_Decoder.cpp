#include <iostream>
#include "Modules/CM_Decoder/CM_Decoder.hpp"

using namespace spu;
using namespace spu::module;
using namespace std;

CM_Decoder:: CM_Decoder(int frame_size, int synd_size, int weight, CM_secret_key& SK) :
    Stateful(),
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

    vec_GF2E s(INIT_SIZE,this->frame_size); s.SetLength(this->frame_size);  /* Actual size : synd_size ... size is extended to match tmp_e for addition */ 
    vec_GF2E e(INIT_SIZE,this->frame_size); e.SetLength(this->frame_size);   

    int input_size_byte = (this->synd_size + 7) /8;

    GF2X tmp = GF2XFromBytes((uint8_t*)input,input_size_byte);

    // GF2X to vec_GF2E (maybe not optimal)
    for(int i = 0; i< this->synd_size;i++){
        s[i] = to_GF2E(tmp[i]);
    }   

    /* decoding */
    int b = CM_syndrome_decoding_bin(e,s,SK.get_alpha(),SK.g,this->frame_size,this->weight);

    e += s;
    
    /* conversion GF2E to int */   
    vec_gf2e_to_int(e,output);


}