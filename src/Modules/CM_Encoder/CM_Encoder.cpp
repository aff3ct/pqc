#include <iostream>
#include "Modules/CM_Encoder/CM_Encoder.hpp"


using namespace spu;
using namespace spu::module;
using namespace std;

CM_Encoder:: CM_Encoder(int frame_size, int out_size, CM_public_key& PK) :
    Stateful(),
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
    vec_GF2 s(INIT_SIZE,this->out_size); s.SetLength(this->out_size); 
    vec_GF2 e(INIT_SIZE,this->frame_size); e.SetLength(this->frame_size); 
    
    /* Sizes for bit-stacking in integers and bytes */ 
    /* Input */
    int frame_size_int = (frame_size + 63) / 64;
    int frame_size_byte = (frame_size + 7) / 8;
    /* Output */
    int out_size_bytes = (this->out_size + 7) /8;

    uint64_t input_int[frame_size_int];
    bit_to_uint(input_int,input,frame_size);
    e = VectorCopy(GF2XFromBytes((uint8_t*)input_int,frame_size_byte),this->frame_size);

    /* Encoding */
    CM_encoding(s,PK.T,e);

    uint8_t output_bytes[out_size_bytes];
    BytesFromGF2X(output_bytes, conv<GF2X>(s), out_size_bytes);
    memcpy(output,output_bytes,out_size_bytes);
}