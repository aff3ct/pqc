#include <iostream>
#include "Modules/Bike_Encoder/Bike_Encoder.hpp"

using namespace spu;
using namespace spu::module;
using namespace std;

Bike_Encoder:: Bike_Encoder(Bike_public_key& PK) :
    Stateful(),
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
    GF2X s(INIT_SIZE,this->output_size);
    s.SetLength(this->output_size);

    Bike_encoding(s,input,PK.h,this->output_size); /* Encoding */

    GF2X_to_int(s,output,this->output_size); /* Converts GF2X polynomial output to int vec */

}
