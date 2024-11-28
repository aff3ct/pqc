#include <random>
#include "RandomVector.hpp"

using namespace spu;
using namespace spu::module;

RandomVector::RandomVector(int frame_size) :
    Module(),
    frame_size(frame_size) {

    this->set_name("RandomVector");
    this->set_short_name("RandomVector");

    auto &t = create_task("random_vector");
    auto input   = create_socket_in<int>(t, "input", frame_size);
    auto output   = create_socket_out<int>(t, "output", frame_size);


    this->create_codelet(t, [input, output](Module &m, runtime::Task &t,
					    const size_t frame_id) -> int {
	static_cast<RandomVector&>(m).random_vector(static_cast<int*>(t[input].get_dataptr()),
							       static_cast<int*>(t[output].get_dataptr()),
							       frame_id);
	return 0;
    });
}

RandomVector::~RandomVector() {
}    

void RandomVector::random_vector(int* input, int* output, const int frame_id) { 
    random_bits(output, this->frame_size);
     
}
