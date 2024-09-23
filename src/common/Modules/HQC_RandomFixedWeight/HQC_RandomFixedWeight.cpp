#include <random>
#include "HQC_RandomFixedWeight.hpp"
using namespace spu;
using namespace spu::module;

HQC_RandomFixedWeight::HQC_RandomFixedWeight(int frame_size, int weight) :
    Module(),
    frame_size(frame_size),
    weight(weight) {

    this->set_name("HQC_RandomFixedWeight");
    this->set_short_name("HQC_RandomFixedWeight");

    auto &t = create_task("random_fixed_weight");
    auto input   = create_socket_in<int>(t, "input", frame_size);
    auto output   = create_socket_out<int>(t, "output", frame_size);


    this->create_codelet(t, [input, output](Module &m, runtime::Task &t,
					    const size_t frame_id) -> int {
	static_cast<HQC_RandomFixedWeight&>(m).random_fixed_weight(static_cast<int*>(t[input].get_dataptr()),
								    static_cast<int*>(t[output].get_dataptr()),
								    frame_id);
	return 0;
    });
}

HQC_RandomFixedWeight::~HQC_RandomFixedWeight() {
}

void HQC_RandomFixedWeight::random_fixed_weight(int* input, int* output, const int frame_id) { 

    int e[frame_size];
    for (auto i = 0; i < frame_size; i++) {
	e[i] = 0;
	output[i] = input[i];
    }
  
    hqc_gen_e(e, frame_size, weight);
  
    for(auto i = 0; i < frame_size; i++) {
	output[i] += e[i];
    }
}
