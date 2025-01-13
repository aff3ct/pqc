#include <random>
#include "Modules/RandomFixedWeight/RandomFixedWeight.hpp"

using namespace spu;
using namespace spu::module;

RandomFixedWeight::RandomFixedWeight(int frame_size, int weight, int N, int tau) :
    Stateful(),
    frame_size(frame_size),
    weight(weight),
    N(N),
    tau(tau) {

    this->set_name("RandomFixedWeight");
    this->set_short_name("RandomFixedWeight");

    auto &t = create_task("random_fixed_weight");
    auto input   = create_socket_in<int>(t, "input", frame_size);
    auto output   = create_socket_out<int>(t, "output", frame_size);


    this->create_codelet(t, [input, output](Module &m, runtime::Task &t,
					    const size_t frame_id) -> int {
	static_cast<RandomFixedWeight&>(m).random_fixed_weight(static_cast<int*>(t[input].get_dataptr()),
							       static_cast<int*>(t[output].get_dataptr()),
							       frame_id);
	return 0;
    });
}

RandomFixedWeight::~RandomFixedWeight() {
}    

void RandomFixedWeight::random_fixed_weight(int* input, int* output, const int frame_id) { 

    // std::random_device rd;
    // std::mt19937 rand_gen(rd());
    // std::uniform_int_distribution<> dis(0, frame_size-1);

    int e[frame_size];
    for (auto i = 0; i < frame_size; i++) e[i] = 0;
  
    cm_gen_e(e, frame_size, weight, N, tau);
  
    for(auto i = 0; i < frame_size; i++)
	{
	    output[i] += e[i];
	}
  
}







// #include <random>
// #include "MySource.hpp"

// using namespace spu;
// using namespace spu::module;

// MySource::MySource(int frame_size) : Stateful(), frame_size(frame_size) {

//   this->set_name("MySource");
//   this->set_short_name("MySource");

//   auto &t = create_task("generate");
//   auto output   = create_socket_out<int>(t, "output", frame_size);

//   this->create_codelet(t, [output](Module &m, runtime::Task &t, const size_t frame_id) -> int {
//     static_cast<MySource&>(m).generate(  static_cast<int*>(t[output].get_dataptr()),
// 					 frame_id);
//     return 0;
//   });

// }

// MySource::~MySource() {
// }    

// void MySource::generate(int *output, const int frame_id) {
        
//   std::random_device rd;
//   std::mt19937 rand_gen(rd());
//   std::uniform_int_distribution<> dis(0, 63);

//   int i = 0;
//   while(i < frame_size)
//     {
//       output[i++] = dis(rand_gen);
//     }        
// }
