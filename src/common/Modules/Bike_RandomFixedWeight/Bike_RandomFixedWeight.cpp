#include <random>
#include "Bike_RandomFixedWeight.hpp"

using namespace spu;
using namespace spu::module;

Bike_RandomFixedWeight::Bike_RandomFixedWeight(int frame_size, int weight) :
    Module(),
    frame_size(frame_size),
    weight(weight) {

    this->set_name("Bike_RandomFixedWeight");
    this->set_short_name("Bike_RandomFixedWeight");

    auto &t = create_task("random_fixed_weight");
    auto input   = create_socket_in<int>(t, "input", frame_size);
    auto output   = create_socket_out<int>(t, "output", frame_size);


    this->create_codelet(t, [input, output](Module &m, runtime::Task &t,
					    const size_t frame_id) -> int {
	static_cast<Bike_RandomFixedWeight&>(m).random_fixed_weight(static_cast<int*>(t[input].get_dataptr()),
							       static_cast<int*>(t[output].get_dataptr()),
							       frame_id);
	return 0;
    });
}

Bike_RandomFixedWeight::~Bike_RandomFixedWeight() {
}    

void Bike_RandomFixedWeight::random_fixed_weight(int* input, int* output, const int frame_id) { 

    // std::random_device rd;
    // std::mt19937 rand_gen(rd());
    // std::uniform_int_distribution<> dis(0, frame_size-1);

    int e[frame_size];
    for (auto i = 0; i < frame_size; i++) {
	e[i] = 0;
	output[i] = input[i];
    }
  
    bike_gen_e(e, frame_size, weight);
  
    for(auto i = 0; i < frame_size; i++)
	{
	    output[i] += e[i];
	}
  
    // for(auto i = 0; i < weight; i++)
    //   {
    //     output[dis(rand_gen)] += 1;
    //   }
}







// #include <random>
// #include "MySource.hpp"

// using namespace spu;
// using namespace spu::module;

// MySource::MySource(int frame_size) : Module(), frame_size(frame_size) {

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
