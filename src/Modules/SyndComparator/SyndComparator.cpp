#include <random>
#include "Modules/SyndComparator/SyndComparator.hpp"

using namespace spu;
using namespace spu::module;

SyndComparator::SyndComparator(int frame_size, int synd_size) :
    Stateful(),
    frame_size(frame_size),
    synd_size(synd_size) {

    this->set_name("SyndComparator");
    this->set_short_name("SyndComparator");

    auto &t = create_task("compare");
    auto input1   = create_socket_in<int>(t, "input1", frame_size);
    auto input2   = create_socket_in<int>(t, "input2", frame_size);
    auto input3   = create_socket_in<int>(t, "input3", synd_size);
    auto output   = create_socket_out<int>(t, "output", frame_size);

    this->create_codelet(t, [input1, input2, input3, output](Module &m, runtime::Task &t, const size_t frame_id) -> int {
	static_cast<SyndComparator&>(m).compare(   static_cast<int*>(t[input1].get_dataptr()),
					       static_cast<int*>(t[input2].get_dataptr()),
					       static_cast<int*>(t[input3].get_dataptr()),
					       static_cast<int*>(t[output].get_dataptr()),
					       frame_id);
	return 0;
    });

}

SyndComparator::~SyndComparator() {
}    

void SyndComparator::compare(int *input1, int *input2, int *input3, int *output,
			 const int frame_id) {
        
    for(auto i = 0; i < frame_size; i++)
	{
	    if (i < synd_size) {
		if(input1[i] == (input2[i] + input3[i]) % 2  ) output[i] = 1;
		else output[i] = 0;            
	    } else {	    
		if(input1[i] == input2[i]) output[i] = 1;
		else output[i] = 0;            
	    }
	}        
}
