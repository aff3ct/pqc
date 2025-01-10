#ifndef SYNDCOMPARATOR_H
#define SYNDCOMPARATOR_H

#include <streampu.hpp>

namespace spu
{
    namespace module
    {

	class SyndComparator  : public Stateful {

	private:
    
	    int frame_size;
	    int synd_size;
	    
	public:

	    SyndComparator(int frame_size, int synd_size);
	    virtual ~SyndComparator();

	protected:

	    virtual void compare(int* input1, int* input2, int* input3, int *output,
				 const int frame_id);

	};
    }
}

#endif // SYNDCOMPARATOR_H
