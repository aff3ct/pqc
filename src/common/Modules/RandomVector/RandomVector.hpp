#ifndef RANDOMVECTOR_H
#define RANDOMVECTOR_H

#include "streampu.hpp"

#include "../../Tools/tools.hpp"


namespace spu
{
    namespace module
    {

	class RandomVector  : public Module {      
  
	private:
    
	    int frame_size;

	public:

	    RandomVector(int frame_size);
	    virtual ~RandomVector();

	protected:

	    virtual void random_vector(int* input, int* output, const int frame_id);

	};
    }
}

#endif // RANDOMVECTOR_H

