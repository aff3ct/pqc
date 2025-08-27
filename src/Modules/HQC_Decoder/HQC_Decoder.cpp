#include <iostream>
#include "Modules/HQC_Decoder/HQC_Decoder.hpp"

using namespace spu;
using namespace spu::module;
using namespace std;

HQC_Decoder::HQC_Decoder(HQC_secret_key &SK, HQC_public_key &PK, int k, int len,
                         int r) : Stateful(),
                                  input_size1(PK.get_n()),
                                  input_size2(len),
                                  output_size(k),
                                  r(r)
{
    int input_size1_64t = (((input_size1) / (8)) + ((input_size1) % (8) == 0 ? 0 : 1));
    int input_size2_64t = (((input_size2) / (8)) + ((input_size2) % (8) == 0 ? 0 : 1));

    this->set_name("HQC_Decoder");
    this->set_short_name("HQC_Decoder");

    auto &t = create_task("hqc_decoder");
    auto input1 = create_socket_in<uint64_t>(t, "input1", input_size1_64t);
    auto input2 = create_socket_in<uint64_t>(t, "input2", input_size2_64t);
    auto output = create_socket_out<int>(t, "output", output_size);

    this->create_codelet(t, [input1, input2, output, &PK, &SK](Module &m, runtime::Task &t, const size_t frame_id) -> int
                         {
	static_cast<HQC_Decoder&>(m).hqc_decoder(static_cast<uint64_t*>(t[input1].get_dataptr()),
						 static_cast<uint64_t*>(t[input2].get_dataptr()),
						 static_cast<int*>(t[output].get_dataptr()),
						 static_cast<HQC_secret_key&>(SK),
						 static_cast<HQC_public_key&>(PK),
						 frame_id);
	return 0; });
}

HQC_Decoder::~HQC_Decoder()
{
}

void HQC_Decoder::hqc_decoder(uint64_t *input1, uint64_t *input2, int *output, const HQC_secret_key &SK,
                              const HQC_public_key &PK, const int frame_id)
{
    /* Encoding result */
    GF2EX res(INIT_SIZE,this->output_size);
    res.SetLength(this->output_size);

    /* u, v as in HQC */
    GF2X u = GF2XFromBytes((uint8_t *)input1, this->input_size1);
    GF2X v = GF2XFromBytes((uint8_t *)input2, this->input_size2);

    /* decoding */
    HQC_decoding(res, u, v, SK.y, PK.alpha, this->input_size1, PK.get_n1(), this->output_size, this->r);
   
    /* Conversion to int vector */
    GF2EX_to_int(res, output, this->output_size);
}
