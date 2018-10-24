#ifndef ERASURE_CODER_HEADER
#define ERASURE_CODER_HEADER

#include <inttypes.h>

struct erasure_coder
{
	uint8_t *tmpmat;
	uint8_t *dm_ids;
	uint8_t *decoding_matrix;
	uint8_t *tmpids;
	uint8_t *cauchy_matrix;

	uint8_t k;
	uint8_t m;
};

uint32_t ec_state_size(uint8_t k_in, uint8_t m_in);
int ec_init(struct erasure_coder *coder, uint8_t k_in, uint8_t m_in, uint8_t *state);
void ec_encode(struct erasure_coder *coder, uint8_t **data_ptrs, uint8_t **coding_ptrs, uintptr_t size);
void ec_decode(struct erasure_coder *coder, uint8_t *erased, uint8_t **data_ptrs, uint8_t **coding_ptrs, uintptr_t size);

#endif /* ERASURE_CODER_HEADER */
