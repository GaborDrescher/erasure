#include "erasure.h"
#include <string.h> // memcpy

#define GF_INVERT_LOOKUP

static uint8_t gf_mul(uint8_t a, uint8_t b)
{
	if (a == 0 || b == 0) return 0;

	uint8_t p = 0;
	while (b) {
		if((b & 1) == 1) {
			p ^= a;
		}
		if(a & 0x80) {
			a = (a << 1) ^ 0x11D;
		}
		else {
			a <<= 1;
		}
		b >>= 1;
	}
	return p;
}

#ifdef GF_INVERT_LOOKUP
	static const uint8_t gf_inv_lookup[256] = {
		0, 1, 142, 244, 71, 167, 122, 186, 173, 157, 221, 152, 61, 170, 93,
		150, 216, 114, 192, 88, 224, 62, 76, 102, 144, 222, 85, 128, 160, 131,
		75, 42, 108, 237, 57, 81, 96, 86, 44, 138, 112, 208, 31, 74, 38, 139,
		51, 110, 72, 137, 111, 46, 164, 195, 64, 94, 80, 34, 207, 169, 171, 12,
		21, 225, 54, 95, 248, 213, 146, 78, 166, 4, 48, 136, 43, 30, 22, 103,
		69, 147, 56, 35, 104, 140, 129, 26, 37, 97, 19, 193, 203, 99, 151, 14,
		55, 65, 36, 87, 202, 91, 185, 196, 23, 77, 82, 141, 239, 179, 32, 236,
		47, 50, 40, 209, 17, 217, 233, 251, 218, 121, 219, 119, 6, 187, 132,
		205, 254, 252, 27, 84, 161, 29, 124, 204, 228, 176, 73, 49, 39, 45, 83,
		105, 2, 245, 24, 223, 68, 79, 155, 188, 15, 92, 11, 220, 189, 148, 172,
		9, 199, 162, 28, 130, 159, 198, 52, 194, 70, 5, 206, 59, 13, 60, 156,
		8, 190, 183, 135, 229, 238, 107, 235, 242, 191, 175, 197, 100, 7, 123,
		149, 154, 174, 182, 18, 89, 165, 53, 101, 184, 163, 158, 210, 247, 98,
		90, 133, 125, 168, 58, 41, 113, 200, 246, 249, 67, 215, 214, 16, 115,
		118, 120, 153, 10, 25, 145, 20, 63, 230, 240, 134, 177, 226, 241, 250,
		116, 243, 180, 109, 33, 178, 106, 227, 231, 181, 234, 3, 143, 211, 201,
		66, 212, 232, 117, 127, 255, 126, 253
	};
	
	static uint8_t galois_single_invert(uint8_t y)
	{
		return gf_inv_lookup[y];
	}
#else
	static uint8_t galois_single_invert(uint8_t y)
	{
		// TODO use extended euclidean algorithm
		for(uint16_t i = 0; i < 256; ++i) {
			if(gf_mul(y, i) == 1) {
				return i;
			}
		}
	
		return 0;
	}
#endif

static void galois_region_multiply(uint8_t *region, uint8_t multby,
	uintptr_t nbytes, uint8_t *r2, uint8_t add)
{
	uint8_t *dest = region;
	if(r2 != 0) {
		dest = r2;
	}
	for(uintptr_t i = 0; i < nbytes; ++i) {
		uint8_t val = region[i];
		uint8_t res = gf_mul(val, multby);
		if(add) {
			dest[i] ^= res;
		}
		else {
			dest[i] = res;
		}
	}
}

static void galois_region_xor(uint8_t *src, uint8_t *dest, uintptr_t nbytes)
{
	for(uintptr_t i = 0; i < nbytes; ++i) {
		dest[i] ^= src[i];
	} 
}

static int ec_check(uint8_t k_in, uint8_t m_in)
{
	uintptr_t k = k_in;
	uintptr_t m = m_in;

	if((k + m) > 256) {
		return 0;
	}
	return 1;
}

uint32_t ec_state_size(uint8_t k_in, uint8_t m_in)
{
	if(!ec_check(k_in, m_in)) {
		return 0;
	}

	uint_fast32_t k = k_in;
	uint_fast32_t m = m_in;
	uint_fast32_t size = k*k + k + k*k + k + k*m;
	return size;
}

int ec_init(struct erasure_coder *coder, uint8_t k_in, uint8_t m_in, uint8_t *state)
{
	if(coder == NULL || state == NULL) {
		return 0;
	}

	coder->k = k_in;
	coder->m = m_in;

	uintptr_t k = k_in;
	uintptr_t m = m_in;

	coder->tmpmat = state;

	state += k * k;
	coder->dm_ids = state;

	state += k;
	coder->decoding_matrix = state;

	state += k * k;
	coder->tmpids = state;

	state += k;
	coder->cauchy_matrix = state;

	uint8_t *cauchy_matrix = coder->cauchy_matrix;
	for(uintptr_t i = 0; i < m; i++) {
		for(uintptr_t j = 0; j < k; j++) {
			*cauchy_matrix = galois_single_invert(i ^ (m+j));
			cauchy_matrix += 1;
		}
	}

	return 1;
}

static void matrix_dot(struct erasure_coder *coder,
	uint8_t *matrix_row, uint8_t *src_ids, uintptr_t dest_id, uint8_t **data_ptrs,
	uint8_t **coding_ptrs, uintptr_t size)
{
	uint8_t *dptr, *sptr;
	uintptr_t k = coder->k;
	uint8_t init = 0;

	dptr = (dest_id < k) ? data_ptrs[dest_id] : coding_ptrs[dest_id-k];

	for(uintptr_t i = 0; i < k; i++) {
		if(matrix_row[i] == 1) {
			if(src_ids == 0) {
				sptr = data_ptrs[i];
			}
			else if(src_ids[i] < k) {
				sptr = data_ptrs[src_ids[i]];
			}
			else {
				sptr = coding_ptrs[src_ids[i]-k];
			}

			if(init == 0) {
				memcpy(dptr, sptr, size);
				init = 1;
			}
			else {
				galois_region_xor(sptr, dptr, size);
			}
		}
	}

	for(uintptr_t i = 0; i < k; i++) {
		if(matrix_row[i] != 0 && matrix_row[i] != 1) {
			if(src_ids == 0) {
				sptr = data_ptrs[i];
			}
			else if(src_ids[i] < k) {
				sptr = data_ptrs[src_ids[i]];
			}
			else {
				sptr = coding_ptrs[src_ids[i]-k];
			}
			galois_region_multiply(sptr, matrix_row[i], size, dptr, init);
			init = 1;
		}
	}
}

static void jerasure_invert_matrix(struct erasure_coder *coder, uint8_t *mat,
	uint8_t *inv)
{
	const uintptr_t k = coder->k;
	const uintptr_t rows = k;
	const uintptr_t cols = k;

	for(uintptr_t i = 0, ki = 0; i < rows; i++) {
		for(uintptr_t j = 0; j < cols; j++) {
			inv[ki] = (i == j) ? 1 : 0;
			ki += 1;
		}
	}

	for(uintptr_t i = 0; i < cols; i++) {
		uintptr_t row_start = cols * i;
		if(mat[row_start + i] == 0) { 
			uintptr_t jtmp = i + 1;
			while(jtmp < rows && mat[cols * jtmp + i] == 0) {
				jtmp += 1;
			}
			uintptr_t rs2 = jtmp * cols;
			for(uintptr_t ki = 0; ki < cols; ki++) {
				uintptr_t tmp = mat[row_start + ki];
				mat[row_start + ki] = mat[rs2 + ki];
				mat[rs2 + ki] = tmp;

				tmp = inv[row_start + ki];
				inv[row_start + ki] = inv[rs2 + ki];
				inv[rs2+ki] = tmp;
			}
		}
		uintptr_t tmp = mat[row_start + i];
		if(tmp != 1) {
			uintptr_t inverse = galois_single_invert(tmp);
			for(uintptr_t j = 0; j < cols; j++) { 
				mat[row_start + j] = gf_mul(mat[row_start + j], inverse);
				inv[row_start + j] = gf_mul(inv[row_start + j], inverse);
			}
		}
		uintptr_t ki = row_start + i;
		for(uintptr_t j = i + 1; j != cols; j++) {
			ki += cols;
			if(mat[ki] != 0) {
				if(mat[ki] == 1) {
					uintptr_t rs2 = cols * j;
					for(uintptr_t x = 0; x < cols; x++) {
						mat[rs2 + x] ^= mat[row_start + x];
						inv[rs2 + x] ^= inv[row_start + x];
					}
				} else {
					tmp = mat[ki];
					uintptr_t rs2 = cols * j;
					for(uintptr_t x = 0; x < cols; x++) {
						mat[rs2 + x] ^= gf_mul(tmp, mat[row_start + x]);
						inv[rs2 + x] ^= gf_mul(tmp, inv[row_start + x]);
					}
				}
			}
		}
	}

	for(uintptr_t i = rows; i > 0; i--) {
		uintptr_t idx = i - 1;
		uintptr_t row_start = idx * cols;
		for(uintptr_t j = 0; j < idx; j++) {
			uintptr_t rs2 = j * cols;
			if (mat[rs2 + idx] != 0) {
				uintptr_t tmp = mat[rs2 + idx];
				mat[rs2 + idx] = 0; 
				for(uintptr_t ki = 0; ki < cols; ki++) {
					inv[rs2 + ki] ^= gf_mul(tmp, inv[row_start + ki]);
				}
			}
		}
	}
}

static void make_decoding_matrix(struct erasure_coder *coder,
	uint8_t *matrix, uint8_t *erased, uint8_t *decoding_matrix,
	uint8_t *dm_ids)
{
	uintptr_t k = coder->k;
	uint8_t *tmpmat = coder->tmpmat;

	for(uintptr_t i = 0, j = 0; j < k; i++) {
		if(erased[i] == 0) {
			dm_ids[j] = i;
			j++;
		}
	}

	for(uintptr_t i = 0; i < k; i++) {
		if(dm_ids[i] < k) {
			for(uintptr_t j = 0; j < k; j++) {
				tmpmat[i*k+j] = 0;
			}
			tmpmat[i*k+dm_ids[i]] = 1;
		}
		else {
			for(uintptr_t j = 0; j < k; j++) {
				tmpmat[i*k+j] = matrix[(dm_ids[i]-k)*k+j];
			}
		}
	}

	jerasure_invert_matrix(coder, tmpmat, decoding_matrix);
}

void ec_encode(struct erasure_coder *coder, uint8_t **data_ptrs,
	uint8_t **coding_ptrs, uintptr_t size)
{
	uintptr_t k = coder->k;
	uintptr_t m = coder->m;
	uint8_t *matrix = coder->cauchy_matrix;
	for(uintptr_t i = 0; i < m; i++) {
		matrix_dot(coder, matrix + (i * k), 0, k+i, data_ptrs, coding_ptrs, size);
	}
}

void ec_decode(struct erasure_coder *coder, uint8_t *erased,
	uint8_t **data_ptrs, uint8_t **coding_ptrs, uintptr_t size)
{
	uintptr_t k = coder->k;
	uintptr_t m = coder->m;

	uint8_t *decoding_matrix = coder->decoding_matrix;
	uint8_t *dm_ids = coder->dm_ids;
	uint8_t *tmpids = coder->tmpids;
	uint8_t *matrix = coder->cauchy_matrix;

	uintptr_t edd = 0;
	for(uintptr_t i = 0; i < k; i++) {
		if(erased[i]) {
			edd += 1;
		}
	}

	if(edd > 0) {
		make_decoding_matrix(coder, matrix, erased, decoding_matrix, dm_ids);
	}

	for(uintptr_t i = 0; edd > 0 && i < k; i++) {
		if(erased[i]) {
			matrix_dot(coder, decoding_matrix+(i*k), dm_ids, i, data_ptrs, coding_ptrs, size);
			edd -= 1;
		}
	}

	if(edd > 0) {
		for(uintptr_t i = 0; i < k; i++) {
			tmpids[i] = i;
		}
		matrix_dot(coder, matrix, tmpids, k, data_ptrs, coding_ptrs, size);
	}

	for(uintptr_t i = 0; i < m; i++) {
		if(erased[k+i]) {
			matrix_dot(coder, matrix+(i*k), 0, i+k, data_ptrs, coding_ptrs, size);
		}
	}
}
