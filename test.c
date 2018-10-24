#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "erasure.h"

static void die(const char *msg)
{
	fprintf(stderr, "ERROR: %s\n", msg);
	exit(EXIT_FAILURE);
}

static void die_errno(const char *msg)
{
	perror(msg);
	exit(EXIT_FAILURE);
}

static void* safe_malloc(uintptr_t size)
{
	void *out = malloc(size);
	if(out == NULL) {
		die_errno("malloc");
	}
	return out;
}

static uint64_t get_nanos()
{
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
	return ((uint64_t)tp.tv_sec) * 1000000000 + ((uint64_t)tp.tv_nsec);
}

/* The state must be seeded so that it is not everywhere zero. */
static uint64_t myrand_s[16];
static int myrand_p = 0;

static uint64_t myrand(void)
{
	const uint64_t s0 = myrand_s[myrand_p];
	uint64_t s1 = myrand_s[myrand_p = (myrand_p + 1) & 15];
	s1 ^= s1 << 31;
	myrand_s[myrand_p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30);
	return myrand_s[myrand_p] * UINT64_C(1181783497276652981);
}

static void myseed(uint64_t x)
{
	if(x == 0) {
		x = (uintptr_t)&x;
	}

	for(uint_fast8_t i = 0; i < 16; ++i) {
		x ^= x >> 12;
		x ^= x << 25;
		x ^= x >> 27;
		myrand_s[i] = x * 0x2545F4914F6CDD1D;
	}
}

static void rnd_fill(uint8_t *buffer, uintptr_t size)
{
	for(uintptr_t i = 0; i < size; ++i) {
		buffer[i] = myrand();
	}
}

static void print_buffer(uint8_t *buffer, uintptr_t packets, uintptr_t size)
{
	printf("[");
	if(packets != 0 && size != 0) {
		uintptr_t idx = 0;
		for(uintptr_t p = 0; p < packets; ++p) {
			for(uintptr_t s = 0; s < size; ++s) {
				printf("%02" PRIx8 "", buffer[idx]);
				idx += 1;
			}
			if(p != (packets - 1)) {
				printf(", ");
			}
		}
	}
	printf("]");
}

int main(int argc, char *argv[])
{
	int k = 7;
	int m = 3;
	uintptr_t msg_size = 3;
	uint64_t seed = (uint64_t)time(NULL);

	if(argc >= 4) {
		k = atoi(argv[1]);
		m = atoi(argv[2]);
		msg_size = atoi(argv[3]);
	}

	uintptr_t overall_data_size = msg_size * k;

	seed += k;
	seed += m;
	seed += msg_size;

	if(argc == 5) {
		seed = atoi(argv[4]);
	}

	printf("packets = %d\n", k);
	printf("redundant packets = %d\n", m);
	printf("packet size = %" PRIuPTR "\n", msg_size);
	printf("seed = %" PRIu64 "\n", seed);
	myseed(seed);

	struct erasure_coder *coder = (struct erasure_coder*)safe_malloc(sizeof(struct erasure_coder));

	uintptr_t size = ec_state_size(k, m);
	if(size == 0) {
		die("ec_state_size");
	}
	printf("error-correction state size is %" PRIuPTR " bytes\n", size);

	uint8_t *state = (uint8_t*)safe_malloc(size);
	int ret = ec_init(coder, k, m, state);
	if(!ret) {
		die("ec_init");
	}

	uint8_t *msg = (uint8_t*)safe_malloc(msg_size * k);
	rnd_fill(msg, msg_size * k);

	uint8_t *add = (uint8_t*)safe_malloc(msg_size * m);
	rnd_fill(add, msg_size * m);

	uint8_t *msg_backup = (uint8_t*)safe_malloc(msg_size * k);
	memcpy(msg_backup, msg, msg_size * k);


	uint8_t **data_ptrs = (uint8_t**)safe_malloc(sizeof(uint8_t*) * k);
	uint8_t **coding_ptrs = (uint8_t**)safe_malloc(sizeof(uint8_t*) * m);
	for(int i = 0; i < k; ++i) {
		data_ptrs[i] = msg + msg_size * i;
	}
	for(int i = 0; i < m; ++i) {
		coding_ptrs[i] = add + msg_size * i;
	}

	uint64_t enc_time = get_nanos();
	ec_encode(coder, data_ptrs, coding_ptrs, msg_size);
	enc_time = get_nanos() - enc_time;
	printf("encoded in %" PRIu64 "ns -> %" PRIu64 "MiB/s\n", enc_time, ((overall_data_size * 1000000000) / enc_time) / (1024 * 1024));

	uint8_t *add_backup = (uint8_t*)safe_malloc(msg_size * m);
	memcpy(add_backup, add, msg_size * m);

	printf("\noriginal data\n");
	printf("data: ");
	print_buffer(msg, k, msg_size);
	printf("\n");

	printf("add:  ");
	print_buffer(add, m, msg_size);
	printf("\n\n");

	if(memcmp(msg, msg_backup, msg_size * k) != 0) {
		die("data changed after encoding");
	}

	// delete up to 'm' packets
	uint8_t *erased = (uint8_t*)safe_malloc(k + m);
	memset(erased, 0, k + m);

	for(int i = 0; i < m; ++i) {
		int idx = myrand() % (k+m);
		erased[idx] = 1;
		if(idx < k) {
			printf("deleting data packet %d\n", idx);
			memset(msg + idx * msg_size, 0, msg_size);
		}
		else {
			printf("deleting add packet %d\n", idx - k);
			memset(add + (idx - k) * msg_size, 0, msg_size);
		}
	}

	printf("\ndata after deletion\n");
	printf("data: ");
	print_buffer(msg, k, msg_size);
	printf("\n");

	printf("add:  ");
	print_buffer(add, m, msg_size);
	printf("\n");

	uint64_t dec_time = get_nanos();
	ec_decode(coder, erased, data_ptrs, coding_ptrs, msg_size);
	dec_time = get_nanos() - dec_time;
	printf("decoded in %" PRIu64 "ns -> %" PRIu64 "MiB/s\n", dec_time, ((overall_data_size * 1000000000) / dec_time) / (1024 * 1024));

	printf("\nrestored data\n");
	printf("data: ");
	print_buffer(msg, k, msg_size);
	printf("\n");

	printf("add:  ");
	print_buffer(add, m, msg_size);
	printf("\n\n");

	if(memcmp(msg, msg_backup, msg_size * k) != 0) {
		die("data not restored");
	}
	if(memcmp(add, add_backup, msg_size * m) != 0) {
		die("add not restored");
	}
	printf("SUCCESS\n");

	free(erased);
	free(add_backup);
	free(coding_ptrs);
	free(data_ptrs);
	free(msg_backup);
	free(add);
	free(msg);
	free(state);
	free(coder);

	return EXIT_SUCCESS;
}
