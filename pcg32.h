#ifndef PCG32_H
#define PCG32_H
#include <stdint.h>
static uint64_t       state      = 0x4d595df4d0f33173;		// Or something seed-dependent
static uint64_t const multiplier = 6364136223846793005u;
static uint64_t const increment  = 1442695040888963407u;	// Or an arbitrary odd constant

static uint32_t rotr32(uint32_t x, unsigned r) {
    return x >> r | x << (-r & 31);
}

uint32_t pcg32(void) {
    uint64_t x = state;
    unsigned count = static_cast<unsigned>(x >> 59);		// 59 = 64 - 5

    state = x * multiplier + increment;
    x ^= x >> 18;								// 18 = (64 - 27)/2
    return rotr32(static_cast<unsigned>(x >> 27), count);	// 27 = 32 - 5
}

void pcg32_init(uint64_t seed) {
    state = seed + increment;
    (void)pcg32();
}

static uint64_t mcg_state = 0xcafef00dd15ea5e5u;	// Must be odd

uint32_t pcg32_fast(void) {
    uint64_t x = mcg_state;
    unsigned count = static_cast<unsigned>(x >> 61);	// 61 = 64 - 3

    mcg_state = x * multiplier;
    x ^= x >> 22;
    return static_cast<unsigned>(x >> (22 + count));	// 22 = 32 - 3 - 7
}

void pcg32_fast_init(uint64_t seed) {
    mcg_state = 2 * seed + 1;
    (void)pcg32_fast();
}

#endif // PCG32_H
