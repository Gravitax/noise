#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WIDTH 128
#define HEIGHT 64
#define LENGTH WIDTH * HEIGHT

#define color(param) printf("\033[%sm",param)


void    get_biome(float e) {
    if (e == -1) // reset
        color("0");
    else if (e < -0.25f) // deep water
        color("0;34");
    else if (e < 0.2f) // water
        color("1;34");
    else if (e < 0.3f) // sand
        color("0;33");
    else if (e < 0.5f) // forest
        color("0;32");
    else if (e < 0.6f) // deep forest
        color("1;32");
    else if (e < 0.7f) // red rock
        color("0;31");
    else if (e < 0.9f) // grey rock
        color("1;30");
    else if (e < 1.0f) // dark rock
        color("0;30");
    else // snow
        color("0;37");
}

/* Function to linearly interpolate between a0 and a1
 * Weight w should be in the range [0.0, 1.0]
 */
float interpolate(float a0, float a1, float w) {
    // You may want clamping by inserting:
    // if (w < 0.0) return a0;
    // if (w > 1.0) return a1;

    // return ((a1 - a0) * w + a0);
    // Use this cubic interpolation [[Smoothstep]] instead, for a smooth appearance:
    return ((a1 - a0) * (3.0f - w * 2.0f) * w * w + a0);

    // Use [[Smootherstep]] for an even smoother result with a second derivative equal to zero on boundaries:
    return ((a1 - a0) * ((w * (w * 6.0f - 15.0f) + 10.0f) * w * w * w) + a0);
    
}

typedef struct {
    float x, y;
} vector2;

/* Create pseudorandom direction vector
 */
vector2 randomGradient(int ix, int iy) {
    // No precomputed gradients mean this works for any number of grid coordinates
    const unsigned w = 8 * sizeof(unsigned);
    const unsigned s = w / 2; // rotation WIDTH
    unsigned a = ix, b = iy;

    a *= 3284157443; b ^= a << s | a >> (w - s);
    b *= 1911520717; a ^= b << s | b >> (w - s);
    a *= 2048419325;

    float random = a * (3.14159265 / ~(~0u >> 1)); // in [0, 2*Pi]
    vector2 v;

    v.x = cos(random); v.y = sin(random);
    return (v);
}

// Computes the dot product of the distance and gradient vectors.
float dotGridGradient(int ix, int iy, float x, float y) {
    // Get gradient from integer coordinates
    vector2 gradient = randomGradient(ix, iy);

    // Compute the distance vector
    float dx = x - (float)ix;
    float dy = y - (float)iy;

    // Compute the dot-product
    return (dx * gradient.x + dy * gradient.y);
}

// Compute Perlin noise at coordinates x, y
float perlin(float x, float y) {
    // Determine grid cell coordinates
    int x0 = (int)floor(x);
    int x1 = x0 + 1;
    int y0 = (int)floor(y);
    int y1 = y0 + 1;

    // Determine interpolation weights
    // Could also use higher order polynomial/s-curve here
    float sx = x - (float)x0;
    float sy = y - (float)y0;

    // Interpolate between grid point gradients
    float n0, n1, ix0, ix1, value;

    n0 = dotGridGradient(x0, y0, x, y);
    n1 = dotGridGradient(x1, y0, x, y);
    ix0 = interpolate(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);
    ix1 = interpolate(n0, n1, sx);

    value = interpolate(ix0, ix1, sy);
    return (value);
}

float	noise(int x, int y) {
	const float frequency = 2.0f;
    const float redistribution = 4.0f;
    const int   octaves = 3;

	float   nx, ny, e, factor, amplitude, amplitude_sum;

	nx = (float)x / WIDTH - 0.5f;
    ny = (float)y / HEIGHT - 0.5f;
    nx *= frequency;
    ny *= frequency;
    e = 0;
    amplitude_sum = 0;
    for (int o = 1; o < octaves + 1; o++) {
        factor = o * 2;
        amplitude = 1 / factor;
        e += amplitude * perlin(nx * factor, ny * factor);
        amplitude_sum += amplitude;
    }
    if (amplitude_sum > 0)
        e /= amplitude_sum;
    return (e * redistribution);
}

void	print_noise(float e) {
	get_biome(e);
    // printf("[%.2f]", e);
    printf("[]");
    get_biome(-1);
}

int     main() {
    float   elevation[LENGTH];
	int		pos;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
			pos = y * WIDTH + x;
            elevation[pos] = noise(x, y);
			print_noise(elevation[pos]);
        }
		printf("\n");
    }
    return (0);
}
