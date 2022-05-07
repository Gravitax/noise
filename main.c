#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define couleur(param) printf("\033[%sm",param)

void    set_color(int color) {
    if (color == -1)
        couleur("0");
    else if (color > 35)
        couleur("37");
    else if (color > 30)
        couleur("30");
    else if (color > 25)
        couleur("31");
    else if (color > 20)
        couleur("33");
    else if (color > 10)
        couleur("32");
    else
        couleur("34");
}

void    perlin_noise1d(int count, char *seed, int octaves, float bias, char *output) {
    int     x = -1;
    int     o, pitch, sample1, sample2;
    float   noise, blend, sample, scale, scaleAcc;

    while (++x < count) {
        noise = 0;
        scale = 1;
        scaleAcc = 0;
        o = -1;
        while (++o < octaves) {
            pitch = count >> o;
            pitch = pitch < 1 ? 1 : pitch;

            sample1 = (x / pitch) * pitch;
            sample2 = (sample1 + pitch) % count;

            blend = (float)(x - sample1) / (float)pitch;

            sample = (1.0f - blend) * seed[sample1] + blend * seed[sample2];

            noise += sample * scale;
            scaleAcc += scale;
            scale = scale / bias;
        }
        output[x] = (char)(noise / scaleAcc);
    }
}

void    perlin_noise2d(int width, int height, char *seed, int octaves, float bias, char *output) {
    int     x = -1;
    int     y = -1;
    int     o, pitch, sampleX1, sampleX2, sampleY1, sampleY2;
    float   noise, blendX, blendY, sampleT, sampleB, scale, scaleAcc;

    while (++x < width) {
        y = -1;
        while (++y < height) {
            noise = 0;
            scale = 1;
            scaleAcc = 0;
            o = -1;
            while (++o < octaves) {
                pitch = width >> o;
                pitch = pitch < 1 ? 1 : pitch;

                sampleX1 = (x / pitch) * pitch;
                sampleY1 = (y / pitch) * pitch;
                sampleX2 = (sampleX1 + pitch) % width;
                sampleY2 = (sampleY1 + pitch) % width;

                blendX = (float)(x - sampleX1) / (float)pitch;
                blendY = (float)(y - sampleY1) / (float)pitch;

                sampleT = (1.0f - blendX) * seed[sampleY1 * width + sampleX1] + blendX * seed[sampleY1 * width + sampleX2];
                sampleB = (1.0f - blendX) * seed[sampleY2 * width + sampleX1] + blendX * seed[sampleY2 * width + sampleX2];

                noise += (blendY * (sampleB - sampleT) + sampleT) * scale;
                scaleAcc += scale;
                scale = scale / bias;
            }
            output[y * width + x] = (char)(noise / scaleAcc);
        }
    }
}

void    javidx9(const int width, const int height, const int length) {
    const int   outsize = 64;
    char        seed1d[outsize];
    char        noise1d[outsize];
    char        seed2d[length];
    char        noise2d[length];
    int         i = -1;

    while (++i < outsize)
        seed1d[i] = (char)(1 + rand() % (outsize / 2));
    i = -1;
    while (++i < length)
        seed2d[i] = (char)(1 + rand() % width);

    printf("========================\n");

    perlin_noise1d(outsize, seed1d, 4, 2, noise1d);
    perlin_noise2d(width, height, seed2d, 4.5, 1.5, noise2d);

    printf("------------------------\n");

    i = -1;
    while (++i < outsize) {
        set_color((int)noise1d[i]);
        printf("[]");
        set_color(-1);
    }
    printf("\n------------------------\n");
    i = -1;
    while (++i < length) {
        set_color((int)noise2d[i]);
        printf("[]");
        if (i > 0 && (i + 1) % width == 0)
            printf("\n");
        set_color(-1);
    }
    printf("------------------------\n");
}

void    get_biome(float e) {
    if (e == -1) // reset
        couleur("0");
    else if (e < 0) // deep water
        couleur("0;34");
    else if (e < 0.1) // water
        couleur("1;34");
    else if (e < 0.2) // sand
        couleur("0;33");
    else if (e < 0.5) // forest
        couleur("0;32");
    else if (e < 0.6) // deep forest
        couleur("1;32");
    else if (e < 0.7) // red rock
        couleur("0;31");
    else if (e < 0.8) // grey rock
        couleur("1;30");
    else if (e < 0.9) // dark rock
        couleur("0;30");
    else // snow
        couleur("0;37");
}

/* Function to linearly interpolate between a0 and a1
 * Weight w should be in the range [0.0, 1.0]
 */
float interpolate(float a0, float a1, float w) {
    // You may want clamping by inserting:
    // if (w < 0.0) return a0;
    // if (w > 1.0) return a1;

    return (a1 - a0) * w + a0;
    // Use this cubic interpolation [[Smoothstep]] instead, for a smooth appearance:
    return (a1 - a0) * (3.0 - w * 2.0) * w * w + a0;

    // Use [[Smootherstep]] for an even smoother result with a second derivative equal to zero on boundaries:
    return (a1 - a0) * ((w * (w * 6.0 - 15.0) + 10.0) * w * w * w) + a0;
    
}

typedef struct {
    float x, y;
} vector2;

/* Create pseudorandom direction vector
 */
vector2 randomGradient(int ix, int iy) {
    // No precomputed gradients mean this works for any number of grid coordinates
    const unsigned w = 8 * sizeof(unsigned);
    const unsigned s = w / 2; // rotation width
    unsigned a = ix, b = iy;
    a *= 3284157443; b ^= a << s | a >> w-s;
    b *= 1911520717; a ^= b << s | b >> w-s;
    a *= 2048419325;
    float random = a * (3.14159265 / ~(~0u >> 1)); // in [0, 2*Pi]
    vector2 v;
    v.x = cos(random); v.y = sin(random);
    return v;
}

// Computes the dot product of the distance and gradient vectors.
float dotGridGradient(int ix, int iy, float x, float y) {
    // Get gradient from integer coordinates
    vector2 gradient = randomGradient(ix, iy);

    // Compute the distance vector
    float dx = x - (float)ix;
    float dy = y - (float)iy;

    // Compute the dot-product
    return (dx*gradient.x + dy*gradient.y);
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
    return value;
}

int     main() {
    
    const int   width = 64;
    const int   height = 32;
    const int   length = width * height;

    // javidx9(width, height, length);

    float   elevation[length];
    float   nx, ny, e, factor, amplitude, amplitude_sum;

    const float frequency = 2.0f;
    const float wavelength = (float)length / frequency;
    const float redistribution = 4.0f;
    const int   octaves = 3;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            nx = (float)x / width - 0.5f;
            ny = (float)y / height - 0.5f;
            nx *= frequency;
            ny *= frequency;
            // nx /= wavelength;
            // ny /= wavelength;
            e = 0;
            amplitude_sum = 0;
            for (int o = 1; o <= octaves; o++) {
                factor = o * 2;
                amplitude = 1 / factor;
                e += amplitude * perlin(nx * factor, ny * factor);
                amplitude_sum += amplitude;
            }
            if (amplitude_sum > 0)
                e /= amplitude_sum;
            elevation[y * width + x] = e * redistribution;
        }
    }
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            get_biome(elevation[y * width + x]);
            // printf("[%.2f]", elevation[y * width + x]);
            printf("[]");
            get_biome(-1);
        }
        printf("\n");
    }
    return (0);
}
