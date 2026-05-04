// robustness_analysis.c
// CS136 - Robustness Analysis
//
// Applies real-world distortions to sample images, then runs:
//   Gaussian filtering, Sobel edges, Canny edges, and texture K-means segmentation.
//
// Compile:
//   gcc robustness_analysis.c ../netpbm.c -o robustness_analysis -lm
//
// Run:
//   ./robustness_analysis

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include "../netpbm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FEATURES 4

typedef struct { double v[FEATURES]; } Feature;
typedef struct { unsigned char r, g, b; } RGB;

static unsigned int rngState = 123456789u;

static int clampInt(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static unsigned char clampByte(double v) {
    if (v < 0.0) return 0;
    if (v > 255.0) return 255;
    return (unsigned char)(v + 0.5);
}

static double randUnit(void) {
    rngState = 1664525u * rngState + 1013904223u;
    return ((double)(rngState & 0x00ffffffu) + 1.0) / 16777217.0;
}

static double randGaussian(void) {
    double u1 = randUnit();
    double u2 = randUnit();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static void setGray(Image img, int y, int x, unsigned char v) {
    img.map[y][x].i = img.map[y][x].r = img.map[y][x].g = img.map[y][x].b = v;
}

static Image cloneImage(Image img) {
    Image out = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            out.map[y][x] = img.map[y][x];
    return out;
}

static void convolve(Matrix src, Matrix kernel, Matrix dst) {
    int ay = kernel.height / 2;
    int ax = kernel.width / 2;
    for (int y = 0; y < src.height; y++) {
        for (int x = 0; x < src.width; x++) {
            double sum = 0.0;
            for (int ky = 0; ky < kernel.height; ky++) {
                int yy = clampInt(y + ky - ay, 0, src.height - 1);
                for (int kx = 0; kx < kernel.width; kx++) {
                    int xx = clampInt(x + kx - ax, 0, src.width - 1);
                    sum += src.map[yy][xx] * kernel.map[ky][kx];
                }
            }
            dst.map[y][x] = sum;
        }
    }
}

static Matrix imageToIntensity(Image img) {
    Matrix mx = createMatrix(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            mx.map[y][x] = img.map[y][x].i;
    return mx;
}

static Image matrixToGray(Matrix mx, int scale) {
    Image out = matrix2Image(mx, scale, 1.0);
    return out;
}

static Matrix gaussianKernel5(void) {
    double data[5][5] = {
        { 2,  4,  5,  4,  2},
        { 4,  9, 12,  9,  4},
        { 5, 12, 15, 12,  5},
        { 4,  9, 12,  9,  4},
        { 2,  4,  5,  4,  2}
    };
    Matrix k = createMatrix(5, 5);
    double total = 0.0;
    for (int y = 0; y < 5; y++)
        for (int x = 0; x < 5; x++)
            total += data[y][x];
    for (int y = 0; y < 5; y++)
        for (int x = 0; x < 5; x++)
            k.map[y][x] = data[y][x] / total;
    return k;
}

static Image gaussianFilter(Image img) {
    Matrix src = imageToIntensity(img);
    Matrix dst = createMatrix(img.height, img.width);
    Matrix k = gaussianKernel5();
    convolve(src, k, dst);
    Image out = matrixToGray(dst, 0);
    deleteMatrix(src);
    deleteMatrix(dst);
    deleteMatrix(k);
    return out;
}

static Image addGaussianNoise(Image img, double sigma) {
    Image out = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            int noise = (int)(sigma * randGaussian());
            out.map[y][x].r = clampByte(img.map[y][x].r + noise);
            out.map[y][x].g = clampByte(img.map[y][x].g + noise);
            out.map[y][x].b = clampByte(img.map[y][x].b + noise);
            out.map[y][x].i = clampByte(img.map[y][x].i + noise);
        }
    }
    return out;
}

static Image lowContrast(Image img, double factor) {
    Image out = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            out.map[y][x].r = clampByte(128.0 + factor * (img.map[y][x].r - 128.0));
            out.map[y][x].g = clampByte(128.0 + factor * (img.map[y][x].g - 128.0));
            out.map[y][x].b = clampByte(128.0 + factor * (img.map[y][x].b - 128.0));
            out.map[y][x].i = clampByte(128.0 + factor * (img.map[y][x].i - 128.0));
        }
    }
    return out;
}

static void sobelMagnitude(Image img, Matrix mag) {
    for (int y = 1; y < img.height - 1; y++) {
        for (int x = 1; x < img.width - 1; x++) {
            double gx =
                -img.map[y-1][x-1].i + img.map[y-1][x+1].i
                -2.0 * img.map[y][x-1].i + 2.0 * img.map[y][x+1].i
                -img.map[y+1][x-1].i + img.map[y+1][x+1].i;
            double gy =
                -img.map[y-1][x-1].i - 2.0 * img.map[y-1][x].i - img.map[y-1][x+1].i
                +img.map[y+1][x-1].i + 2.0 * img.map[y+1][x].i + img.map[y+1][x+1].i;
            mag.map[y][x] = sqrt(gx * gx + gy * gy);
        }
    }
}

static Image sobelEdges(Image img) {
    Matrix mag = createMatrix(img.height, img.width);
    sobelMagnitude(img, mag);
    Image out = matrixToGray(mag, 1);
    deleteMatrix(mag);
    return out;
}

static Image cannyEdges(Image img, double highRatio, double lowRatio) {
    int H = img.height, W = img.width;
    Image smoothedImg = gaussianFilter(img);
    Matrix mag = createMatrix(H, W);
    Matrix dir = createMatrix(H, W);

    for (int y = 1; y < H - 1; y++) {
        for (int x = 1; x < W - 1; x++) {
            double gx =
                -smoothedImg.map[y-1][x-1].i + smoothedImg.map[y-1][x+1].i
                -2.0 * smoothedImg.map[y][x-1].i + 2.0 * smoothedImg.map[y][x+1].i
                -smoothedImg.map[y+1][x-1].i + smoothedImg.map[y+1][x+1].i;
            double gy =
                -smoothedImg.map[y-1][x-1].i - 2.0 * smoothedImg.map[y-1][x].i - smoothedImg.map[y-1][x+1].i
                +smoothedImg.map[y+1][x-1].i + 2.0 * smoothedImg.map[y+1][x].i + smoothedImg.map[y+1][x+1].i;
            mag.map[y][x] = sqrt(gx * gx + gy * gy);
            dir.map[y][x] = atan2(gy, gx);
        }
    }

    Matrix nms = createMatrix(H, W);
    for (int y = 1; y < H - 1; y++) {
        for (int x = 1; x < W - 1; x++) {
            double angle = dir.map[y][x] * 180.0 / M_PI;
            if (angle < 0.0) angle += 180.0;
            double n1, n2;
            if      (angle < 22.5 || angle >= 157.5) { n1 = mag.map[y][x-1]; n2 = mag.map[y][x+1]; }
            else if (angle < 67.5)                   { n1 = mag.map[y-1][x+1]; n2 = mag.map[y+1][x-1]; }
            else if (angle < 112.5)                  { n1 = mag.map[y-1][x]; n2 = mag.map[y+1][x]; }
            else                                     { n1 = mag.map[y-1][x-1]; n2 = mag.map[y+1][x+1]; }
            nms.map[y][x] = (mag.map[y][x] >= n1 && mag.map[y][x] >= n2) ? mag.map[y][x] : 0.0;
        }
    }

    double maxMag = 0.0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            if (nms.map[y][x] > maxMag) maxMag = nms.map[y][x];

    double high = maxMag * highRatio;
    double low = high * lowRatio;
    Matrix edges = createMatrix(H, W);
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            if (nms.map[y][x] >= high) edges.map[y][x] = 255.0;
            else if (nms.map[y][x] >= low) edges.map[y][x] = 128.0;
            else edges.map[y][x] = 0.0;
        }
    }

    int changed = 1;
    while (changed) {
        changed = 0;
        for (int y = 1; y < H - 1; y++) {
            for (int x = 1; x < W - 1; x++) {
                if (edges.map[y][x] != 128.0) continue;
                int strong = 0;
                for (int dy = -1; dy <= 1 && !strong; dy++)
                    for (int dx = -1; dx <= 1 && !strong; dx++)
                        if (edges.map[y+dy][x+dx] == 255.0) strong = 1;
                if (strong) { edges.map[y][x] = 255.0; changed = 1; }
            }
        }
    }

    Image out = createImage(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            setGray(out, y, x, (edges.map[y][x] == 255.0) ? 255 : 0);

    deleteImage(smoothedImg);
    deleteMatrix(mag); deleteMatrix(dir); deleteMatrix(nms); deleteMatrix(edges);
    return out;
}

static long countEdges(Image edgeImg, int threshold) {
    long count = 0;
    for (int y = 0; y < edgeImg.height; y++)
        for (int x = 0; x < edgeImg.width; x++)
            if (edgeImg.map[y][x].i > threshold) count++;
    return count;
}

static Feature *textureFeatures(Image img, int window) {
    int H = img.height, W = img.width, radius = window / 2;
    Feature *features = (Feature *)malloc((size_t)H * W * sizeof(Feature));
    Matrix grad = createMatrix(H, W);
    sobelMagnitude(img, grad);

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            double sum = 0.0, sumSq = 0.0, gradSum = 0.0, entropy = 0.0;
            int hist[16] = {0};
            int n = 0;
            for (int dy = -radius; dy <= radius; dy++) {
                int yy = clampInt(y + dy, 0, H - 1);
                for (int dx = -radius; dx <= radius; dx++) {
                    int xx = clampInt(x + dx, 0, W - 1);
                    double v = img.map[yy][xx].i;
                    sum += v;
                    sumSq += v * v;
                    gradSum += grad.map[yy][xx];
                    hist[img.map[yy][xx].i / 16]++;
                    n++;
                }
            }
            for (int b = 0; b < 16; b++) {
                if (hist[b] == 0) continue;
                double p = (double)hist[b] / n;
                entropy -= p * log(p) / log(2.0);
            }
            double mean = sum / n;
            double variance = sumSq / n - mean * mean;
            int idx = y * W + x;
            features[idx].v[0] = mean;
            features[idx].v[1] = sqrt((variance > 0.0) ? variance : 0.0);
            features[idx].v[2] = gradSum / n;
            features[idx].v[3] = entropy;
        }
    }

    deleteMatrix(grad);
    return features;
}

static void normalizeFeatures(Feature *features, int count) {
    double mn[FEATURES], mx[FEATURES];
    for (int f = 0; f < FEATURES; f++) { mn[f] = DBL_MAX; mx[f] = -DBL_MAX; }
    for (int i = 0; i < count; i++)
        for (int f = 0; f < FEATURES; f++) {
            if (features[i].v[f] < mn[f]) mn[f] = features[i].v[f];
            if (features[i].v[f] > mx[f]) mx[f] = features[i].v[f];
        }
    for (int i = 0; i < count; i++)
        for (int f = 0; f < FEATURES; f++)
            features[i].v[f] = (mx[f] - mn[f] > 1e-10) ? (features[i].v[f] - mn[f]) / (mx[f] - mn[f]) : 0.0;
}

static double distSq(Feature a, Feature b) {
    double d = 0.0;
    for (int f = 0; f < FEATURES; f++) {
        double x = a.v[f] - b.v[f];
        d += x * x;
    }
    return d;
}

static int *kmeans(Feature *features, int count, int K, int iterations, int *clusterCounts) {
    int *labels = (int *)malloc((size_t)count * sizeof(int));
    Feature *centers = (Feature *)malloc((size_t)K * sizeof(Feature));
    for (int k = 0; k < K; k++) centers[k] = features[(k * (count - 1)) / (K - 1)];
    for (int i = 0; i < count; i++) labels[i] = -1;

    for (int it = 0; it < iterations; it++) {
        for (int i = 0; i < count; i++) {
            int best = 0;
            double bestD = distSq(features[i], centers[0]);
            for (int k = 1; k < K; k++) {
                double d = distSq(features[i], centers[k]);
                if (d < bestD) { bestD = d; best = k; }
            }
            labels[i] = best;
        }

        for (int k = 0; k < K; k++) {
            clusterCounts[k] = 0;
            for (int f = 0; f < FEATURES; f++) centers[k].v[f] = 0.0;
        }
        for (int i = 0; i < count; i++) {
            int k = labels[i];
            clusterCounts[k]++;
            for (int f = 0; f < FEATURES; f++) centers[k].v[f] += features[i].v[f];
        }
        for (int k = 0; k < K; k++) {
            if (clusterCounts[k] == 0) centers[k] = features[(k * 9973) % count];
            else for (int f = 0; f < FEATURES; f++) centers[k].v[f] /= clusterCounts[k];
        }
    }
    free(centers);
    return labels;
}

static Image textureSegmentation(Image img, int K, int window, int iterations, int *clusterCounts) {
    int count = img.height * img.width;
    Feature *features = textureFeatures(img, window);
    normalizeFeatures(features, count);
    int *labels = kmeans(features, count, K, iterations, clusterCounts);
    RGB palette[] = {
        {230, 57, 70}, {42, 157, 143}, {69, 123, 157}, {244, 162, 97},
        {131, 56, 236}, {255, 190, 11}
    };
    Image out = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            RGB c = palette[labels[y * img.width + x] % 6];
            out.map[y][x].r = c.r; out.map[y][x].g = c.g; out.map[y][x].b = c.b;
            out.map[y][x].i = (unsigned char)((c.r + c.g + c.b) / 3);
        }
    }
    free(features);
    free(labels);
    return out;
}

static void saveName(char *dst, int size, const char *imageName, const char *caseName, const char *suffix) {
    snprintf(dst, size, "output_images/%s_%s_%s", imageName, caseName, suffix);
}

static void runCase(FILE *csv, Image img, const char *imageName, const char *caseName) {
    char path[384];
    int total = img.height * img.width;
    printf("Processing %s / %s\n", imageName, caseName);

    saveName(path, sizeof(path), imageName, caseName, "distorted.ppm");
    writeImage(img, path);

    Image gauss = gaussianFilter(img);
    saveName(path, sizeof(path), imageName, caseName, "gaussian.ppm");
    writeImage(gauss, path);

    Image sobelRaw = sobelEdges(img);
    Image sobelGauss = sobelEdges(gauss);
    Image cannyRaw = cannyEdges(img, 0.15, 0.50);
    Image cannyGauss = cannyEdges(gauss, 0.15, 0.50);

    saveName(path, sizeof(path), imageName, caseName, "sobel_raw.pgm");
    writeImage(sobelRaw, path);
    saveName(path, sizeof(path), imageName, caseName, "sobel_after_gaussian.pgm");
    writeImage(sobelGauss, path);
    saveName(path, sizeof(path), imageName, caseName, "canny_raw.pgm");
    writeImage(cannyRaw, path);
    saveName(path, sizeof(path), imageName, caseName, "canny_after_gaussian.pgm");
    writeImage(cannyGauss, path);

    int clusterCounts[4] = {0, 0, 0, 0};
    Image seg = textureSegmentation(img, 4, 7, 12, clusterCounts);
    saveName(path, sizeof(path), imageName, caseName, "texture_segmentation.ppm");
    writeImage(seg, path);

    long sobelRawCount = countEdges(sobelRaw, 128);
    long sobelGaussCount = countEdges(sobelGauss, 128);
    long cannyRawCount = countEdges(cannyRaw, 128);
    long cannyGaussCount = countEdges(cannyGauss, 128);

    fprintf(csv, "%s,%s,%ld,%.3f,%ld,%.3f,%ld,%.3f,%ld,%.3f,%d,%d,%d,%d\n",
            imageName, caseName,
            sobelRawCount, 100.0 * sobelRawCount / total,
            sobelGaussCount, 100.0 * sobelGaussCount / total,
            cannyRawCount, 100.0 * cannyRawCount / total,
            cannyGaussCount, 100.0 * cannyGaussCount / total,
            clusterCounts[0], clusterCounts[1], clusterCounts[2], clusterCounts[3]);

    deleteImage(gauss);
    deleteImage(sobelRaw); deleteImage(sobelGauss);
    deleteImage(cannyRaw); deleteImage(cannyGauss);
    deleteImage(seg);
}

static void processImage(FILE *csv, const char *path, const char *name) {
    Image original = readImage((char *)path);
    Image noisy = addGaussianNoise(original, 24.0);
    Image blurred = gaussianFilter(original);
    Image contrast = lowContrast(original, 0.45);

    runCase(csv, original, name, "original");
    runCase(csv, noisy, name, "gaussian_noise");
    runCase(csv, blurred, name, "blur");
    runCase(csv, contrast, name, "low_contrast");

    deleteImage(original);
    deleteImage(noisy);
    deleteImage(blurred);
    deleteImage(contrast);
}

int main(void) {
#ifdef _WIN32
    mkdir("output_images");
#else
    mkdir("output_images", 0755);
#endif

    FILE *csv = fopen("output_images/metrics.csv", "w");
    if (csv == NULL) {
        fprintf(stderr, "Could not write output_images/metrics.csv\n");
        return 1;
    }
    fprintf(csv, "image,case,sobel_raw_edges,sobel_raw_pct,sobel_after_gaussian_edges,sobel_after_gaussian_pct,canny_raw_edges,canny_raw_pct,canny_after_gaussian_edges,canny_after_gaussian_pct,seg_cluster0,seg_cluster1,seg_cluster2,seg_cluster3\n");

    processImage(csv, "../input_images/coastline.ppm", "coastline");
    processImage(csv, "../input_images/tahoe.ppm", "tahoe");

    fclose(csv);
    printf("Done. Metrics -> output_images/metrics.csv\n");
    return 0;
}
