// texture_segmentation.c
// CS136 - Texture-Based Image Segmentation
//
// Pipeline:
//   1. Load a NetPBM PGM/PPM image and use its grayscale intensity channel.
//   2. Compute local texture features in a square window:
//        local mean, local standard deviation, local gradient energy, entropy.
//   3. Normalize all features and cluster pixels with K-means.
//   4. Save a grayscale label mask, false-color region map, overlay, and
//      feature visualization images.
//
// Compile:
//   gcc texture_segmentation.c ../netpbm.c -o texture_segmentation -lm
//
// Run:
//   ./texture_segmentation ../input_images/coastline.pgm 4 9 30
//   ./texture_segmentation ../input_images/tahoe.pgm 5 11 40
//
// Arguments:
//   argv[1] input image path            default: ../input_images/coastline.pgm
//   argv[2] number of clusters K        default: 4
//   argv[3] local texture window size   default: 9
//   argv[4] K-means iterations          default: 30

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <sys/stat.h>
#include "../netpbm.h"

#define FEATURE_COUNT 4

typedef struct {
    double v[FEATURE_COUNT];
} Feature;

typedef struct {
    unsigned char r, g, b;
} RGB;

static int clampInt(int value, int lo, int hi) {
    if (value < lo) return lo;
    if (value > hi) return hi;
    return value;
}

static void writeMatrixAsPGM(Matrix mx, char *filename, int scale) {
    Image out = matrix2Image(mx, scale, 1.0);
    writeImage(out, filename);
    deleteImage(out);
}

static void getOutputPrefix(char *inputPath, char *prefix, int prefixSize) {
    char temp[256];
    char *base = strrchr(inputPath, '/');
    base = (base == NULL) ? inputPath : base + 1;

    snprintf(temp, sizeof(temp), "%s", base);
    char *dot = strrchr(temp, '.');
    if (dot != NULL) *dot = '\0';

    snprintf(prefix, prefixSize, "output_images/%s", temp);
}

static Matrix computeGradientMagnitude(Image img) {
    int H = img.height, W = img.width;
    Matrix grad = createMatrix(H, W);

    for (int y = 1; y < H - 1; y++) {
        for (int x = 1; x < W - 1; x++) {
            double gx =
                -img.map[y-1][x-1].i + img.map[y-1][x+1].i
                -2.0 * img.map[y][x-1].i + 2.0 * img.map[y][x+1].i
                -img.map[y+1][x-1].i + img.map[y+1][x+1].i;
            double gy =
                -img.map[y-1][x-1].i - 2.0 * img.map[y-1][x].i - img.map[y-1][x+1].i
                +img.map[y+1][x-1].i + 2.0 * img.map[y+1][x].i + img.map[y+1][x+1].i;
            grad.map[y][x] = sqrt(gx * gx + gy * gy);
        }
    }

    return grad;
}

static Feature *computeTextureFeatures(Image img, int windowSize,
                                       Matrix *meanImg, Matrix *stdImg,
                                       Matrix *gradImg, Matrix *entropyImg) {
    int H = img.height, W = img.width;
    int radius = windowSize / 2;
    Feature *features = (Feature *)malloc((size_t)H * W * sizeof(Feature));
    Matrix gradient = computeGradientMagnitude(img);

    *meanImg = createMatrix(H, W);
    *stdImg = createMatrix(H, W);
    *gradImg = createMatrix(H, W);
    *entropyImg = createMatrix(H, W);

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            double sum = 0.0;
            double sumSq = 0.0;
            double gradSum = 0.0;
            int hist[16] = {0};
            int count = 0;

            for (int dy = -radius; dy <= radius; dy++) {
                int yy = clampInt(y + dy, 0, H - 1);
                for (int dx = -radius; dx <= radius; dx++) {
                    int xx = clampInt(x + dx, 0, W - 1);
                    double val = (double)img.map[yy][xx].i;
                    sum += val;
                    sumSq += val * val;
                    gradSum += gradient.map[yy][xx];
                    hist[img.map[yy][xx].i / 16]++;
                    count++;
                }
            }

            double mean = sum / count;
            double variance = sumSq / count - mean * mean;
            double stddev = sqrt((variance > 0.0) ? variance : 0.0);
            double gradEnergy = gradSum / count;
            double entropy = 0.0;

            for (int b = 0; b < 16; b++) {
                if (hist[b] == 0) continue;
                double p = (double)hist[b] / count;
                entropy -= p * log(p) / log(2.0);
            }

            int idx = y * W + x;
            features[idx].v[0] = mean;
            features[idx].v[1] = stddev;
            features[idx].v[2] = gradEnergy;
            features[idx].v[3] = entropy;

            meanImg->map[y][x] = mean;
            stdImg->map[y][x] = stddev;
            gradImg->map[y][x] = gradEnergy;
            entropyImg->map[y][x] = entropy;
        }
    }

    deleteMatrix(gradient);
    return features;
}

static void normalizeFeatures(Feature *features, int count) {
    double minVal[FEATURE_COUNT];
    double maxVal[FEATURE_COUNT];

    for (int f = 0; f < FEATURE_COUNT; f++) {
        minVal[f] = DBL_MAX;
        maxVal[f] = -DBL_MAX;
    }

    for (int i = 0; i < count; i++) {
        for (int f = 0; f < FEATURE_COUNT; f++) {
            if (features[i].v[f] < minVal[f]) minVal[f] = features[i].v[f];
            if (features[i].v[f] > maxVal[f]) maxVal[f] = features[i].v[f];
        }
    }

    for (int i = 0; i < count; i++) {
        for (int f = 0; f < FEATURE_COUNT; f++) {
            double range = maxVal[f] - minVal[f];
            features[i].v[f] = (range > 1e-10) ? (features[i].v[f] - minVal[f]) / range : 0.0;
        }
    }
}

static double featureDistanceSq(Feature a, Feature b) {
    double total = 0.0;
    for (int f = 0; f < FEATURE_COUNT; f++) {
        double d = a.v[f] - b.v[f];
        total += d * d;
    }
    return total;
}

static int *kmeansCluster(Feature *features, int count, int K, int iterations) {
    int *labels = (int *)malloc((size_t)count * sizeof(int));
    Feature *centers = (Feature *)malloc((size_t)K * sizeof(Feature));
    int *clusterCounts = (int *)calloc((size_t)K, sizeof(int));

    for (int k = 0; k < K; k++) {
        int idx = (count <= 1 || K <= 1) ? 0 : (k * (count - 1)) / (K - 1);
        centers[k] = features[idx];
    }

    for (int i = 0; i < count; i++) labels[i] = -1;

    for (int iter = 0; iter < iterations; iter++) {
        int changed = 0;

        for (int i = 0; i < count; i++) {
            int best = 0;
            double bestDist = featureDistanceSq(features[i], centers[0]);
            for (int k = 1; k < K; k++) {
                double dist = featureDistanceSq(features[i], centers[k]);
                if (dist < bestDist) {
                    bestDist = dist;
                    best = k;
                }
            }
            if (labels[i] != best) {
                labels[i] = best;
                changed = 1;
            }
        }

        for (int k = 0; k < K; k++) {
            clusterCounts[k] = 0;
            for (int f = 0; f < FEATURE_COUNT; f++) centers[k].v[f] = 0.0;
        }

        for (int i = 0; i < count; i++) {
            int k = labels[i];
            clusterCounts[k]++;
            for (int f = 0; f < FEATURE_COUNT; f++) centers[k].v[f] += features[i].v[f];
        }

        for (int k = 0; k < K; k++) {
            if (clusterCounts[k] == 0) {
                centers[k] = features[(k * 9973) % count];
                continue;
            }
            for (int f = 0; f < FEATURE_COUNT; f++) centers[k].v[f] /= clusterCounts[k];
        }

        printf("  K-means iteration %d/%d%s\n",
               iter + 1, iterations, changed ? "" : " converged");
        if (!changed) break;
    }

    for (int k = 0; k < K; k++)
        printf("  Cluster %d pixels: %d\n", k, clusterCounts[k]);

    free(centers);
    free(clusterCounts);
    return labels;
}

static Image makeLabelMask(int *labels, int H, int W, int K) {
    Image mask = createImage(H, W);
    int denom = (K > 1) ? K - 1 : 1;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            unsigned char v = (unsigned char)(labels[y * W + x] * 255 / denom);
            mask.map[y][x].i = mask.map[y][x].r =
            mask.map[y][x].g = mask.map[y][x].b = v;
        }
    }

    return mask;
}

static Image makeColorSegmentation(int *labels, int H, int W) {
    RGB palette[] = {
        {230,  57,  70}, { 42, 157, 143}, { 69, 123, 157}, {244, 162,  97},
        {131,  56, 236}, {255, 190,  11}, {  6, 214, 160}, { 17, 138, 178}
    };
    int paletteSize = (int)(sizeof(palette) / sizeof(palette[0]));
    Image out = createImage(H, W);

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            RGB c = palette[labels[y * W + x] % paletteSize];
            out.map[y][x].r = c.r;
            out.map[y][x].g = c.g;
            out.map[y][x].b = c.b;
            out.map[y][x].i = (unsigned char)((c.r + c.g + c.b) / 3);
        }
    }

    return out;
}

static Image makeOverlay(Image img, Image colorSeg) {
    Image out = createImage(img.height, img.width);

    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            out.map[y][x].r = (unsigned char)(0.55 * img.map[y][x].r + 0.45 * colorSeg.map[y][x].r);
            out.map[y][x].g = (unsigned char)(0.55 * img.map[y][x].g + 0.45 * colorSeg.map[y][x].g);
            out.map[y][x].b = (unsigned char)(0.55 * img.map[y][x].b + 0.45 * colorSeg.map[y][x].b);
            out.map[y][x].i = (unsigned char)((out.map[y][x].r + out.map[y][x].g + out.map[y][x].b) / 3);
        }
    }

    return out;
}

int main(int argc, char *argv[]) {
    char *inputFile = (argc > 1) ? argv[1] : "../input_images/coastline.pgm";
    int K = (argc > 2) ? atoi(argv[2]) : 4;
    int windowSize = (argc > 3) ? atoi(argv[3]) : 9;
    int iterations = (argc > 4) ? atoi(argv[4]) : 30;

    if (K < 2) K = 2;
    if (K > 8) K = 8;
    if (windowSize < 3) windowSize = 3;
    if (windowSize % 2 == 0) windowSize++;
    if (iterations < 1) iterations = 1;

    printf("Reading '%s' ...\n", inputFile);
    Image img = readImage(inputFile);
    printf("Image: %d x %d  K=%d  window=%d  iterations=%d\n",
           img.width, img.height, K, windowSize, iterations);

#ifdef _WIN32
    mkdir("output_images");
#else
    mkdir("output_images", 0755);
#endif

    Matrix meanImg, stdImg, gradImg, entropyImg;
    Feature *features = computeTextureFeatures(img, windowSize,
                                               &meanImg, &stdImg,
                                               &gradImg, &entropyImg);
    int pixelCount = img.height * img.width;
    normalizeFeatures(features, pixelCount);
    int *labels = kmeansCluster(features, pixelCount, K, iterations);

    Image labelMask = makeLabelMask(labels, img.height, img.width, K);
    Image colorSeg = makeColorSegmentation(labels, img.height, img.width);
    Image overlay = makeOverlay(img, colorSeg);

    char prefix[320];
    char labelFile[384], regionFile[384], overlayFile[384];
    char meanFile[384], stdFile[384], gradFile[384], entropyFile[384];
    getOutputPrefix(inputFile, prefix, sizeof(prefix));
    snprintf(labelFile, sizeof(labelFile), "%s_texture_labels.pgm", prefix);
    snprintf(regionFile, sizeof(regionFile), "%s_texture_regions.ppm", prefix);
    snprintf(overlayFile, sizeof(overlayFile), "%s_texture_overlay.ppm", prefix);
    snprintf(meanFile, sizeof(meanFile), "%s_feature_mean.pgm", prefix);
    snprintf(stdFile, sizeof(stdFile), "%s_feature_stddev.pgm", prefix);
    snprintf(gradFile, sizeof(gradFile), "%s_feature_gradient_energy.pgm", prefix);
    snprintf(entropyFile, sizeof(entropyFile), "%s_feature_entropy.pgm", prefix);

    writeImage(labelMask, labelFile);
    writeImage(colorSeg, regionFile);
    writeImage(overlay, overlayFile);
    writeMatrixAsPGM(meanImg, meanFile, 1);
    writeMatrixAsPGM(stdImg, stdFile, 1);
    writeMatrixAsPGM(gradImg, gradFile, 1);
    writeMatrixAsPGM(entropyImg, entropyFile, 1);

    printf("Label mask        -> %s\n", labelFile);
    printf("Color regions     -> %s\n", regionFile);
    printf("Overlay           -> %s\n", overlayFile);
    printf("Feature images    -> %s_feature_*.pgm\n", prefix);

    free(features);
    free(labels);
    deleteMatrix(meanImg);
    deleteMatrix(stdImg);
    deleteMatrix(gradImg);
    deleteMatrix(entropyImg);
    deleteImage(img);
    deleteImage(labelMask);
    deleteImage(colorSeg);
    deleteImage(overlay);

    printf("Done.\n");
    return 0;
}
