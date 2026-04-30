// gaussian.c
// CS136 – Gaussian / Smoothing Filter
//
// Contains:
//   smoothing_filter() uniform (averaging) kernel convolution
//   median_filter() median over an arbitrary kernel-shaped window
//   gaussian_filter() true Gaussian kernel convolution
//
// Compile:
//   gcc gaussian.c ../netpbm.c -o gaussian -lm
// Run:
//   gaussian.exe ../input_images/coastline.ppm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "../netpbm.h"

static void matrixToImageScaled(Matrix *m, Image *img) {
    int h = m->height, w = m->width;
    double minVal = m->map[0][0], maxVal = m->map[0][0];
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            if (m->map[y][x] < minVal) minVal = m->map[y][x];
            if (m->map[y][x] > maxVal) maxVal = m->map[y][x];
        }
    double range = (maxVal - minVal > 1e-10) ? (maxVal - minVal) : 1.0;
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            unsigned char v = (unsigned char)(
                ((m->map[y][x] - minVal) / range) * 255.0 + 0.5);
            img->map[y][x].i = img->map[y][x].r =
            img->map[y][x].g = img->map[y][x].b = v;
        }
}

static void _convolve(Matrix *m1, Matrix *m2, Matrix *result) {
    int H = m1->height, W = m1->width;
    int anchorY = (m2->height - 1) / 2;
    int anchorX = (m2->width  - 1) / 2;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int top   = y - anchorY;
            int bot   = y + (m2->height - 1 - anchorY);
            int left  = x - anchorX;
            int right = x + (m2->width  - 1 - anchorX);
            if (top < 0 || bot >= H || left < 0 || right >= W) {
                result->map[y][x] = 0.0;
                continue;
            }
            double sum = 0.0;
            for (int ky = 0; ky < m2->height; ky++)
                for (int kx = 0; kx < m2->width; kx++)
                    sum += m1->map[top + ky][left + kx] * m2->map[ky][kx];
            result->map[y][x] = sum;
        }
    }
}

static int cmpDouble(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

void smoothing_filter(Matrix m1, Matrix m2, Matrix *result) {
    int H = m1.height, W = m1.width;
    int anchorY = (m2.height - 1) / 2, anchorX = (m2.width - 1) / 2;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int top   = y - anchorY, bot   = y + (m2.height - 1 - anchorY);
            int left  = x - anchorX, right = x + (m2.width  - 1 - anchorX);
            if (top < 0 || bot >= H || left < 0 || right >= W) {
                result->map[y][x] = 0.0;
                continue;
            }
            double sum = 0.0;
            for (int ky = 0; ky < m2.height; ky++)
                for (int kx = 0; kx < m2.width; kx++)
                    sum += m1.map[top + ky][left + kx] * m2.map[ky][kx];
            result->map[y][x] = sum;
        }
    }
}

void median_filter(Matrix m1, Matrix m2, Matrix *result) {
    int H = m1.height, W = m1.width;
    int kH = m2.height, kW = m2.width;
    int anchorY = (kH - 1) / 2, anchorX = (kW - 1) / 2;
    int kernelSize = kH * kW;
    double *window = (double *)malloc(kernelSize * sizeof(double));
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int top   = y - anchorY, bot   = y + (kH - 1 - anchorY);
            int left  = x - anchorX, right = x + (kW - 1 - anchorX);
            if (top < 0 || bot >= H || left < 0 || right >= W) {
                result->map[y][x] = 0.0;
                continue;
            }
            int idx = 0;
            for (int ky = 0; ky < kH; ky++)
                for (int kx = 0; kx < kW; kx++)
                    window[idx++] = m1.map[top + ky][left + kx];
            qsort(window, kernelSize, sizeof(double), cmpDouble);
            result->map[y][x] = window[kernelSize / 2];
        }
    }
    free(window);
}

void gaussian_filter(Matrix m1, Matrix *result) {
    // Unnormalised weights from main.c canny():
    double gaussData[5][5] = {
        { 2,  4,  5,  4,  2},
        { 4,  9, 12,  9,  4},
        { 5, 12, 15, 12,  5},
        { 4,  9, 12,  9,  4},
        { 2,  4,  5,  4,  2}
    };
    Matrix kGauss = createMatrix(5, 5);
    double gaussSum = 0.0;
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) gaussSum += gaussData[i][j];
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            kGauss.map[i][j] = gaussData[i][j] / gaussSum;

    _convolve(&m1, &kGauss, result);
    deleteMatrix(kGauss);
}

// Apply uniform averaging filter
Image applyAveragingFilter(Image img) {
    int H = img.height, W = img.width;

    Matrix kernel5 = createMatrix(5, 5);
    for (int ky = 0; ky < 5; ky++)
        for (int kx = 0; kx < 5; kx++)
            kernel5.map[ky][kx] = 1.0 / 25.0;

    Matrix imgMx = createMatrix(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            imgMx.map[y][x] = (double)img.map[y][x].i;

    Matrix result = createMatrix(H, W);
    smoothing_filter(imgMx, kernel5, &result);

    Image out = createImage(H, W);
    matrixToImageScaled(&result, &out);

    deleteMatrix(kernel5);
    deleteMatrix(imgMx);
    deleteMatrix(result);
    return out;
}

// Apply true Gaussian smoothing
Image applyGaussianFilter(Image img) {
    int H = img.height, W = img.width;

    Matrix imgMx = createMatrix(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            imgMx.map[y][x] = (double)img.map[y][x].i;

    Matrix result = createMatrix(H, W);
    gaussian_filter(imgMx, &result);

    Image out = createImage(H, W);
    matrixToImageScaled(&result, &out);

    deleteMatrix(imgMx);
    deleteMatrix(result);
    return out;
}

// Apply median filter
Image applyMedianFilter(Image img) {
    int H = img.height, W = img.width;

    Matrix kernel5 = createMatrix(5, 5);
    for (int ky = 0; ky < 5; ky++)
        for (int kx = 0; kx < 5; kx++)
            kernel5.map[ky][kx] = 1.0 / 25.0;

    Matrix imgMx = createMatrix(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            imgMx.map[y][x] = (double)img.map[y][x].i;

    Matrix result = createMatrix(H, W);
    median_filter(imgMx, kernel5, &result);

    Image out = createImage(H, W);
    matrixToImageScaled(&result, &out);

    deleteMatrix(kernel5);
    deleteMatrix(imgMx);
    deleteMatrix(result);
    return out;
}

int main(int argc, char *argv[]) {
    char *inputFile = (argc > 1) ? argv[1] : "input.ppm";

    printf("Reading '%s' ...\n", inputFile);
    Image img = readImage(inputFile);
    printf("Image: %d x %d\n", img.width, img.height);

#ifdef _WIN32
    mkdir("output_images");
#else
    mkdir("output_images", 0755);
#endif

    // 1. Uniform averaging filter
    Image avgOut = applyAveragingFilter(img);
    writeImage(avgOut, "output_images/averaging_smoothed.ppm");
    printf("Averaging filter  -> output_images/averaging_smoothed.ppm\n");

    // 2. True Gaussian filter
    Image gaussOut = applyGaussianFilter(img);
    writeImage(gaussOut, "output_images/gaussian_smoothed.ppm");
    printf("Gaussian filter   -> output_images/gaussian_smoothed.ppm\n");

    // 3. Median filter
    Image medOut = applyMedianFilter(img);
    writeImage(medOut, "output_images/median_filtered.ppm");
    printf("Median filter     -> output_images/median_filtered.ppm\n");

    deleteImage(img);
    deleteImage(avgOut);
    deleteImage(gaussOut);
    deleteImage(medOut);

    printf("Done. All outputs in output_images/\n");
    return 0;
}

// compile: gcc gaussian.c ../netpbm.c -o gaussian -lm
// run:     gaussian.exe ../input_images/coastline.ppm
