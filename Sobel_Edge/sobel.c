// sobel.c
// CS136 – Sobel Edge Detection
//
// Contains:
//   sobel() gradient magnitude edge detector
//
// The algorithm:
//   1. Load intensity channel into a Matrix
//   2. Convolve with 3×3 Gx and Gy kernels
//   3. Compute magnitude = sqrt(Gx²+Gy²) at each pixel
//   4. Scale magnitude into [0,255] and write output image
//
// Compile:
//   gcc sobel.c ../netpbm.c -o sobel -lm
// Run:
//   ./sobel <input.ppm|pgm>

#include <stdio.h>
#include <stdlib.h>
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

Image sobel(Image img) {
    int H = img.height, W = img.width;

    // 3×3 Sobel kernels
    double gxData[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    double gyData[3][3] = {{-1,-2,-1}, { 0, 0, 0}, { 1, 2, 1}};

    Matrix kGx = createMatrix(3, 3);
    Matrix kGy = createMatrix(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            kGx.map[i][j] = gxData[i][j];
            kGy.map[i][j] = gyData[i][j];
        }

    Matrix imgMx = createMatrix(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            imgMx.map[y][x] = (double)img.map[y][x].i;

    Matrix Gx  = createMatrix(H, W);
    Matrix Gy  = createMatrix(H, W);
    Matrix mag = createMatrix(H, W);

    _convolve(&imgMx, &kGx, &Gx);
    _convolve(&imgMx, &kGy, &Gy);

    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            mag.map[y][x] = sqrt(Gx.map[y][x] * Gx.map[y][x] +
                                 Gy.map[y][x] * Gy.map[y][x]);

    Image out = createImage(H, W);
    matrixToImageScaled(&mag, &out);

    deleteMatrix(kGx); deleteMatrix(kGy); deleteMatrix(imgMx);
    deleteMatrix(Gx);  deleteMatrix(Gy);  deleteMatrix(mag);
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

    Image edgeImg = sobel(img);
    writeImage(edgeImg, "output_images/sobel_edges.ppm");
    printf("Sobel edges -> output_images/sobel_edges.ppm\n");

    // Also save a grayscale (PGM) version for quantitative comparison
    writeImage(edgeImg, "output_images/sobel_edges.pgm");
    printf("Sobel edges -> output_images/sobel_edges.pgm\n");

    deleteImage(img);
    deleteImage(edgeImg);

    printf("Done.\n");
    return 0;
}

// compile: gcc sobel.c ../netpbm.c -o sobel -lm
// run:     sobel.exe ../coastline.ppm
