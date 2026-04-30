// canny.c
// CS136 – Canny Edge Detection
//
// Pipeline:
//   Stage 1  – 5×5 Gaussian smoothing (sigma ≈ 1.4)
//   Stage 2  – Sobel gradient (magnitude + direction)
//   Stage 3  – Non-maximum suppression (8 gradient directions → 4 bins)
//   Stage 4  – Double threshold (high = 15 % of max, low = 50 % of high)
//   Stage 5  – Hysteresis edge linking (forward + backward passes)
//
// The high/low threshold ratio can be overridden via command-line:
//   ./canny <input> [highRatio] [lowRatio]
//   Defaults: highRatio=0.15  lowRatio=0.50
//
// Compile:
//   gcc canny.c ../netpbm.c -o canny -lm
// Run:
//   canny.exe ../input_images/coastline.ppm
//   canny.exe ../coastline.ppm 0.20 0.40

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "../netpbm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

Image canny(Image img, double highRatio, double lowRatio) {
    int H = img.height, W = img.width;

    // Stage 1: Gaussian smoothing
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

    Matrix imgMx = createMatrix(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            imgMx.map[y][x] = (double)img.map[y][x].i;

    Matrix smoothed = createMatrix(H, W);
    _convolve(&imgMx, &kGauss, &smoothed);

    // Stage 2: Sobel gradients
    double gxData[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    double gyData[3][3] = {{-1,-2,-1}, { 0, 0, 0}, { 1, 2, 1}};
    Matrix kGx = createMatrix(3, 3);
    Matrix kGy = createMatrix(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            kGx.map[i][j] = gxData[i][j];
            kGy.map[i][j] = gyData[i][j];
        }

    Matrix Gx  = createMatrix(H, W);
    Matrix Gy  = createMatrix(H, W);
    Matrix mag = createMatrix(H, W);
    Matrix dir = createMatrix(H, W);

    _convolve(&smoothed, &kGx, &Gx);
    _convolve(&smoothed, &kGy, &Gy);

    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            mag.map[y][x] = sqrt(Gx.map[y][x] * Gx.map[y][x] +
                                 Gy.map[y][x] * Gy.map[y][x]);
            dir.map[y][x] = atan2(Gy.map[y][x], Gx.map[y][x]);
        }

    // Stage 3: Non-maximum suppression
    Matrix nms = createMatrix(H, W);
    for (int y = 1; y < H - 1; y++) {
        for (int x = 1; x < W - 1; x++) {
            double angle = dir.map[y][x] * 180.0 / M_PI;
            if (angle < 0) angle += 180.0;
            double n1, n2;
            if      (angle < 22.5  || angle >= 157.5) { n1 = mag.map[y][x-1];   n2 = mag.map[y][x+1]; }
            else if (angle < 67.5)                     { n1 = mag.map[y-1][x+1]; n2 = mag.map[y+1][x-1]; }
            else if (angle < 112.5)                    { n1 = mag.map[y-1][x];   n2 = mag.map[y+1][x]; }
            else                                       { n1 = mag.map[y-1][x-1]; n2 = mag.map[y+1][x+1]; }
            nms.map[y][x] = (mag.map[y][x] >= n1 && mag.map[y][x] >= n2)
                            ? mag.map[y][x] : 0.0;
        }
    }

    // Stage 4: Double threshold
    double maxMag = 0.0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            if (nms.map[y][x] > maxMag) maxMag = nms.map[y][x];

    double highThresh = maxMag * highRatio;
    double lowThresh  = highThresh * lowRatio;

    printf("  NMS max=%.2f  highThresh=%.2f  lowThresh=%.2f\n",
           maxMag, highThresh, lowThresh);

    Matrix edges = createMatrix(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            if      (nms.map[y][x] >= highThresh) edges.map[y][x] = 255.0;
            else if (nms.map[y][x] >= lowThresh)  edges.map[y][x] = 128.0;
            else                                   edges.map[y][x] = 0.0;
        }

    // Stage 5: Hysteresis thresholding
    int changed = 1;
    while (changed) {
        changed = 0;
        // Forward pass
        for (int y = 1; y < H - 1; y++)
            for (int x = 1; x < W - 1; x++)
                if (edges.map[y][x] == 128.0) {
                    int hasStrong = 0;
                    for (int dy = -1; dy <= 1 && !hasStrong; dy++)
                        for (int dx = -1; dx <= 1 && !hasStrong; dx++)
                            if (edges.map[y+dy][x+dx] == 255.0) hasStrong = 1;
                    if (hasStrong) { edges.map[y][x] = 255.0; changed = 1; }
                }
        // Backward pass
        for (int y = H - 2; y >= 1; y--)
            for (int x = W - 2; x >= 1; x--)
                if (edges.map[y][x] == 128.0) {
                    int hasStrong = 0;
                    for (int dy = -1; dy <= 1 && !hasStrong; dy++)
                        for (int dx = -1; dx <= 1 && !hasStrong; dx++)
                            if (edges.map[y+dy][x+dx] == 255.0) hasStrong = 1;
                    if (hasStrong) { edges.map[y][x] = 255.0; changed = 1; }
                }
    }
    // Suppress remaining weak edges
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            if (edges.map[y][x] == 128.0) edges.map[y][x] = 0.0;

    // Build output image
    Image out = createImage(H, W);
    matrixToImageScaled(&edges, &out);

    deleteMatrix(kGauss); deleteMatrix(imgMx);  deleteMatrix(smoothed);
    deleteMatrix(kGx);    deleteMatrix(kGy);
    deleteMatrix(Gx);     deleteMatrix(Gy);
    deleteMatrix(mag);    deleteMatrix(dir);
    deleteMatrix(nms);    deleteMatrix(edges);
    return out;
}

int main(int argc, char *argv[]) {
    char   *inputFile = (argc > 1) ? argv[1] : "../input_images/input.ppm";
    double  highRatio = (argc > 2) ? atof(argv[2]) : 0.15;
    double  lowRatio  = (argc > 3) ? atof(argv[3]) : 0.50;

    printf("Reading '%s' ...\n", inputFile);
    Image img = readImage(inputFile);
    printf("Image: %d x %d\n", img.width, img.height);
    printf("Thresholds: high=%.2f * max,  low=%.2f * high\n", highRatio, lowRatio);

#ifdef _WIN32
    mkdir("output_images");
#else
    mkdir("output_images", 0755);
#endif

    Image edgeImg = canny(img, highRatio, lowRatio);
    writeImage(edgeImg, "output_images/canny_edges.ppm");
    printf("Canny edges -> output_images/canny_edges.ppm\n");

    writeImage(edgeImg, "output_images/canny_edges.pgm");
    printf("Canny edges -> output_images/canny_edges.pgm\n");

    // Count detected edge pixels if you want to know the percentage of edge pixels
    long edgeCount = 0;
    for (int y = 0; y < edgeImg.height; y++)
        for (int x = 0; x < edgeImg.width; x++)
            if (edgeImg.map[y][x].i > 128) edgeCount++;
    printf("Edge pixel count: %ld  (%.2f%% of image)\n",
           edgeCount,
           100.0 * edgeCount / (edgeImg.height * edgeImg.width));

    deleteImage(img);
    deleteImage(edgeImg);

    return 0;
}

// compile: gcc canny.c ../netpbm.c -o canny -lm
// run:     canny.exe ../input_images/coastline.ppm
// run:     canny.exe ../input_images/coastline.ppm 0.20 0.40   (custom thresholds)
