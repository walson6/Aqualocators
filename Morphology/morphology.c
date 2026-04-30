// morphology.c
// CS136 – Morphological Operations
//
// Contains:
//   function_imageBlackWhite() intensity threshold to binary image (0/255)
//   function_expandImage() morphological dilation (4-connectivity)
//   function_shrinkImage() morphological erosion (4-connectivity)
//   morphologicalOpen() erosion  then dilation (removes thin protrusions and noise)
//   morphologicalClose() dilation then erosion (fills small gaps and connects fragmented edges)
//
// Compile:
//   gcc morphology.c netpbm.c -o morphology -lm
// Run:
//   ./morphology <input.ppm|pgm> [threshold]
//   ./morphology coastline.ppm 128

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include "../netpbm.h"

Image function_imageBlackWhite(Image img, int threshold) {
    Image bw = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++) {
            unsigned char v = (img.map[y][x].i < threshold) ? 0 : 255;
            bw.map[y][x].i = bw.map[y][x].r =
            bw.map[y][x].g = bw.map[y][x].b = v;
        }
    return bw;
}

Image function_expandImage(Image img) {
    Image expanded = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            expanded.map[y][x] = img.map[y][x];

    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            if (img.map[y][x].i == 0) {
                int hit = (y > 0              && img.map[y-1][x].i == 255) ||
                          (y < img.height - 1 && img.map[y+1][x].i == 255) ||
                          (x > 0              && img.map[y][x-1].i == 255) ||
                          (x < img.width  - 1 && img.map[y][x+1].i == 255);
                if (hit)
                    expanded.map[y][x].i = expanded.map[y][x].r =
                    expanded.map[y][x].g = expanded.map[y][x].b = 255;
            }
    return expanded;
}

Image function_shrinkImage(Image img) {
    Image shrunk = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            shrunk.map[y][x] = img.map[y][x];

    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            if (img.map[y][x].i == 255) {
                int hit = (y > 0              && img.map[y-1][x].i == 0) ||
                          (y < img.height - 1 && img.map[y+1][x].i == 0) ||
                          (x > 0              && img.map[y][x-1].i == 0) ||
                          (x < img.width  - 1 && img.map[y][x+1].i == 0);
                if (hit)
                    shrunk.map[y][x].i = shrunk.map[y][x].r =
                    shrunk.map[y][x].g = shrunk.map[y][x].b = 0;
            }
    return shrunk;
}

Image morphologicalOpen(Image img, int passes) {
    Image current = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            current.map[y][x] = img.map[y][x];

    for (int p = 0; p < passes; p++) {
        Image eroded  = function_shrinkImage(current);
        Image dilated = function_expandImage(eroded);
        deleteImage(eroded);
        deleteImage(current);
        current = dilated;
    }
    return current;
}

Image morphologicalClose(Image img, int passes) {
    Image current = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            current.map[y][x] = img.map[y][x];

    for (int p = 0; p < passes; p++) {
        Image dilated = function_expandImage(current);
        Image eroded  = function_shrinkImage(dilated);
        deleteImage(dilated);
        deleteImage(current);
        current = eroded;
    }
    return current;
}

typedef struct { unsigned char r, g, b; } RGB;

static int findRoot(int *labels, int x) {
    while (labels[x] != x) x = labels[x];
    return x;
}

static void unionLabels(int *labels, int a, int b) {
    int ra = findRoot(labels, a), rb = findRoot(labels, b);
    if (ra != rb) labels[rb] = ra;
}

int labelAndColorImage(Image img, int sizeThreshold, Image *outputImage) {
    int h = img.height, w = img.width;
    Matrix labels = createMatrix(h, w);
    int nextLabel = 1;
    int *equiv    = (int *)malloc(sizeof(int) * h * w);
    for (int i = 0; i < h * w; i++) equiv[i] = i;

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            if (img.map[y][x].i != 255) continue;
            int L = (x > 0) ? (int)labels.map[y][x-1] : 0;
            int U = (y > 0) ? (int)labels.map[y-1][x]  : 0;
            if      (L == 0 && U == 0) { labels.map[y][x] = nextLabel++; }
            else if (L != 0 && U == 0) { labels.map[y][x] = L; }
            else if (L == 0 && U != 0) { labels.map[y][x] = U; }
            else { labels.map[y][x] = L; if (L != U) unionLabels(equiv, L, U); }
        }
    }

    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
            if ((int)labels.map[y][x] != 0)
                labels.map[y][x] = findRoot(equiv, (int)labels.map[y][x]);

    int *labelSizes = (int *)calloc(nextLabel, sizeof(int));
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
            if ((int)labels.map[y][x] != 0)
                labelSizes[(int)labels.map[y][x]]++;

    srand((unsigned int)time(NULL));
    RGB *labelColor = (RGB *)malloc(nextLabel * sizeof(RGB));
    for (int i = 0; i < nextLabel; i++) {
        if (labelSizes[i] >= sizeThreshold) {
            labelColor[i].r = 50 + rand() % 206;
            labelColor[i].g = 50 + rand() % 206;
            labelColor[i].b = 50 + rand() % 206;
        } else {
            labelColor[i].r = labelColor[i].g = labelColor[i].b = 0;
        }
    }

    *outputImage = createImage(h, w);
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
            (*outputImage).map[y][x] = img.map[y][x];

    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            int lbl = (int)labels.map[y][x];
            if (lbl != 0 && labelSizes[lbl] >= sizeThreshold) {
                (*outputImage).map[y][x].r = labelColor[lbl].r;
                (*outputImage).map[y][x].g = labelColor[lbl].g;
                (*outputImage).map[y][x].b = labelColor[lbl].b;
            }
        }

    int count = 0;
    for (int i = 1; i < nextLabel; i++)
        if (labelSizes[i] >= sizeThreshold) count++;

    deleteMatrix(labels);
    free(equiv); free(labelSizes); free(labelColor);
    return count;
}

int main(int argc, char *argv[]) {
    char *inputFile  = (argc > 1) ? argv[1] : "input.ppm";
    int   threshold  = (argc > 2) ? atoi(argv[2]) : 128;
    int   sizeThresh = (argc > 3) ? atoi(argv[3]) : 50;

    printf("Reading '%s' ...\n", inputFile);
    Image img = readImage(inputFile);
    printf("Image: %d x %d  bw_threshold=%d  size_threshold=%d\n",
           img.width, img.height, threshold, sizeThresh);

#ifdef _WIN32
    mkdir("output_images");
#else
    mkdir("output_images", 0755);
#endif

    // 1. Binary image
    Image bw = function_imageBlackWhite(img, threshold);
    writeImage(bw, "output_images/binary.pbm");
    printf("Binary image       -> output_images/binary.pbm\n");

    // 2. Dilation
    Image dilated = function_expandImage(bw);
    writeImage(dilated, "output_images/dilated.pbm");
    printf("Dilation           -> output_images/dilated.pbm\n");

    // 3. Erosion
    Image eroded = function_shrinkImage(bw);
    writeImage(eroded, "output_images/eroded.pbm");
    printf("Erosion            -> output_images/eroded.pbm\n");

    // 4. Opening (erode then dilate)
    Image opened = morphologicalOpen(bw, 1);
    writeImage(opened, "output_images/opened.pbm");
    printf("Morphological open -> output_images/opened.pbm\n");

    // 5. Closing (dilate then erode)
    Image closed = morphologicalClose(bw, 1);
    writeImage(closed, "output_images/closed.pbm");
    printf("Morphological close-> output_images/closed.pbm\n");

    // 6. Connected components
    Image colored;
    int numComponents = labelAndColorImage(closed, sizeThresh, &colored);
    writeImage(colored, "output_images/components.ppm");
    printf("Connected components (>=%d px): %d  -> output_images/components.ppm\n",
           sizeThresh, numComponents);

    deleteImage(img);
    deleteImage(bw);
    deleteImage(dilated);
    deleteImage(eroded);
    deleteImage(opened);
    deleteImage(closed);
    deleteImage(colored);

    printf("Done.\n");
    return 0;
}

// compile: gcc morphology.c netpbm.c -o morphology -lm
// run:     ./morphology coastline.ppm 128 50
