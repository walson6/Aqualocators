// creative_morphology.c
// CS136 - Creative Exploration: Morphology-Assisted Segmentation
//
// This module adapts the Morphology assignment code and combines it with
// automatic thresholding and Sobel edge detection. The goal is to use
// morphological operations as a creative segmentation tool instead of only
// showing dilation/erosion independently.
//
// Pipeline:
//   1. Read PGM/PPM image using the shared NetPBM loader.
//   2. Build an Otsu binary mask from grayscale intensity.
//   3. Build a Sobel edge mask using a threshold ratio.
//   4. Combine the threshold mask with edge support.
//   5. Apply closing to connect gaps, then opening to remove small noise.
//   6. Color connected components and blend them over the original image.
//
// Compile:
//   gcc creative_morphology.c ../netpbm.c -o creative_morphology -lm
//
// Run:
//   ./creative_morphology ../input_images/coastline.ppm 0.25 2 1 800
//   ./creative_morphology ../input_images/tahoe.ppm 0.20 3 1 1200
//
// Arguments:
//   argv[1] input image path                   default: ../input_images/coastline.ppm
//   argv[2] Sobel edge threshold ratio         default: 0.25
//   argv[3] closing passes                     default: 2
//   argv[4] opening passes                     default: 1
//   argv[5] minimum component size in pixels   default: 800

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "../netpbm.h"

typedef struct {
    unsigned char r, g, b;
} RGB;

static void getOutputPrefix(char *inputPath, char *prefix, int prefixSize) {
    char temp[256];
    char *base = strrchr(inputPath, '/');
    base = (base == NULL) ? inputPath : base + 1;

    snprintf(temp, sizeof(temp), "%s", base);
    char *dot = strrchr(temp, '.');
    if (dot != NULL) *dot = '\0';

    snprintf(prefix, prefixSize, "output_images/%s", temp);
}

static Image cloneImage(Image img) {
    Image out = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            out.map[y][x] = img.map[y][x];
    return out;
}

static void setGrayPixel(Image img, int y, int x, unsigned char v) {
    img.map[y][x].i = img.map[y][x].r =
    img.map[y][x].g = img.map[y][x].b = v;
}

static int otsuThreshold(Image img) {
    long hist[256] = {0};
    int total = img.height * img.width;

    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            hist[img.map[y][x].i]++;

    double sum = 0.0;
    for (int t = 0; t < 256; t++) sum += t * hist[t];

    double sumB = 0.0;
    long wB = 0;
    long wF = 0;
    double maxBetween = -1.0;
    int threshold = 128;

    for (int t = 0; t < 256; t++) {
        wB += hist[t];
        if (wB == 0) continue;

        wF = total - wB;
        if (wF == 0) break;

        sumB += t * hist[t];
        double mB = sumB / wB;
        double mF = (sum - sumB) / wF;
        double between = (double)wB * (double)wF * (mB - mF) * (mB - mF);

        if (between > maxBetween) {
            maxBetween = between;
            threshold = t;
        }
    }

    return threshold;
}

static Image thresholdImage(Image img, int threshold) {
    Image out = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++) {
            unsigned char v = (img.map[y][x].i >= threshold) ? 255 : 0;
            setGrayPixel(out, y, x, v);
        }
    return out;
}

static Image sobelEdgeMask(Image img, double thresholdRatio) {
    int H = img.height, W = img.width;
    Matrix mag = createMatrix(H, W);
    double maxMag = 0.0;

    for (int y = 1; y < H - 1; y++) {
        for (int x = 1; x < W - 1; x++) {
            double gx =
                -img.map[y-1][x-1].i + img.map[y-1][x+1].i
                -2.0 * img.map[y][x-1].i + 2.0 * img.map[y][x+1].i
                -img.map[y+1][x-1].i + img.map[y+1][x+1].i;
            double gy =
                -img.map[y-1][x-1].i - 2.0 * img.map[y-1][x].i - img.map[y-1][x+1].i
                +img.map[y+1][x-1].i + 2.0 * img.map[y+1][x].i + img.map[y+1][x+1].i;
            mag.map[y][x] = sqrt(gx * gx + gy * gy);
            if (mag.map[y][x] > maxMag) maxMag = mag.map[y][x];
        }
    }

    double threshold = maxMag * thresholdRatio;
    Image edges = createImage(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            unsigned char v = (mag.map[y][x] >= threshold) ? 255 : 0;
            setGrayPixel(edges, y, x, v);
        }

    deleteMatrix(mag);
    return edges;
}

static Image dilate(Image img) {
    Image out = cloneImage(img);
    int H = img.height, W = img.width;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            if (img.map[y][x].i == 255) continue;
            int hit = 0;
            for (int dy = -1; dy <= 1 && !hit; dy++)
                for (int dx = -1; dx <= 1 && !hit; dx++) {
                    int yy = y + dy, xx = x + dx;
                    if (yy >= 0 && yy < H && xx >= 0 && xx < W && img.map[yy][xx].i == 255)
                        hit = 1;
                }
            if (hit) setGrayPixel(out, y, x, 255);
        }
    }

    return out;
}

static Image erode(Image img) {
    Image out = cloneImage(img);
    int H = img.height, W = img.width;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            if (img.map[y][x].i == 0) continue;
            int touchesBackground = 0;
            for (int dy = -1; dy <= 1 && !touchesBackground; dy++)
                for (int dx = -1; dx <= 1 && !touchesBackground; dx++) {
                    int yy = y + dy, xx = x + dx;
                    if (yy < 0 || yy >= H || xx < 0 || xx >= W || img.map[yy][xx].i == 0)
                        touchesBackground = 1;
                }
            if (touchesBackground) setGrayPixel(out, y, x, 0);
        }
    }

    return out;
}

static Image morphClose(Image img, int passes) {
    Image current = cloneImage(img);
    for (int p = 0; p < passes; p++) {
        Image expanded = dilate(current);
        Image shrunk = erode(expanded);
        deleteImage(current);
        deleteImage(expanded);
        current = shrunk;
    }
    return current;
}

static Image morphOpen(Image img, int passes) {
    Image current = cloneImage(img);
    for (int p = 0; p < passes; p++) {
        Image shrunk = erode(current);
        Image expanded = dilate(shrunk);
        deleteImage(current);
        deleteImage(shrunk);
        current = expanded;
    }
    return current;
}

static Image combineThresholdAndEdges(Image thresholdMask, Image edgeMask) {
    Image combined = cloneImage(thresholdMask);

    for (int y = 0; y < thresholdMask.height; y++) {
        for (int x = 0; x < thresholdMask.width; x++) {
            if (thresholdMask.map[y][x].i == 255 || edgeMask.map[y][x].i == 255)
                setGrayPixel(combined, y, x, 255);
            else
                setGrayPixel(combined, y, x, 0);
        }
    }

    return combined;
}

static int findRoot(int *labels, int x) {
    while (labels[x] != x) {
        labels[x] = labels[labels[x]];
        x = labels[x];
    }
    return x;
}

static void unionLabels(int *labels, int a, int b) {
    int ra = findRoot(labels, a);
    int rb = findRoot(labels, b);
    if (ra != rb) labels[rb] = ra;
}

static int labelAndColor(Image binary, int sizeThreshold, Image *colored) {
    int H = binary.height, W = binary.width;
    Matrix labels = createMatrix(H, W);
    int maxLabels = H * W + 1;
    int *equiv = (int *)malloc((size_t)maxLabels * sizeof(int));
    int nextLabel = 1;

    for (int i = 0; i < maxLabels; i++) equiv[i] = i;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            if (binary.map[y][x].i != 255) continue;

            int L = (x > 0) ? (int)labels.map[y][x - 1] : 0;
            int U = (y > 0) ? (int)labels.map[y - 1][x] : 0;
            int UL = (y > 0 && x > 0) ? (int)labels.map[y - 1][x - 1] : 0;
            int UR = (y > 0 && x < W - 1) ? (int)labels.map[y - 1][x + 1] : 0;
            int chosen = L ? L : (U ? U : (UL ? UL : UR));

            if (chosen == 0) {
                labels.map[y][x] = nextLabel++;
            } else {
                labels.map[y][x] = chosen;
                if (L && L != chosen) unionLabels(equiv, chosen, L);
                if (U && U != chosen) unionLabels(equiv, chosen, U);
                if (UL && UL != chosen) unionLabels(equiv, chosen, UL);
                if (UR && UR != chosen) unionLabels(equiv, chosen, UR);
            }
        }
    }

    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            if ((int)labels.map[y][x] != 0)
                labels.map[y][x] = findRoot(equiv, (int)labels.map[y][x]);

    int *sizes = (int *)calloc((size_t)nextLabel, sizeof(int));
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            if ((int)labels.map[y][x] != 0)
                sizes[(int)labels.map[y][x]]++;

    RGB palette[] = {
        {230,  57,  70}, { 42, 157, 143}, { 69, 123, 157}, {244, 162,  97},
        {131,  56, 236}, {255, 190,  11}, {  6, 214, 160}, { 17, 138, 178},
        {239,  71, 111}, {118, 200, 147}, { 29,  53,  87}, {233, 196, 106}
    };
    int paletteSize = (int)(sizeof(palette) / sizeof(palette[0]));
    RGB *colors = (RGB *)calloc((size_t)nextLabel, sizeof(RGB));
    int componentCount = 0;

    for (int i = 1; i < nextLabel; i++) {
        if (sizes[i] >= sizeThreshold) {
            colors[i] = palette[componentCount % paletteSize];
            componentCount++;
        }
    }

    *colored = createImage(H, W);
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int lbl = (int)labels.map[y][x];
            if (lbl != 0 && sizes[lbl] >= sizeThreshold) {
                (*colored).map[y][x].r = colors[lbl].r;
                (*colored).map[y][x].g = colors[lbl].g;
                (*colored).map[y][x].b = colors[lbl].b;
                (*colored).map[y][x].i = (unsigned char)((colors[lbl].r + colors[lbl].g + colors[lbl].b) / 3);
            } else {
                setGrayPixel(*colored, y, x, 0);
            }
        }
    }

    deleteMatrix(labels);
    free(equiv);
    free(sizes);
    free(colors);
    return componentCount;
}

static Image makeOverlay(Image original, Image components) {
    Image out = createImage(original.height, original.width);

    for (int y = 0; y < original.height; y++) {
        for (int x = 0; x < original.width; x++) {
            if (components.map[y][x].r || components.map[y][x].g || components.map[y][x].b) {
                out.map[y][x].r = (unsigned char)(0.55 * original.map[y][x].r + 0.45 * components.map[y][x].r);
                out.map[y][x].g = (unsigned char)(0.55 * original.map[y][x].g + 0.45 * components.map[y][x].g);
                out.map[y][x].b = (unsigned char)(0.55 * original.map[y][x].b + 0.45 * components.map[y][x].b);
            } else {
                out.map[y][x] = original.map[y][x];
            }
            out.map[y][x].i = (unsigned char)((out.map[y][x].r + out.map[y][x].g + out.map[y][x].b) / 3);
        }
    }

    return out;
}

int main(int argc, char *argv[]) {
    char *inputFile = (argc > 1) ? argv[1] : "../input_images/coastline.ppm";
    double edgeRatio = (argc > 2) ? atof(argv[2]) : 0.25;
    int closePasses = (argc > 3) ? atoi(argv[3]) : 2;
    int openPasses = (argc > 4) ? atoi(argv[4]) : 1;
    int sizeThreshold = (argc > 5) ? atoi(argv[5]) : 800;

    if (edgeRatio < 0.01) edgeRatio = 0.01;
    if (edgeRatio > 1.00) edgeRatio = 1.00;
    if (closePasses < 0) closePasses = 0;
    if (openPasses < 0) openPasses = 0;
    if (sizeThreshold < 1) sizeThreshold = 1;

    printf("Reading '%s' ...\n", inputFile);
    Image img = readImage(inputFile);
    printf("Image: %d x %d  edgeRatio=%.2f  close=%d  open=%d  minSize=%d\n",
           img.width, img.height, edgeRatio, closePasses, openPasses, sizeThreshold);

#ifdef _WIN32
    mkdir("output_images");
#else
    mkdir("output_images", 0755);
#endif

    int threshold = otsuThreshold(img);
    printf("Otsu threshold: %d\n", threshold);

    Image otsu = thresholdImage(img, threshold);
    Image edges = sobelEdgeMask(img, edgeRatio);
    Image combined = combineThresholdAndEdges(otsu, edges);
    Image closed = morphClose(combined, closePasses);
    Image cleaned = morphOpen(closed, openPasses);

    Image components;
    int count = labelAndColor(cleaned, sizeThreshold, &components);
    Image overlay = makeOverlay(img, components);

    char prefix[320];
    char otsuFile[384], edgeFile[384], combinedFile[384], cleanedFile[384];
    char componentFile[384], overlayFile[384];
    getOutputPrefix(inputFile, prefix, sizeof(prefix));
    snprintf(otsuFile, sizeof(otsuFile), "%s_otsu_mask.pbm", prefix);
    snprintf(edgeFile, sizeof(edgeFile), "%s_sobel_edge_mask.pbm", prefix);
    snprintf(combinedFile, sizeof(combinedFile), "%s_combined_mask.pbm", prefix);
    snprintf(cleanedFile, sizeof(cleanedFile), "%s_morph_cleaned.pbm", prefix);
    snprintf(componentFile, sizeof(componentFile), "%s_morph_components.ppm", prefix);
    snprintf(overlayFile, sizeof(overlayFile), "%s_morph_overlay.ppm", prefix);

    writeImage(otsu, otsuFile);
    writeImage(edges, edgeFile);
    writeImage(combined, combinedFile);
    writeImage(cleaned, cleanedFile);
    writeImage(components, componentFile);
    writeImage(overlay, overlayFile);

    printf("Otsu mask          -> %s\n", otsuFile);
    printf("Sobel edge mask    -> %s\n", edgeFile);
    printf("Combined mask      -> %s\n", combinedFile);
    printf("Cleaned mask       -> %s\n", cleanedFile);
    printf("Components (%d)    -> %s\n", count, componentFile);
    printf("Overlay            -> %s\n", overlayFile);

    deleteImage(img);
    deleteImage(otsu);
    deleteImage(edges);
    deleteImage(combined);
    deleteImage(closed);
    deleteImage(cleaned);
    deleteImage(components);
    deleteImage(overlay);

    printf("Done.\n");
    return 0;
}
