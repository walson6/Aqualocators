// hough.c
// CS136 – Hough Transform Module
// Coastline & Lake Detection via Line and Circle Hough Transforms
//
// Contains:
//   buildEdgeMatrix() Gaussian smooth → Sobel magnitude → threshold
//   houghTransformLines() accumulates votes in (alpha, d) parameter space.
//   isLocalMaximum() 8-neighbour local-max check.
//   findHoughMaxima() returns top-N well-separated peaks.
//   drawDetectedLines() renders detected lines onto an image.
//   houghTransformCircles() 3D accumulator indexed (cy, cx, radiusIndex).
//   projectHough3D() collapses radius axis (max-vote projection).
//   findCircleMaxima() 3-D local-max peak picking.
//   drawDetectedCircles() renders detected circles onto an image.
//   runLinePipeline() full line detection pipeline for one image.
//   runCirclePipeline() full circle detection pipeline for one image.
//
// Compile:
//   gcc hough.c ../netpbm.c -o hough -lm
//
// Run:
//   hough.exe ../input_images/coastline.ppm [mode] [numPeaks] [edgeThreshRatio]
//      mode: "lines" (default) or "circles"
//      numPeaks: how many lines/circles to detect  (default: 5)
//      edgeThreshRatio: edge threshold as fraction of max gradient (default: 0.15)
//
//   Example: detect 8 coastline lines:
//     hough.exe ../input_images/coastline.ppm lines 8 0.15
//
//   Example: detect lake boundary circles:
//     hough.exe ../input_images/tahoe.ppm circles 3 0.12

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "../netpbm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Hough parameter-space resolution
#define LINE_MAP_H  360   // number of angle bins  (alpha rows)
#define LINE_MAP_W  480   // number of distance bins (d columns)

// Circle radius search range (pixels)
#define CIRCLE_R_MIN   15
#define CIRCLE_R_MAX  300
#define CIRCLE_R_STEP   3

// Minimum separation between accepted peaks (parameter-space pixels for lines, spatial pixels for circles)
#define LINE_MIN_SEP   20.0
#define CIRCLE_MIN_SEP 40.0

typedef struct { int H, W, R; int *data; } Hough3D;

static Hough3D createHough3D(int H, int W, int R) {
    Hough3D h; h.H = H; h.W = W; h.R = R;
    h.data = (int *)calloc((size_t)H * W * R, sizeof(int));
    return h;
}
static void deleteHough3D(Hough3D h) { free(h.data); }
static inline int h3Idx(Hough3D *h, int y, int x, int ri) {
    return ri * (h->H * h->W) + y * h->W + x;
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
                result->map[y][x] = 0.0; continue;
            }
            double sum = 0.0;
            for (int ky = 0; ky < m2->height; ky++)
                for (int kx = 0; kx < m2->width; kx++)
                    sum += m1->map[top + ky][left + kx] * m2->map[ky][kx];
            result->map[y][x] = sum;
        }
    }
}

static void matrixToImageScaled(Matrix *m, Image *img) {
    int h = m->height, w = m->width;
    double mn = m->map[0][0], mx = m->map[0][0];
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            if (m->map[y][x] < mn) mn = m->map[y][x];
            if (m->map[y][x] > mx) mx = m->map[y][x];
        }
    double range = (mx - mn > 1e-10) ? (mx - mn) : 1.0;
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            unsigned char v = (unsigned char)(((m->map[y][x] - mn) / range) * 255.0 + 0.5);
            img->map[y][x].i = img->map[y][x].r =
            img->map[y][x].g = img->map[y][x].b = v;
        }
}

// Preprocessing: Build binary edge matrix
// Pipeline: Gaussian 3×3 smooth -> Sobel magnitude -> threshold at edgeThreshRatio × maxMagnitude -> binary {0, 255}
// We use a 3×3 Gaussian (faster) because the Hough accumulator already integrates over all edge pixels; fine detail matters less than coverage
Matrix buildEdgeMatrix(Image img, double edgeThreshRatio) {
    int H = img.height, W = img.width;

    // Load intensity
    Matrix imgMx = createMatrix(H, W);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            imgMx.map[y][x] = (double)img.map[y][x].i;

    // 3×3 Gaussian smoothing (sigma approximately 1.0)
    double gk[3][3] = {{1,2,1},{2,4,2},{1,2,1}};
    Matrix kG = createMatrix(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) kG.map[i][j] = gk[i][j] / 16.0;

    Matrix smoothed = createMatrix(H, W);
    _convolve(&imgMx, &kG, &smoothed);
    deleteMatrix(kG); deleteMatrix(imgMx);

    // Sobel Gx and Gy matrices
    double gxd[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
    double gyd[3][3] = {{-1,-2,-1},{0,0,0},{1,2,1}};
    Matrix kGx = createMatrix(3,3), kGy = createMatrix(3,3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            kGx.map[i][j] = gxd[i][j];
            kGy.map[i][j] = gyd[i][j];
        }
    Matrix Gx = createMatrix(H,W), Gy = createMatrix(H,W);
    _convolve(&smoothed, &kGx, &Gx);
    _convolve(&smoothed, &kGy, &Gy);
    deleteMatrix(kGx); deleteMatrix(kGy); deleteMatrix(smoothed);

    // Magnitude + threshold
    Matrix edges = createMatrix(H, W);
    double maxMag = 0.0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            double mag = sqrt(Gx.map[y][x]*Gx.map[y][x] + Gy.map[y][x]*Gy.map[y][x]);
            edges.map[y][x] = mag;
            if (mag > maxMag) maxMag = mag;
        }
    double thresh = maxMag * edgeThreshRatio;
    // Exclude a 10-pixel border strip to prevent image-frame edges from flooding the Hough accumulator with spurious near-d=0 votes
    int edgePx = 0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            int onBorder = (y < 10 || y >= H-10 || x < 10 || x >= W-10);
            edges.map[y][x] = (!onBorder && edges.map[y][x] >= thresh) ? 255.0 : 0.0;
            if (edges.map[y][x] > 0) edgePx++;
        }
    printf("  Edge pixels: %d  (%.2f%% of image)  threshold=%.1f\n",
           edgePx, 100.0*edgePx/(H*W), thresh);

    deleteMatrix(Gx); deleteMatrix(Gy);
    return edges;
}

Matrix houghTransformLines(Matrix mxSpatial, int mapHeight, int mapWidth) {
    int m, n, angle, dist;
    double alpha;
    double maxD = sqrt((double)(mxSpatial.height * mxSpatial.height
                               + mxSpatial.width  * mxSpatial.width));
    Matrix mxParam = createMatrix(mapHeight, mapWidth);

    Matrix sincos = createMatrix(mapHeight, 2);
    for (angle = 0; angle < mapHeight; angle++) {
        alpha = -0.5*M_PI + M_PI*(double)angle/(double)mapHeight;
        sincos.map[angle][0] = sin(alpha);
        sincos.map[angle][1] = cos(alpha);
    }

    // Accumulate votes
    for (m = 0; m < mxSpatial.height; m++) {
        for (n = 0; n < mxSpatial.width; n++) {
            if (mxSpatial.map[m][n] < 128.0) continue;
            for (angle = 0; angle < mapHeight; angle++) {
                // Normal form: d = row*sin(alpha) + col*cos(alpha)
                double d = m * sincos.map[angle][0] + n * sincos.map[angle][1];
                dist = (int)(d / maxD * (mapWidth - 1) + 0.5);
                if (dist >= 0 && dist < mapWidth)
                    mxParam.map[angle][dist]++;
            }
        }
    }
    deleteMatrix(sincos);
    return mxParam;
}

// isLocalMaximum returns 1 if mx[m][n] is not exceeded by any of its (up to 8) neighbouring pixels
int isLocalMaximum(Matrix mx, int m, int n) {
    double strength = mx.map[m][n];
    int iMin = (m == 0) ? 0 : m-1;
    int iMax = (m == mx.height-1) ? m : m+1;
    int jMin = (n == 0) ? 0 : n-1;
    int jMax = (n == mx.width-1)  ? n : n+1;
    for (int i = iMin; i <= iMax; i++)
        for (int j = jMin; j <= jMax; j++)
            if (mx.map[i][j] > strength) return 0;
    return 1;
}

// insertMaxEntry maintains a sorted list of the top-N peaks
static void insertMaxEntry(Matrix mx, int vPos, int hPos, double strength) {
    int m, n = mx.width - 1;
    while (n > 0 && mx.map[2][n-1] < strength) {
        for (m = 0; m < 3; m++) mx.map[m][n] = mx.map[m][n-1];
        n--;
    }
    mx.map[0][n] = (double)vPos;
    mx.map[1][n] = (double)hPos;
    mx.map[2][n] = strength;
}

// findHoughMaxima returns the <number> highest well-separated local maxima in the Hough accumulator
Matrix findHoughMaxima(Matrix mx, int number, double minSeparation) {
    double minSepSq = minSeparation * minSeparation;
    Matrix maxima = createMatrix(3, number);
    for (int i = 0; i < number; i++) maxima.map[2][i] = -1.0;

    for (int m = 0; m < mx.height; m++) {
        for (int n = 0; n < mx.width; n++) {
            double strength = mx.map[m][n];
            if (strength <= 0) continue;
            if (!isLocalMaximum(mx, m, n)) continue;

            int tooClose = 0;
            for (int i = 0; i < number && !tooClose; i++) {
                if (maxima.map[2][i] < 0) break;
                double dm = m - maxima.map[0][i];
                double dn = n - maxima.map[1][i];
                if (dm*dm + dn*dn < minSepSq) tooClose = 1;
            }
            if (tooClose) continue;

            // Replace weakest entry if this peak is stronger
            int weakest = 0;
            for (int i = 1; i < number; i++)
                if (maxima.map[2][i] < maxima.map[2][weakest]) weakest = i;
            if (strength > maxima.map[2][weakest])
                insertMaxEntry(maxima, m, n, strength);
        }
    }
    return maxima;
}

// Line Rendering
// Convert each peak back to image-space and draw the line
void drawDetectedLines(Image resultImg, Matrix maxima,
                       int mapH, int mapW, double maxD) {
    int H = resultImg.height, W = resultImg.width;
    int numLines = maxima.width;

    unsigned char lineColours[5][3] = {
        {255,   0,   0},   // red
        {  0, 220,   0},   // green
        {  0,  80, 255},   // blue
        {255, 200,   0},   // yellow
        {200,   0, 200},   // magenta
    };

    for (int li = 0; li < numLines; li++) {
        if (maxima.map[2][li] < 0) continue;   // empty slot

        int    alphaIdx = (int)maxima.map[0][li];
        int    dIdx     = (int)maxima.map[1][li];
        double alpha    = -0.5*M_PI + M_PI*(double)alphaIdx/(double)mapH;
        double d        = (double)dIdx / (double)(mapW - 1) * maxD;
        double sa = sin(alpha), ca = cos(alpha);

        unsigned char r = lineColours[li % 5][0];
        unsigned char g = lineColours[li % 5][1];
        unsigned char b = lineColours[li % 5][2];

        printf("  Line %d: alpha=%.1f°  d=%.1f px  votes=%.0f\n",
               li+1, alpha*180.0/M_PI, d, maxima.map[2][li]);

        // Sweep column-by-column when |cos(alpha)| is large (nearly horizontal lines), otherwise sweep row-by-row (nearly vertical lines)
        if (fabs(ca) >= fabs(sa)) {
            for (int col = 0; col < W; col++) {
                if (fabs(sa) < 1e-10) {
                    // Horizontal line: draw the entire row at d
                    int row = (int)(d + 0.5);
                    for (int t = -1; t <= 1; t++) {
                        int rr = row + t;
                        if (rr < 0 || rr >= H) continue;
                        resultImg.map[rr][col].r = r;
                        resultImg.map[rr][col].g = g;
                        resultImg.map[rr][col].b = b;
                    }
                    break;
                }
                int row = (int)((d - col * ca) / sa + 0.5);
                for (int t = -1; t <= 1; t++) {
                    int rr = row + t;
                    if (rr < 0 || rr >= H || col < 0 || col >= W) continue;
                    resultImg.map[rr][col].r = r;
                    resultImg.map[rr][col].g = g;
                    resultImg.map[rr][col].b = b;
                }
            }
        } else {
            // row is the free variable: col = (d - row*sin) / cos
            for (int row = 0; row < H; row++) {
                if (fabs(ca) < 1e-10) {
                    // Vertical line: draw entire column at d
                    int col = (int)(d + 0.5);
                    for (int t = -1; t <= 1; t++) {
                        int cc = col + t;
                        if (cc < 0 || cc >= W) continue;
                        resultImg.map[row][cc].r = r;
                        resultImg.map[row][cc].g = g;
                        resultImg.map[row][cc].b = b;
                    }
                    break;
                }
                int col = (int)((d - row * sa) / ca + 0.5);
                for (int t = -1; t <= 1; t++) {
                    int cc = col + t;
                    if (row < 0 || row >= H || cc < 0 || cc >= W) continue;
                    resultImg.map[row][cc].r = r;
                    resultImg.map[row][cc].g = g;
                    resultImg.map[row][cc].b = b;
                }
            }
        }
    }
}

// For each edge pixel and each candidate radius, vote for all possible circle centres on the circumference
Hough3D houghTransformCircles(Matrix edges, int rMin, int rMax, int rStep) {
    int H = edges.height, W = edges.width;
    int numR = (rMax - rMin) / rStep + 1;
    Hough3D acc = createHough3D(H, W, numR);

    // Precompute sin/cos table at 1 degree resolution
    double *sinA = (double *)malloc(360 * sizeof(double));
    double *cosA = (double *)malloc(360 * sizeof(double));
    for (int a = 0; a < 360; a++) {
        double theta = a * M_PI / 180.0;
        sinA[a] = sin(theta); cosA[a] = cos(theta);
    }

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            if (edges.map[y][x] < 128.0) continue;
            for (int ri = 0; ri < numR; ri++) {
                int r = rMin + ri * rStep;
                for (int a = 0; a < 360; a++) {
                    int cy = (int)(y - r * sinA[a] + 0.5);
                    int cx = (int)(x - r * cosA[a] + 0.5);
                    if (cy >= 0 && cy < H && cx >= 0 && cx < W)
                        acc.data[h3Idx(&acc, cy, cx, ri)]++;
                }
            }
        }
    }
    free(sinA); free(cosA);
    return acc;
}

// projectHough3D collapses the radius axis by taking the max-vote at each (cy, cx)
Matrix projectHough3D(Hough3D *acc) {
    Matrix proj = createMatrix(acc->H, acc->W);
    for (int y = 0; y < acc->H; y++)
        for (int x = 0; x < acc->W; x++) {
            double best = 0;
            for (int ri = 0; ri < acc->R; ri++) {
                int v = acc->data[h3Idx(acc, y, x, ri)];
                if (v > best) best = v;
            }
            proj.map[y][x] = best;
        }
    return proj;
}

// findCircleMaxima returns the <number> highest well-separated local maxima in the Hough accumulator
Matrix findCircleMaxima(Hough3D *acc, int rMin, int rStep,
                        int number, double minSep) {
    Matrix maxima = createMatrix(4, number);
    for (int i = 0; i < number; i++) maxima.map[3][i] = -1.0;
    double minSepSq = minSep * minSep;

    for (int y = 0; y < acc->H; y++) {
        for (int x = 0; x < acc->W; x++) {
            for (int ri = 0; ri < acc->R; ri++) {
                double v = (double)acc->data[h3Idx(acc, y, x, ri)];
                if (v <= 0) continue;

                // 3×3×3 local maximum check
                int isMax = 1;
                for (int dy=-1; dy<=1 && isMax; dy++)
                for (int dx=-1; dx<=1 && isMax; dx++)
                for (int dr=-1; dr<=1 && isMax; dr++) {
                    int ny=y+dy, nx=x+dx, nr=ri+dr;
                    if (ny<0||ny>=acc->H||nx<0||nx>=acc->W||nr<0||nr>=acc->R) continue;
                    if (acc->data[h3Idx(acc,ny,nx,nr)] > (int)v) isMax = 0;
                }
                if (!isMax) continue;

                int tooClose = 0;
                for (int i = 0; i < number && !tooClose; i++) {
                    if (maxima.map[3][i] < 0) break;
                    double dy2 = y - maxima.map[0][i];
                    double dx2 = x - maxima.map[1][i];
                    if (dy2*dy2 + dx2*dx2 < minSepSq) tooClose = 1;
                }
                if (tooClose) continue;

                int weakest = 0;
                for (int i = 1; i < number; i++)
                    if (maxima.map[3][i] < maxima.map[3][weakest]) weakest = i;
                if (v > maxima.map[3][weakest]) {
                    maxima.map[0][weakest] = y;
                    maxima.map[1][weakest] = x;
                    maxima.map[2][weakest] = rMin + ri * rStep;
                    maxima.map[3][weakest] = v;
                }
            }
        }
    }

    // Sort descending by strength (bubble sort: small N)
    for (int i = 0; i < number-1; i++)
        for (int j = i+1; j < number; j++)
            if (maxima.map[3][j] > maxima.map[3][i])
                for (int row = 0; row < 4; row++) {
                    double tmp = maxima.map[row][i];
                    maxima.map[row][i] = maxima.map[row][j];
                    maxima.map[row][j] = tmp;
                }
    return maxima;
}

// Draw ellipses onto resultImg using ellipse()
void drawDetectedCircles(Image resultImg, Matrix maxima) {
    unsigned char circColours[5][3] = {
        {255,   0,   0},
        {  0, 220,   0},
        {  0,  80, 255},
        {255, 200,   0},
        {200,   0, 200},
    };
    for (int i = 0; i < maxima.width; i++) {
        if (maxima.map[3][i] < 0) continue;
        int cy = (int)maxima.map[0][i];
        int cx = (int)maxima.map[1][i];
        int r  = (int)maxima.map[2][i];
        unsigned char rc = circColours[i % 5][0];
        unsigned char gc = circColours[i % 5][1];
        unsigned char bc = circColours[i % 5][2];
        printf("  Circle %d: center=(%d, %d)  radius=%d  votes=%.0f\n",
               i+1, cx, cy, r, maxima.map[3][i]);
        ellipse(resultImg, cy, cx, r, r, 2, 0, 0, rc, gc, bc, 0);
    }
}

// runLinePipeline complete line detection pipeline for one image file
//   1. Load image and build edge matrix
//   2. Compute line Hough accumulator
//   3. Visualise the accumulator as a scaled image
//   4. Pick the top numLines peaks
//   5. Render detected lines onto a copy of the input
//   6. Write all outputs to outDir
void runLinePipeline(const char *inputFile, const char *outDir,
                     int numLines, double edgeRatio) {
    printf("\n=== Line Hough: %s ===\n", inputFile);

    Image img = readImage((char *)inputFile);
    printf("  Image: %d × %d\n", img.width, img.height);

    // Edge map
    printf("  Building edge matrix (ratio=%.2f)...\n", edgeRatio);
    Matrix edges = buildEdgeMatrix(img, edgeRatio);

    // Edge image
    Image edgeImg = createImage(img.height, img.width);
    matrixToImageScaled(&edges, &edgeImg);
    char path[512];
    snprintf(path, sizeof(path), "%s/edges.pgm", outDir);
    writeImage(edgeImg, path);
    printf("  Edge map saved -> %s\n", path);

    // Hough accumulator
    double maxD = sqrt((double)(img.height*img.height + img.width*img.width));
    printf("  Computing line Hough transform (%d × %d bins)...\n",
           LINE_MAP_H, LINE_MAP_W);
    Matrix houghMx = houghTransformLines(edges, LINE_MAP_H, LINE_MAP_W);

    // Accumulator image
    Image houghImg = createImage(LINE_MAP_H, LINE_MAP_W);
    matrixToImageScaled(&houghMx, &houghImg);
    snprintf(path, sizeof(path), "%s/hough_accumulator.pgm", outDir);
    writeImage(houghImg, path);
    printf("  Accumulator saved -> %s\n", path);

    // Peak picking
    printf("  Finding top %d peaks (minSep=%.0f)...\n", numLines, LINE_MIN_SEP);
    Matrix maxima = findHoughMaxima(houghMx, numLines, LINE_MIN_SEP);

    // Accumulator image with peaks marked (white dots)
    Image houghAnnotated = createImage(LINE_MAP_H, LINE_MAP_W);
    matrixToImageScaled(&houghMx, &houghAnnotated);
    for (int i = 0; i < numLines; i++) {
        if (maxima.map[2][i] < 0) continue;
        int ar = (int)maxima.map[0][i];
        int ac = (int)maxima.map[1][i];
        // Draw a small 5×5 white square around each peak
        for (int dy = -2; dy <= 2; dy++)
            for (int dx = -2; dx <= 2; dx++) {
                int rr = ar+dy, cc = ac+dx;
                if (rr<0||rr>=LINE_MAP_H||cc<0||cc>=LINE_MAP_W) continue;
                houghAnnotated.map[rr][cc].r = 255;
                houghAnnotated.map[rr][cc].g = 0;
                houghAnnotated.map[rr][cc].b = 0;
            }
    }
    snprintf(path, sizeof(path), "%s/hough_peaks.ppm", outDir);
    writeImage(houghAnnotated, path);
    printf("  Annotated accumulator -> %s\n", path);

    // Result image with lines drawn
    Image resultImg = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            resultImg.map[y][x] = img.map[y][x];

    printf("  Drawing %d detected lines:\n", numLines);
    drawDetectedLines(resultImg, maxima, LINE_MAP_H, LINE_MAP_W, maxD);

    snprintf(path, sizeof(path), "%s/detected_lines.ppm", outDir);
    writeImage(resultImg, path);
    printf("  Result image -> %s\n", path);

    // Cleanup
    deleteMatrix(edges); deleteMatrix(houghMx); deleteMatrix(maxima);
    deleteImage(edgeImg); deleteImage(houghImg);
    deleteImage(houghAnnotated); deleteImage(resultImg);
    deleteImage(img);
}

// runCirclePipeline complete circle detection pipeline for one image file
void runCirclePipeline(const char *inputFile, const char *outDir,
                       int numCircles, double edgeRatio) {
    printf("\n=== Circle Hough: %s ===\n", inputFile);

    Image img = readImage((char *)inputFile);
    printf("  Image: %d × %d\n", img.width, img.height);

    // Edge map (slightly lower threshold: circles need more edge coverage)
    double circEdgeRatio = edgeRatio * 0.8;
    printf("  Building edge matrix (ratio=%.2f)...\n", circEdgeRatio);
    Matrix edges = buildEdgeMatrix(img, circEdgeRatio);

    Image edgeImg = createImage(img.height, img.width);
    matrixToImageScaled(&edges, &edgeImg);
    char path[512];
    snprintf(path, sizeof(path), "%s/circle_edges.pgm", outDir);
    writeImage(edgeImg, path);

    // 3-D Hough accumulator
    printf("  Computing circle Hough transform"
           " (r=%d..%d step=%d)...\n", CIRCLE_R_MIN, CIRCLE_R_MAX, CIRCLE_R_STEP);
    Hough3D acc = houghTransformCircles(edges, CIRCLE_R_MIN, CIRCLE_R_MAX, CIRCLE_R_STEP);

    // 2-D projection for visualisation (heatmap)
    Matrix proj = projectHough3D(&acc);
    Image projImg = createImage(img.height, img.width);
    matrixToImageScaled(&proj, &projImg);
    snprintf(path, sizeof(path), "%s/circle_accumulator.pgm", outDir);
    writeImage(projImg, path);
    printf("  Circle accumulator projection -> %s\n", path);

    // Peak picking
    printf("  Finding top %d circles (minSep=%.0f px)...\n",
           numCircles, CIRCLE_MIN_SEP);
    Matrix maxima = findCircleMaxima(&acc, CIRCLE_R_MIN, CIRCLE_R_STEP,
                                     numCircles, CIRCLE_MIN_SEP);

    // Result image with circles drawn
    Image resultImg = createImage(img.height, img.width);
    for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++)
            resultImg.map[y][x] = img.map[y][x];

    printf("  Drawing %d detected circles:\n", numCircles);
    drawDetectedCircles(resultImg, maxima);

    snprintf(path, sizeof(path), "%s/detected_circles.ppm", outDir);
    writeImage(resultImg, path);
    printf("  Result image -> %s\n", path);

    deleteMatrix(edges); deleteMatrix(proj); deleteMatrix(maxima);
    deleteImage(edgeImg); deleteImage(projImg); deleteImage(resultImg);
    deleteHough3D(acc);
    deleteImage(img);
}

int main(int argc, char *argv[]) {
    // Single-image mode
    if (argc > 1) {
        char  *inputFile = argv[1];
        char  *mode      = (argc > 2) ? argv[2] : "lines";
        int    numPeaks  = (argc > 3) ? atoi(argv[3]) : 5;
        double edgeRatio = (argc > 4) ? atof(argv[4]) : 0.15;
#ifdef _WIN32
        mkdir("output_images");
#else
        mkdir("output_images", 0755);
#endif
        if (strcmp(mode, "circles") == 0)
            runCirclePipeline(inputFile, "output_images", numPeaks, edgeRatio);
        else
            runLinePipeline  (inputFile, "output_images", numPeaks, edgeRatio);
        printf("\nDone. All outputs in output_images/\n");
        return 0;
    }

    // Default demo mode: run all scenarios, separate output dirs
#ifdef _WIN32
    mkdir("output_images");
    mkdir("output_images/coastline_lines");
    mkdir("output_images/tahoe_lines");
    mkdir("output_images/tahoe_circles");
#else
    mkdir("output_images",                 0755);
    mkdir("output_images/coastline_lines", 0755);
    mkdir("output_images/tahoe_lines",     0755);
    mkdir("output_images/tahoe_circles",   0755);
#endif

    // Coastline: line detection (8 most voted lines, edge ratio 0.15)
    runLinePipeline("coastline.ppm",
                    "output_images/coastline_lines", 8, 0.15);

    // Tahoe: line detection (5 lines, slightly higher ratio for dense texture: circles need more edge coverage)
    runLinePipeline("tahoe.ppm",
                    "output_images/tahoe_lines", 5, 0.18);

    // Tahoe: circle detection (3 circles: looking for lake boundary)
    runCirclePipeline("tahoe.ppm",
                      "output_images/tahoe_circles", 3, 0.12);

    printf("\nDone. All outputs in output_images/\n");
    return 0;
}

// compile: gcc hough.c netpbm.c -o hough -lm
//
// Example: detect 8 coastline lines:
//   hough.exe ../input_images/coastline.ppm lines 8 0.15
//
// Example: detect lake boundary circles:
//   hough.exe ../input_images/tahoe.ppm circles 3 0.12
