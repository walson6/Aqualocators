# Canny Edge Detection Report

## Goal

This module detects edges in satellite images using the Canny edge detector. The goal is to produce thin, well-localized, continuous edge maps that accurately follow coastlines and land-water boundaries. Canny is used here as the primary edge detector because it is designed to suppress noise, reduce false edges, and produce single-pixel-wide edge chains, which makes it far more suitable for coastline tracing than a simple gradient filter.

## Method

The program reads a PPM image and processes the grayscale intensity channel through five stages.

Stage 1 - Gaussian Smoothing

A 5x5 Gaussian kernel with sigma approximately 1.4 is applied to the image before any gradient computation. This removes high-frequency noise that would otherwise create spurious edge responses in later stages. The kernel weights are derived from the standard Gaussian function and normalized.

Stage 2 - Sobel Gradient Computation

The smoothed image is convolved with the standard 3x3 Sobel kernels in both the horizontal and vertical directions. From these two result maps, the gradient magnitude and gradient direction are computed at every pixel using the square root of the sum of squares and the arctangent of the vertical over horizontal result, respectively.

Stage 3 - Non-Maximum Suppression

The gradient direction at each pixel is quantized into one of four bins covering the angles 0, 45, 90, and 135 degrees. For each pixel, the gradient magnitude is compared with the two neighbors that lie along the gradient direction. If the pixel is not a local maximum in that direction, its magnitude is set to zero. This step thins the edges from several pixels wide down to a single pixel.

Stage 4 - Thresholding

Two thresholds are applied to the suppressed gradient map. Any pixel above the high threshold is immediately classified as a strong edge. Any pixel between the low and high thresholds is classified as a weak edge. Any pixel below the low threshold is suppressed entirely. The default high threshold is 15 percent of the maximum suppressed magnitude and the default low threshold is 50 percent of the high threshold. Both values can be overridden on the command line.

Stage 5 - Hysteresis Edge Linking

Weak edge pixels are kept only if they are connected to at least one strong edge pixel through an 8-connected neighborhood. The algorithm alternates between forward and backward passes through the image until no more changes occur. Any remaining weak pixels that were never connected to a strong edge are discarded.

## Parameters

Default parameters:

- High threshold ratio = 0.15 (15 percent of the maximum NMS gradient magnitude)
- Low threshold ratio = 0.50 (50 percent of the high threshold)
- Gaussian kernel size = 5x5
- Gaussian sigma = approximately 1.4

Suggested starting commands:

```bash
gcc canny.c ../netpbm.c -o canny -lm
./canny ../input_images/coastline.ppm
./canny ../input_images/coastline.ppm 0.20 0.40
./canny ../input_images/tahoe.ppm 0.12 0.50
```

Lowering the high threshold detects more edges but risks including weak noise. Raising it keeps only the strongest edges. Lowering the low ratio allows more weak edges to be linked to strong ones, which helps close gaps in the coastline boundary.

## Outputs

The program saves all results in output_images/:

- canny_edges.ppm - color copy of the edge map for visual inspection
- canny_edges.pgm - grayscale edge map for quantitative comparison and further processing

Edge pixels in both outputs are set to intensity 255 (white). Non-edge pixels are set to 0 (black). The edge map is binary because the hysteresis step either keeps or discards each pixel with no intermediate values surviving.

## Observations

On the coastline image the detector produces a clean, nearly continuous boundary along the water-land transition. The surf zone, where white foam meets the darker water and brown earth, generates the strongest gradient response and is detected even at conservative threshold settings. Interior field boundaries and road edges are also detected, though they appear as shorter, more fragmented chains than the coastline itself.

On the tahoe image the rocky terrain generates a very high density of edge responses because the surface has many closely spaced brightness transitions. At the default threshold the edge map is dominated by rock texture rather than the lake boundary. Increasing the threshold ratio to around 0.20 and adding Gaussian blur before running the detector can help suppress the texture edges while preserving the strong lake shoreline.

The threshold and hysteresis combination is the main practical advantage over the simpler Sobel detector. Sobel outputs a continuous gradient magnitude map and uses a single cutoff, which means it cannot distinguish between a strong coherent edge and a cluster of moderate noise responses. Canny links weak edges only when they are connected to confirmed strong edges, so isolated noise peaks are discarded even when their magnitude would pass a simple threshold.

## Strengths and Weaknesses

Strengths:

- Produces thin, single-pixel-wide edges that are easier to trace and measure
- Hysteresis linking closes small gaps in coastline edges caused by surf or shadows
- The built-in Gaussian smoother suppresses noise before gradient computation, making the detector more stable than raw Sobel under moderate noise levels
- Threshold separates confirmed edges from uncertain ones, reducing false positives significantly compared to a single-threshold approach
- Both threshold values can be tuned on the command line without recompiling

Weaknesses:

- The internal Gaussian smoother adds to any external blur already applied, which can over-smooth fine boundary detail when used together with preprocessing
- Performance degrades significantly under low-contrast conditions because the hysteresis step relies on the absolute difference between strong and weak edge magnitudes shrinking to near zero as contrast is compressed
- On high-texture images like tahoe the edge map is dense with rock-texture responses unless thresholds are raised or preprocessing blur is added first
- Threshold selection is image-dependent and requires manual tuning or an automatic method such as Otsu thresholding to work reliably across different scenes
- The zero-padding boundary condition in the convolution step creates a faint artificial edge at the 2-pixel image border, though this is typically suppressed by the high threshold

## Takeaway

- Canny is a strong edge detector because it gives clean, thin lines which makes it great for tracing coastlines.
- Canny works better than simple methods like Sobel because it filters out noise first and only keeps edges that are meaningful.
- The threshold system helps a lot because it keeps strong edges and only keeps weaker ones if they’re connected, which removes random noise while still filling small gaps.
- On detailed images like Tahoe's rocky terrain, Canny can pick up too many edges, so you may need to increase thresholds or add more smoothing to focus on the important boundaries.
- Canny is good, but getting the best results depends on choosing the right thresholds for your image.