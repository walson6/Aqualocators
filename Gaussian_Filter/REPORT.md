# Gaussian Filter Report

## Goal

This module smooths satellite images before edge detection or texture segmentation is applied. The goal is to reduce high-frequency noise and fine surface detail that would otherwise cause false edge responses or unstable texture features in later processing stages. Three filter variants are implemented and compared: a uniform averaging filter, a true Gaussian filter, and a median filter. All three share the same interface and output format so their results can be compared directly.

## Method

The program reads a PPM image and processes the grayscale intensity channel. Each filter operates by sliding a kernel window across the image and replacing each pixel value with a weighted or ordered combination of its neighbors.

Uniform Averaging Filter

A 5x5 kernel filled with the value 1/25 is convolved with the image. Every neighbor within the window contributes equally to the output. This is the simplest possible smoothing operation and serves as a baseline for comparison. It reduces noise effectively but also blurs edges more aggressively than the other two filters because it does not weight closer neighbors more heavily.

True Gaussian Filter

A 5x5 kernel whose weights follow the two-dimensional Gaussian function with sigma approximately 1.4 is convolved with the image. The unnormalized weights are taken from the pattern used inside the Canny pipeline so that the standalone filter and the Canny pre-smoother produce identical results when given the same input. The weights are normalized to sum to 1 before convolution so that the mean intensity of the image is preserved. Pixels closer to the center of the kernel receive higher weight than pixels at the edges, which means the filter smooths noise while preserving relatively more edge contrast than the uniform average.

Median Filter

A 5x5 window is placed at each pixel and all 25 values within the window are sorted. The output at that pixel is set to the median value. The median filter does not convolve a kernel at all. Because it outputs an actual pixel value that existed in the neighborhood rather than a weighted average, it preserves edges more sharply than linear filters. It is effective at removing salt-and-pepper noise, which the Gaussian and averaging filters tend to spread rather than eliminate.

Convolution Boundary Handling

All three filters use zero-padding at the image border. Pixels within half the kernel size of the image edge are set to zero in the output rather than attempting to compute a partial sum. This creates a border in the output images.

## Parameters

Default parameters:

- Kernel size = 5x5 for all three filters
- Gaussian sigma = approximately 1.4
- Median window = 5x5 (25 values, median at index 12)
- Averaging weight = 1/25 per cell

Suggested starting commands:

```bash
gcc gaussian.c ../netpbm.c -o gaussian -lm
./gaussian ../input_images/coastline.ppm
./gaussian ../input_images/tahoe.ppm
```

The kernel size is fixed at 5x5 in the current implementation. To experiment with stronger smoothing, the blur passes inside the robustness module can be used to apply the Gaussian kernel multiple times in sequence, which approximates a larger effective sigma without changing the kernel definition.

## Outputs

The program saves all results in output_images/:

- averaging_smoothed.ppm - result of the 5x5 uniform averaging filter
- gaussian_smoothed.ppm - result of the 5x5 Gaussian filter
- median_filtered.ppm - result of the 5x5 median filter

## Observations

On the coastline image all three filters reduce the fine surface detail visible in the farm fields and surf zone. The uniform averaging filter produces a visibly blurrier result than the Gaussian, with the coastline boundary becoming slightly less sharp. The Gaussian filter preserves the water-land boundary more cleanly because the center-weighted kernel applies less smoothing to the transition region. The median filter produces the sharpest coastline edge of the three while still reducing the texture noise in the inland fields, which makes it the most attractive preprocessing step when edge localization accuracy matters.

The median filter is slower than the Gaussian because it requires sorting 25 values per pixel rather than performing a simple weighted sum. The median filter's edge-preserving property makes it useful as a denoising step when the future detector is Sobel, which lacks any built-in smoothing. When the future detector is Canny, the internal 5x5 Gaussian smoother already provides most of the noise suppression needed, so adding an external median filter before Canny provides diminishing returns but can be good for images like Tahoe where there are lots of sharp rock textures.

## Strengths and Weaknesses

Strengths:

- All three filters reduce noise that would otherwise cause false edge responses in Sobel and false cluster assignments in K-means texture segmentation
- The Gaussian filter matches the internal smoother used by the Canny pipeline, so results are predictable and consistent when comparing preprocessed and unprocessed inputs
- The median filter preserves edges more sharply than either linear filter, which benefits edge localization accuracy

Weaknesses:

- Zero-padding at the image border creates a border each output, which can inflate edge pixel counts under heavy blur if not masked
- Repeated Gaussian passes double-smooth the image when used before Canny, since Canny applies its own internal Gaussian, which can over-blur weak but useful edge features
- All three filters use a fixed 5x5 kernel with no command-line option to change the kernel size or sigma, requiring direct code edits to experiment with larger windows
- Linear filters spread salt-and-pepper noise across neighboring pixels rather than eliminating it, making the median filter the better choice when impulse noise is present