# Texture-Based Image Segmentation Report

## Goal

This module segments an image into regions based on local texture. It follows the assignment direction to start with grayscale images, then uses K-means clustering to group pixels with similar local texture statistics.

## Method

The program reads a PGM or PPM image through the shared NetPBM loader and uses the grayscale intensity channel. For each pixel, it computes four local texture features inside an odd square window:

- Local mean intensity
- Local standard deviation
- Local Sobel gradient energy
- Local entropy using 16 gray-level bins

Each feature is normalized to `[0, 1]` before clustering so that high-range features do not dominate the K-means distance. K-means then assigns each pixel to one of `K` texture regions.

## Parameters

Default parameters:

- `K = 4` clusters
- `windowSize = 9`
- `iterations = 30`

Suggested starting commands:

```bash
gcc texture_segmentation.c ../netpbm.c -o texture_segmentation -lm
./texture_segmentation ../input_images/coastline.pgm 4 9 30
./texture_segmentation ../input_images/tahoe.pgm 5 11 40
```

Larger windows create smoother, broader texture regions. Smaller windows preserve local detail but can produce noisier region boundaries. Increasing `K` separates more texture types, but too many clusters may split one meaningful region into several visually similar labels.

## Outputs

The program saves all results in `output_images/`:

- `<input>_texture_labels.pgm`: grayscale label mask
- `<input>_texture_regions.ppm`: false-color segmentation
- `<input>_texture_overlay.ppm`: segmentation blended with the original image
- `<input>_feature_mean.pgm`: local average brightness visualization
- `<input>_feature_stddev.pgm`: local contrast/roughness visualization
- `<input>_feature_gradient_energy.pgm`: edge-density texture visualization
- `<input>_feature_entropy.pgm`: local intensity complexity visualization

## Observations

This approach is useful for separating regions that differ by repeated local structure rather than only by brightness. For example, water, shoreline, sky, and land textures can cluster differently even when some regions have overlapping intensity values.

The grayscale version is the best first step for this assignment because it isolates texture behavior from color effects. A later extension could add color features such as local average red, green, and blue values, but that would change the method from texture-first segmentation to combined color-and-texture segmentation.
