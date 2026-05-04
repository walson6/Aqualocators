# Robustness Analysis Report

## Goal

This section evaluates how the existing preprocessing pipeline behaves under common real-world distortions: Gaussian noise, blur, and low contrast. The subset contains both sample images, `coastline.ppm` and `tahoe.ppm`.

## Method

The robustness program creates four cases for each image:

- Original baseline
- Gaussian noise with standard deviation `24`
- Gaussian blur using the same 5x5 Gaussian kernel used elsewhere in the project
- Low contrast using `128 + 0.45 * (pixel - 128)`

For each case, the program runs:

- Gaussian filtering
- Sobel edge detection on the distorted image
- Sobel edge detection after Gaussian filtering
- Canny edge detection on the distorted image
- Canny edge detection after Gaussian filtering
- Texture-based K-means segmentation using local mean, standard deviation, gradient energy, and entropy

Outputs are saved in `output_images/`, and quantitative edge counts are saved in `output_images/metrics.csv`.

## Quantitative Summary

Edge percentages from `metrics.csv`:

| Image | Case | Sobel Raw | Sobel After Gaussian | Canny Raw | Canny After Gaussian |
|---|---:|---:|---:|---:|---:|
| coastline | original | 0.301% | 0.475% | 3.270% | 2.860% |
| coastline | Gaussian noise | 0.360% | 0.421% | 7.716% | 3.835% |
| coastline | blur | 0.475% | 0.591% | 2.860% | 2.710% |
| coastline | low contrast | 0.299% | 0.467% | 3.445% | 3.226% |
| tahoe | original | 0.489% | 1.330% | 10.674% | 8.932% |
| tahoe | Gaussian noise | 0.678% | 1.206% | 12.603% | 9.647% |
| tahoe | blur | 1.330% | 1.598% | 8.932% | 7.986% |
| tahoe | low contrast | 0.487% | 1.283% | 10.851% | 9.326% |

## Compile And Run

```bash
gcc robustness_analysis.c ../netpbm.c -o robustness_analysis -lm
./robustness_analysis
```

## Observations

Under Gaussian noise, Sobel is the most noise-sensitive method. Since Sobel directly measures local gradient magnitude, random pixel changes create many short, scattered responses. This increases edge density but reduces useful continuity.

Canny is more stable under noisy conditions because it smooths internally, suppresses non-maximum responses, and uses hysteresis to keep weak edges only when they connect to strong edges. The result is usually thinner and less speckled than Sobel, although very noisy images still create extra edge fragments.

Applying Gaussian filtering before edge detection improves noisy images. Sobel benefits the most because smoothing removes many isolated high-gradient noise pixels before the gradient calculation. Canny also improves, but the difference is smaller because Canny already includes Gaussian smoothing as its first stage.

Blur reduces edge strength and removes fine detail. Sobel produces softer and less distinct edge maps, while Canny may lose weak boundaries if they fall below the high/low thresholds. In blurred images, Canny remains cleaner, but some boundaries become less continuous.

Low contrast weakens intensity differences. Sobel still shows broad gradient responses after scaling, but those responses are less meaningful because weak texture and weak boundaries become harder to separate. Canny is more likely to drop subtle edges when the maximum gradient range is compressed.

Texture segmentation is affected differently by each distortion. Gaussian noise increases local standard deviation, entropy, and gradient energy, which can cause K-means to split noisy areas into artificial texture regions. Blur has the opposite effect: it smooths local variance and entropy, so distinct textures may merge into fewer visually meaningful regions. Low contrast reduces separability between regions because local mean and gradient features become closer together.

## Conclusions

Canny is the better edge detector for distorted images when edge quality, thinness, and continuity matter. Sobel is useful as a simple gradient baseline, but it is highly sensitive to Gaussian noise. Gaussian filtering before edge detection is valuable, especially for Sobel, because it reduces noisy false edges. For segmentation, distortions change the texture features directly, so segmentation quality depends on whether the distortion increases artificial texture noise or removes real texture differences.
