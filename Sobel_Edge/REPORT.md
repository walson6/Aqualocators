# Sobel Edge Detection Report

## Goal

This module detects edges in images using the Sobel gradient operator. The goal is to produce a gradient magnitude map that highlights regions of rapid intensity change, which correspond to boundaries between land and water, field edges, roads, and rocky terrain features. Sobel is used here as the edge detector and as a reference baseline for comparing against the Canny edge detector. It is also used internally by the Hough transform preprocessing stage.

## Method

The program reads a PPM image and processes the grayscale intensity channel through three stages:

Stage 1 - Intensity Extraction

The intensity channel of the input image is copied into a floating-point matrix so that intermediate gradient values are not truncated during computation. For color PPM inputs the intensity field already holds the average of the red, green, and blue channels as computed by the NetPBM loader, so no additional channel mixing is needed.

Stage 2 - Gradient Computation

Two 3x3 kernels are convolved with the intensity matrix separately. The horizontal kernel Gx detects vertical edges by computing the difference between the right and left columns of the neighborhood, with the center row weighted by a factor of two to reduce sensitivity to diagonal structure. The vertical kernel Gy detects horizontal edges in the same way along the column direction. The kernels used are the standard Sobel operators:

Gx:
-1  0  1
-2  0  2
-1  0  1

Gy:
-1 -2 -1
 0  0  0
 1  2  1

The gradient magnitude at each pixel is computed as the square root of Gx squared plus Gy squared. This gives a single value at every pixel that represents how strongly the intensity is changing in any direction.

Stage 3 - Output Scaling

The gradient magnitude matrix contains values ranging from zero up to several hundred, depending on the contrast of the input image. These values are normalized to the range 0 to 255 so that the weakest gradient in the image maps to black and the strongest maps to white. The normalized values are written to all four pixel channels (r, g, b, i) so the output displays correctly as a grayscale image.

## Parameters

Default parameters:

- Kernel size = 3x3
- Gradient magnitude = sqrt(Gx squared plus Gy squared)
- Output scaling = linear normalization to 0 to 255 range

Suggested starting commands:

```bash
gcc sobel.c ../netpbm.c -o sobel -lm
./sobel ../input_images/coastline.ppm
./sobel ../input_images/tahoe.ppm
```

## Outputs

The program saves all results in output_images/:

- sobel_edges.ppm - color copy of the gradient magnitude map
- sobel_edges.pgm - grayscale copy

## Observations

On the coastline image the Sobel output shows a strong bright line along the water-land boundary where the darker water meets the brighter land. The surf, where white foam produces sharp local contrast against both the water and the land, generates some of the highest gradient values in the entire image. The inland farm fields produce moderate gradient responses at their boundaries, and the roads visible in the upper portion of the image show as thin bright lines. 

On the tahoe image the gradient magnitude map is far denser. The rocky terrain covering most of the image produces a high density of moderate-strength gradient responses because the rock surface has many small brightness transitions at the scale of the 3x3 kernel. The lake boundary generates a strong response where the dark water meets the lighter rock, but this response is mixed in with the surrounding rock texture rather than standing out cleanly. The measured baseline edge pixel count on tahoe is 7198, which is higher than the coastline count of 4275 despite tahoe having a larger image area, reflecting the higher texture density of the rocky landscape.

The global output scaling behavior of Sobel has an important practical consequence for robustness. Because the entire gradient magnitude range is always mapped to 0 to 255 regardless of the absolute scale of the gradients, reducing the image contrast by a factor of five also reduces all gradient magnitudes by the same factor, but after rescaling the relative distribution of pixel values in the output is identical to the original. The robustness measurements confirm this directly: the Sobel edge pixel count stays within 0.993 to 1.013 times the baseline across all four low-contrast levels tested, including an 80 percent dynamic range compression. This makes Sobel almost completely immune to global contrast changes, which is a useful property for satellite imagery affected by haze or overcast lighting. However, this invariance is a property of the normalization step rather than of the gradient computation itself, and it means that Sobel cannot distinguish between a high-contrast image and a low-contrast image once the output has been scaled.

Comparing Sobel with Canny on the same inputs, Sobel detects fewer edge pixels on the clean coastline image (4275 versus 13394) because it has no hysteresis linking step to extend strong edges into adjacent weak-response regions. On the other hand, Sobel produces a more stable count under noise because the global rescaling absorbs changes in the overall gradient energy level, whereas Canny's hysteresis chains can break or extend unpredictably when noise shifts individual pixel values across the threshold boundaries.

## Strengths and Weaknesses

Strengths:

- Simple and fast to compute, requiring only two 3x3 convolutions and a magnitude calculation per pixel
- The global output normalization makes the edge map nearly invariant to changes in overall image brightness and contrast, which is valuable for images taken under variable lighting conditions
- Produces a continuous gradient magnitude output rather than a binary map, which preserves information about edge strength that binary detectors discard
- The output can be thresholded at any level after the fact, giving flexibility in downstream processing without requiring recomputation
- Serves as a reliable and interpretable baseline for comparing against more complex detectors

Weaknesses:

- Produces thick edges that are typically two to three pixels wide because the 3x3 kernel responds to a gradient region rather than a precise boundary location. This makes it unsuitable for applications requiring sub-pixel edge localization
- Has no built-in noise suppression. Without a preceding Gaussian smoothing step, salt-and-pepper noise and fine texture produce numerous false edge responses that cannot be distinguished from true boundaries in the output
- Has no thresholding mechanism. The gradient magnitude output requires a separate thresholding step to produce a binary edge map, and the appropriate threshold is not automatically determined
- The global output scaling means that a noisy image and a clean image can produce edge maps that look similar in intensity distribution even when the actual edge quality is very different
- On high-texture images like tahoe the gradient map is dominated by texture noise, and the true boundary of interest such as the lake shoreline is not visually separable from the background without additional processing
- The 3x3 kernel is fixed and operates at a single spatial scale. Fine texture at the kernel scale produces the same response magnitude as a true object boundary, with no mechanism to distinguish between them based on scale alone

## Takeaway

- Sobel is a fast and simple way to detect edges by finding where brightness changes quickly in an image.
- It works well for clear boundaries like coastlines, but it also picks up a lot of unwanted texture and noise in detailed scenes like rocky terrain. But this can be good as water is generally less detailed and appear as dark areas in grayscale images.
- The output shows edge strength, not clean edges, so it often looks thick and messy compared to more advanced methods like Canny. But identifying which areas are water is better identified with Sobel.