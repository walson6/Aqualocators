# Hough Transform Report

## Goal

This module detects dominant geometric shapes in images using the Hough Transform. Two modes are implemented: line detection and circle detection. The line detector is used to find straight coastline segments and field boundaries in the coastline image. The circle detector is used to find the roughly circular boundary of the lake in the tahoe image. Both modes operate on a binary edge map produced by Gaussian smoothing followed by Sobel thresholding.

## Method

The Hough Transform works by converting the edge detection problem into a voting problem. Instead of asking which pixels belong to a line, it asks each edge pixel to vote for every line or circle that could pass through it. Shapes that pass through many edge pixels accumulate many votes, and the peaks in the vote accumulator correspond to the most likely shapes in the image.

Edge Preprocessing

Before Hough transform runs, the input image is processed through a shared preprocessing function. A 3x3 Gaussian kernel is first convolved with the intensity channel to remove pixel-level noise. The Sobel magnitude is then computed and pixels above a fraction of the maximum gradient are marked as edge pixels. A 10-pixel border strip around the image is excluded from the edge map to prevent the image frame itself from flooding the accumulator with votes near distance zero. The fraction is controlled by the edgeThreshRatio parameter, with a default of 0.15 for lines and 0.12 for circles.

Line Detection

Lines are represented in normal form as d = row times sin(alpha) plus col times cos(alpha). This parameteric form avoids the infinite slope problem that occurs when using slope-intercept form for vertical lines.

The accumulator is a 2D matrix with 360 rows for the angle axis and 480 columns for the distance axis. For each edge pixel, the distance d is computed at every angle and the corresponding accumulator cell is incremented. After all edge pixels have voted, the accumulator is searched for local maxima using an 8-neighbor comparison. The top N peaks are selected with a minimum separation constraint to prevent the same physical line from being detected multiple times at slightly different parameter values. Each detected line is then projected back into image space and drawn across the full image width or height in a distinct color.

Circle Detection

Circles are parameterized by their center coordinates (cy, cx) and radius r. The accumulator is a 3D array indexed by these three values. For each edge pixel and each candidate radius in the search range, the pixel votes for all possible circle centers that would place it on the circumference. This is done by sweeping 360 angles and computing the corresponding center position for each angle at the current radius.

After accumulation, the 3D array is projected down to a 2D heatmap by taking the maximum vote across all radii at each center position. This projection is saved as a visualization. Peak picking then searches the full 3D accumulator for local maxima in a 3x3x3 neighborhood. The top N peaks with sufficient spatial separation are kept and drawn onto the result image as ellipses.

Peak Picking

Both line and circle modes use the same general peak picking strategy. A candidate peak must be a strict local maximum in its neighborhood. It must also be farther than a minimum separation distance from any previously accepted peak. If a new candidate is stronger than the weakest currently accepted peak and is not too close to any existing peak, it replaces the weakest entry. This produces the N strongest well-separated detections.

## Parameters

Default parameters for line detection:

- Accumulator size = 360 x 480 (angle x distance bins)
- Number of peaks = 5
- Edge threshold ratio = 0.15
- Minimum peak separation = 20 accumulator bins

Default parameters for circle detection:

- Radius search range = 15 to 300 pixels, step 3
- Number of circles = 3
- Edge threshold ratio = 0.12 (slightly lower to capture more edge coverage)
- Minimum peak separation = 40 spatial pixels

Suggested starting commands:

```bash
gcc hough.c ../netpbm.c -o hough -lm
./hough ../input_images/coastline.ppm lines 8 0.15
./hough ../input_images/tahoe.ppm circles 3 0.12
./hough ../input_images/coastline.ppm lines 5 0.20
```

The number of peaks and the edge threshold ratio are the two most useful parameters to adjust. Raising the threshold ratio reduces the number of edge pixels that vote, which cleans up the accumulator but risks missing weaker boundaries. Raising the number of peaks allows more lines or circles to be detected but increases the chance of including spurious detections.

## Outputs

The program saves all results in output_images/:

For line detection:
- edges.pgm - binary edge map used as input to the accumulator
- hough_accumulator.pgm - grayscale visualization of the 2D vote accumulator
- hough_peaks.ppm - accumulator image with the detected peaks marked as small red squares
- detected_lines.ppm - copy of the original image with the detected lines drawn in color

For circle detection:
- circle_edges.pgm - binary edge map used for circle voting
- circle_accumulator.pgm - 2D projection of the 3D accumulator (max vote per center position)
- detected_circles.ppm - copy of the original image with the detected circles drawn in color

## Observations

On the coastline image the line detector identifies the dominant horizontal boundary between the dark water region and the lighter land. At 8 peaks with a threshold ratio of 0.15 the top detections include the main surf-zone boundary, the beach edge, the upper field boundary running roughly parallel to the coast, and several of the diagonal road or fence lines in the upper portion of the image.

On the tahoe image the line detector has difficulty finding useful lines. These are less useful for lake boundary detection than the circle detector.

The circle detector on the tahoe image does not produce useful results. Because the lake outline is not a perfect circle, the accumulator peak is broader and lower than it would be for a manufactured circular object. Lowering the edge threshold ratio to 0.10 or 0.12 increases the number of edge pixels on the lake shore and raises the peak vote count, which can improve detection confidence. The smaller ponds visible in the upper portion of the tahoe image are also detected as secondary peaks when the number of circles is set to 3 or more.

## Strengths and Weaknesses

Strengths:

- The normal-form line parameterization handles all line orientations including vertical without any special cases or infinite slope issues
- The minimum separation constraint in peak picking prevents duplicate detections of the same physical line from nearby accumulator cells
- The circle detector searches a continuous radius range, which means it does not require prior knowledge of the object size and can detect the lake boundary even when its exact radius is unknown
- Separate output images for the accumulator, the annotated accumulator, and the final result make it easy to diagnose why a particular detection succeeded or failed

Weaknesses:

- Circle detection is slow because every edge pixel votes for every radius and every angle
- The line detector performs best on images with long straight edges. The irregular, curved coastline in the coastline image is not easily identified by a small number of lines, and the Hough line result captures only the approximate orientation rather than the true boundary shape
- Both detectors are sensitive to the edge threshold ratio. A ratio that is too low floods the accumulator with weak edges and produces unreliable peaks. A ratio that is too high misses important boundary pixels and drops the vote count below the detection threshold

## Takeaway
- The line detector works well for finding straight features like roads and field boundaries, but struggles with uneven water outlines.
- The circle detector in practice struggles when the shape is not perfectly circular, like the irregular outline of Lake Tahoe.
- Results depend heavily on tuning parameters like edge threshold and number of peaks. Lots of tuning needs to be done in order to achieve the desired results.
- Line detection is generally more reliable than circle detection in real-world satellite images because natural features are rarely perfectly circle.