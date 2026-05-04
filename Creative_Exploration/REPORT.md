# Creative Exploration Report: Morphology-Assisted Segmentation

## Goal

This creative exploration combines multiple computer vision techniques to improve segmentation. The existing morphology assignment demonstrated dilation, erosion, opening, closing, and connected components. This version adapts those ideas into a fuller segmentation pipeline.

## Method

The program uses:

- Otsu thresholding to automatically choose a grayscale cutoff.
- Sobel edge detection to capture boundaries that thresholding alone may miss.
- Mask combination to merge intensity regions with edge-supported structure.
- Morphological closing to connect broken boundaries and fill small gaps.
- Morphological opening to remove isolated noise.
- Connected-component labeling to color large segmented regions.
- Overlay blending to show the segmentation on top of the original image.

## Parameters

Compile:

```bash
gcc creative_morphology.c ../netpbm.c -o creative_morphology -lm
```

Commands used for the sample images:

```bash
./creative_morphology ../input_images/coastline.ppm 0.25 2 1 800
./creative_morphology ../input_images/tahoe.ppm 0.20 3 1 1200
```

Parameter meanings:

- Sobel edge threshold ratio: fraction of maximum gradient magnitude.
- Closing passes: number of dilation-then-erosion cleanup cycles.
- Opening passes: number of erosion-then-dilation noise-removal cycles.
- Minimum component size: small components below this pixel count are removed from the colored segmentation.

## Outputs

For each input image, outputs are saved in `output_images/`:

- `<input>_otsu_mask.pbm`
- `<input>_sobel_edge_mask.pbm`
- `<input>_combined_mask.pbm`
- `<input>_morph_cleaned.pbm`
- `<input>_morph_components.ppm`
- `<input>_morph_overlay.ppm`

## Observations

The combined approach is more useful than thresholding alone because Sobel edges preserve important boundaries around shoreline and terrain changes. Closing helps join fragmented structures into larger regions, while opening removes small speckles that would otherwise become many tiny components.

For coastline, a moderate edge threshold and two closing passes work well because the image has long, broken boundaries and texture variation. For Tahoe, a slightly lower edge threshold and more closing help connect lake and mountain boundary structures over a larger image area.

The limitation is that morphology depends strongly on scale. Too many closing passes can merge separate regions, while too much opening can erase thin structures. The saved intermediate masks make those tradeoffs visible.
