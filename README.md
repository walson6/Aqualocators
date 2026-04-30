# Aqualocators Computer Vision Project

## Overview
This project applies computer vision techniques to real-world image datasets. The pipeline includes image filtering, edge detection, segmentation, and robustness analysis.

Each algorithm is implemented independently in its own folder with corresponding source code and outputs.

## Authors

- Arthur Wong  
- Jimmy Vu

---

## Project Structure

Each algorithm is stored in its own folder:

```text
Aqualocators/
├── Gaussian_Filter/
│   ├── output_images/
│   ├── gaussian.c
├── Sobel_Edge/
│   ├── output_images/
│   ├── sobel.c
├── Canny_Edge/
│   ├── output_images/
│   ├── canny.c
├── input_images/
│   ├── input1.jpg
│   ├── input1.pgm
│   ├── input1.ppm
└── README.md
```

---

## Preprocessing

The project uses the Netpbm image formats (PPM/PGM). Therefore, input images in formats such as `.jpg` or `.png` should be converted before processing.  
  
This can be done using ImageMagick:  

```bash
magick input.jpg output.ppm   # for color images
magick input.jpg output.pgm   # for grayscale images
```

---

## Compilation

Ensure you are inside the Algorithm_Name_Module/ directory before compiling

Each module compiles independently:
```bash
gcc code.c netpbm.c -o program -lm

# example compilation for gaussian.c
gcc gaussian.c ../netpbm.c -o gaussian -lm
```

The instructions for compiling each algorithm are included in comments within their respective .c files

---

## Example Execution

### Gaussian Filter
```bash
cd Gaussian_Filter gcc gaussian.c ../netpbm.c -o gaussian -lm  
gaussian.exe ../input_images/coastline.ppm
```

The instructions for running each algorithm are included in comments within their respective .c files

---

## Outputs

Each module generates:
- output_images/ → processed images
