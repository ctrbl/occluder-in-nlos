# Non-Line-of-Sight Imaging with Occluders (MATLAB Implementation)
This repository implements a non-line-of-sight (NLoS) imaging pipeline using MATLAB, focusing on estimating the occluder position and reconstructing the hidden scene using mathematical modeling and inverse problem-solving techniques.

## Overview 
The project explores computational reconstruction of an obscured object by analyzing indirect illumination on an imaging wall. Using camera measurements, we estimate:
1. The occluder's position
2. Hidden scene using light transport analysis

## Methodology
### 1. Estimating the Occluder Position
- Compute light transport matrix from camera measurements y.
- Use projection operator to refine the occluder's position by minimizng the Euclidean squared error.
- Background noise b is assumed negligible for reconstruction accuracy.

### 2. Reconstructing the Hidden Scene
- Use Moore-Penrose pseudoinverse to approximate the inverse of the light transport matrix.
- Compute pixel differences in neighboring 16x16 blocks for robustness against ambient light.
- Reconstruction methods:
  - Ideal conditions: Least-squares estimation for RGB reconstruction.
  - Real-world conditions: Total variation (TV) regularization to handle noise and model mismatch.

## Implementation Details
- Language: MATLAB

## Applications
- Seeing around corners (robotics, autonomous driving, etc.).
- Medical imaging, security, surveillance.
