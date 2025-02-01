# Non-Line-of-Sight Imaging with Occluders (MATLAB Implementation)
This repository implements a non-line-of-sight (NLoS) imaging pipeline using MATLAB, focusing on estimating the occluder position and reconstructing the hidden scene from a single indirect camera measurement. The setup simulates real-world conditions and uses inverse problem-solving and light transport analysis to estimate both the occluder's position and the hidden scene.

## Overview 
The project explores computational reconstruction of an obscured object by analyzing indirect illumination on an imaging wall. Using camera measurements, we estimate:
1. The occluder's position
2. Hidden scene using light transport analysis
![nlos](https://github.com/user-attachments/assets/6066e03b-1d02-4681-b5fb-f83a6fbfd094)

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
![image](https://github.com/user-attachments/assets/6d7a2704-879c-4a80-87f4-68918c44b641)


## Implementation Details
- Language: MATLAB
- Setup: A high-resolution camera captures indirect measurements, simulating real-world conditions with an occluder blocking direct line-of-sight.
- Light Transport: The light transport matrix is computed from camera measurements, accounting for occlusion effects.
- Occluder Position Estimation: The occluder position is refined by minimizing the Euclidean squared error between measured and predicted light transport.
- Scene Reconstruction: The Moore-Penrose pseudoinverse is used for scene reconstruction, with total variation (TV) regularization applied to handle noise and model mismatch.
- Noise Handling: Experiments include varying occluder positions and noise levels to simulate real-world conditions.

## Applications
- Seeing around corners (robotics, autonomous driving, etc.).
- Medical imaging, security, surveillance.

## Contributions
This project was a part of CIS 6930 - Computer Methods for Imaging taught by Professor John Murray-Bruce at the University of South Florida based on his previous works with non-light-of-sight imaging. 
- References:
  - S. W. Seidel, J. Murray-Bruce, Y. Ma, C. Yu, W. T. Freeman and V. K. Goyal, "Two-Dimensional Non-Line-of-Sight Scene Estimation From a Single Edge Occluder," in IEEE Transactions on Computational Imaging, vol. 7, pp. 58-72, 2021, doi: 10.1109/TCI.2020.3037405
