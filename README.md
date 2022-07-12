# Spatial symmetries in multipolar metasurfaces

Computes the scattering matrix for oblique wave propagation and the material parameter tensors (dipolar and quadrupolar responses) that correspond to the spatial symmetries of a metasurface made of a periodic square lattice.

## Requirements

Python >= 3.9

Python packages: tkinter, numpy, sympy, scipy, matplotlib


## Example

Consider the following metasurface unit cell composed of a flat T-shape structure 

<img src="/images/T-struct.png" width="250">

## List of implemented symmetries

```
1. px: symmetry along the x-axis
2. py: symmetry along the y-axis
3. pz: symmetry along the z-axis
4. i: inversion symmetry
5. c2x: 180° symmetry around the x-axis
6. c2y: 180° symmetry around the y-axis
7. c2z: 180° symmetry around the z-axis
8. c4z: 90° symmetry around the z-axis
```

## References
[1] Bernal Arango, Felipe, Toon Coenen, and A. Femius Koenderink. "Underpinning hybridization intuition for complex nanoantennas by magnetoelectric quadrupolar polarizability retrieval." ACS Photonics 1.5 (2014): 444-453.

[2] Dmitriev, Victor. "Symmetry properties of electromagnetic planar arrays in transfer matrix description." IEEE transactions on antennas and propagation 61.1 (2012): 185-194.
