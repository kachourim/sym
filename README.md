# Spatial symmetries in multipolar metasurfaces

## Requirements

Python >= 3.9

Python packages: tkinter, numpy, sympy, scipy, matplotlib

## Description

Based on [1,2], the provided scripts computes the scattering matrix for oblique wave propagation and the material parameter tensors (dipolar and quadrupolar responses) that correspond to the spatial symmetries of a metasurface made of a periodic square lattice.

### Example

Consider the following reciprocal metasurface unit cell composed of a flat T-shaped structure 

<img src="/images/T-struct.png" width="250">

The structure is aligned with the metasurface lattice and its symmetries (Px, Pz, C2y) are emphasized (descriptions given in the next section).

To use the provided script, enter the symmetries of the structure in the first box, as follows

<img src="/images/interface.png" width="500">

Then, with your cursor placed in the first box, press Enter to obtain the symmetry-allowed material parameters:

<img src="/images/MaterialTensors.png" width="500">

Note that the numbers and colors are arbitrary and are only useful to see which of the components are equal to each other.

The scattering matrix can then be obtained by pressed Enter in the second box. By default, the plane of incidence is defined for $\phi = 0$ (xz-plane)

<img src="/images/Smatrix.png" width="400">

For the same metasurface, the scattering matrix is, for $\phi = 90^\circ$, given by 

<img src="/images/Smatrix2.png" width="400">

## List of implemented symmetries

```
1. px: symmetry along the x-axis
2. py: symmetry along the y-axis
3. pz: symmetry along the z-axis
4. i: inversion symmetry
5. c2x: 180° rotation symmetry around the x-axis
6. c2y: 180° rotation symmetry around the y-axis
7. c2z: 180° rotation symmetry around the z-axis
8. c4z: 90° rotation symmetry around the z-axis
```

## References
[1] Bernal Arango, Felipe, Toon Coenen, and A. Femius Koenderink. "Underpinning hybridization intuition for complex nanoantennas by magnetoelectric quadrupolar polarizability retrieval." ACS Photonics 1.5 (2014): 444-453.

[2] Dmitriev, Victor. "Symmetry properties of electromagnetic planar arrays in transfer matrix description." IEEE transactions on antennas and propagation 61.1 (2012): 185-194.
