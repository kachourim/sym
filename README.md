# Spatial symmetries in multipolar metasurfaces

## Requirements

Python >= 3.9

Python packages: tkinter, numpy, sympy, scipy, matplotlib

## Description

Based on [1,2], the provided script computes the scattering matrix for oblique wave propagation and the material parameter tensors (dipolar and quadrupolar responses) that correspond to the spatial symmetries of a metasurface made of a periodic square lattice.

The script works by specifying the list of symmetries that the metasurface unit cell possesses, those that are available are defined below

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

The script then computes the components of the material parameter tensors that are allowed to exist due to the spatial symmetries of the system. These components are relating the electric dipole ($p$), the magnetic dipole ($m$), the electric quadrupole ($Q$) and the magnetic quadrupole ($S$) to the electric field ($E$), the magnetic field ($H$), the gradient of the electric field ($\diamond E$) and the gradient of the magnetic field ($\diamond H$) such that [3]

<img src="/images/PMQS.png" width="400">

where 

$\diamond E = \[(\partial_x E_y+\partial_y E_x)/2,(\partial_x E_z+\partial_z E_x)/2 (\partial_y E_z+\partial_z E_y)/2, \partial_xE_x, \partial_yE_y, \partial_zE_z  \]$

Using reciprocity [3] and considering that the quadrupolar tensors are irreducible, we obtain a material parameter matrix similar to those shown in [1].

The scattering matrix is obtained following the procudure described in [2]. It connects the input ports (1, 2, 3, 4) to the output ports (1, 2, 3, 4), where even numbers correspond to TE polarization whereas odd numbers correspond to TM polarization, as such

<img src="/images/S-description.png" width="800">


### Example

Consider the following reciprocal metasurface unit cell composed of a flat T-shaped structure 

<img src="/images/T-struct.png" width="250">

The structure is aligned with the metasurface lattice and its symmetries (Px, Pz, C2y) are emphasized.

To use the provided script, enter the symmetries of the structure in the first box, as follows

<img src="/images/interface.png" width="500">

Then, with your cursor placed in the first box, press Enter to obtain the symmetry-allowed material parameters:

<img src="/images/MaterialTensors.png" width="500">

Note that the numbers and colors are arbitrary and are only useful to see which of the components are equal to each other.

The scattering matrix [2] can be obtained by pressing Enter in the second box. By default, the plane of incidence is defined for $\phi = 0$ (xz-plane)

<img src="/images/Smatrix.png" width="400">

For the same metasurface, the scattering matrix is, for $\phi = 90^\circ$, given by 

<img src="/images/Smatrix2.png" width="400">




## References
[1] Bernal Arango, Felipe, Toon Coenen, and A. Femius Koenderink. "Underpinning hybridization intuition for complex nanoantennas by magnetoelectric quadrupolar polarizability retrieval." ACS Photonics 1.5 (2014): 444-453.

[2] Dmitriev, Victor. "Symmetry properties of electromagnetic planar arrays in transfer matrix description." IEEE transactions on antennas and propagation 61.1 (2012): 185-194.

[3] Achouri, Karim, and Olivier JF Martin. "Extension of Lorentz reciprocity and Poynting theorems for spatially dispersive media with quadrupolar responses." Physical Review B 104.16 (2021): 165426.
