# T0-homog-strain-grad 


Topological optimization of homogenized strain-gradient materials


## Description 


This code computes the topological optimization of the unit-cell of a 2D periodic material, for which the criterion **s** can be chosen as a smooth function depending on the components of the Cauchy and strain-gradient homogenized tensors.
The rectangular unit-cell (sides **a**, **b**) is made of a stiff phase and a soft phase, the distribution of which is the design variable. 
The Young modulii of two phases differ from a contrast parameter **gamma**.


The topological optimization procedure is based on the algorithm proposed in:

*  _[Amstutz S., Andrä, H., A new algorithm for topology optimization using a
level-set method. J. Comput. Phys., 2006](https://www.sciencedirect.com/science/article/pii/S0021999105005656)_
  

For a presentation of the homogenization framework, the computation of the sensitivies and convention used in the code, and numerical results, we refer to:

* _[Calisti V., Synthèse de microstructures par optimisation topologique, et optimisation de forme d’un problème d’interaction fluide-structure, 2021](https://hal.univ-lorraine.fr/tel-03598154)_

* _[Calisti V., Lebée A., Novotny, A. A., Sokolowski J., Sensitivity of the second order homogenized elasticity tensor to topological microstructural changes, Journal of Elasticity, 2021](https://link.springer.com/article/10.1007/s10659-021-09836-6)_

* _[Calisti V., Lebée A., Novotny, A. A., Sokolowski J., Emergence of elastostatic strain-gradient effects from topological optimization, European Journal of Mechanics. A. Solids, 2023](https://link.springer.com/article/10.1007/s10659-021-09836-6)_



## Settings

* Definition: **Ch**= Cauchy homogenized tensor and **Dh**= strain-gradient homogenized tensor

* The file to configure to launch the optimization calculation is initialisation/initfile.m
  - Default: initial cell is a 1x1 square, and the initial distribution of stiff and soft is given by a centered disk **psi0**
  - Size of the rectangular cell can be chosen (**a**, **b**) 
  - Initial size of the mesh is given by **ni**
  - Young modulus (**E0**), material contrast (**gamma**), poisson ratio (**nu**)
  - The shape functioni (criterion) to be optimized:
      **s**= string expression of the functional, is written with the components
      of the selected tensors Chij, Dhij. 
      These components have to be surrounded with empty spaces. 
      Also, tensors need to be written Chij or Dhij with i <= j. 
      This is not restrictive because all tensors are symetric.
      _Example: **s** = ' Dh22 / Ch11  '_
  - The compliance tensor **Sh** can also be used


## Authors

* Valentin Calisti and Antonio André Novotny  _(TO of strain-gradient homogenization)_
* Samuel Amstutz and Antonio André Novotny _(TO of first gradient homogenization) with contribution of Sebastian Miguel Giusti_ 
  

