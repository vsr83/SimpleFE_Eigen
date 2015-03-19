This is a quick implementation of a very simple FE solver using Eigen.

Mesh files generated with GMsh are imported with mesh.cc, mesh_element.cc,
mesh_file.cc and mesh.cc. Material parameters are specified with Region-
objects and assembled into a map according to the "physical numbers" 
(see region.cc and region.h). The element stiffness and mass matrices for 
first order basis functions are computed with Gaussian quadrature and 
assembled into global stiffness and mass matrices in element.cc and 
assembly.cc. 

Sparse LU solver of Eigen is used for the solution of the resulting system
of equations. This code is applicable to large problems even with ~1M 
degrees of freedom.

The code is applied to a simple magnetostatic problem in example_valve.cc.