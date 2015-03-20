/* A Simple FE solver.
   Copyright (C) 2015 Ville Räisänen <vsr at vsr.name>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "assembly.h"
#include "partition.h"
#include <Eigen/SparseLU>
#include <map>

int
main(int argc, char **argv) {
    double mu0 = 4*M_PI*1e-7;

    // Region data
    Region region_Fe ("Armature", 1/(1000*mu0), 0, 0);
    Region region_Cu ("Coil", 1/(mu0), 0, 1e7);
    Region region_Air("Air", 1/(mu0), 0, 0);
    std::map <int, Region> regions;

    regions[201] = region_Fe;
    regions[202] = region_Fe;
    regions[203] = region_Cu;
    regions[204] = region_Air;
    regions[205] = region_Air;
    regions[206] = region_Air;

    MeshFile *meshfile = new MeshFile("valve.msh");
    Mesh *mesh = new Mesh(meshfile);

    // Assemble the global stiffness and mass matrices.
    Assembly ass(mesh, 6, regions);

    // Partition nodes into free and bounded nodes (see partition.h).
    std::map <int, int> phys_partition;
    phys_partition[206]=1;
    Partition part(mesh, phys_partition);
    int num_free = part.parts_size[0];
    int num_bnd  = part.parts_size[1];
    Eigen::SparseMatrix<double> S_FI(num_free, num_bnd),
                                S_FF(num_free, num_free),
                                F_F(num_free, 1);
    S_FF = part.part_left(0, 0)*(*ass.global_stiff)*part.part_right(0, 0);
    S_FI = part.part_left(0, 1)*(*ass.global_stiff)*part.part_right(0, 1);
    F_F  = part.part_left(0, 0) * ass.excitation;

    // Solve the system of equations.
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    Eigen::SparseMatrix<double> sol_sparse;
    Eigen::MatrixXd sol_dense;
    solver.compute(S_FF);
    if (solver.info() != Eigen::Success) {
        std::cerr << "SparseLU decomposition failed." << std::endl;
        exit(-1);
    }

    sol_sparse = solver.solve(F_F);
    sol_dense = Eigen::MatrixXd(sol_sparse);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solution failed." << std::endl;
        exit(-1);
    }

    // Export the DoF's and the coordinates of the corresponding nodes.

    for (int ind_triangle=0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        Mesh_Element *elem = mesh->elements[ind_triangle];
        for (int ind_node=0; ind_node < 3; ind_node++) {
            int node = elem->nodes[ind_node]-1;

            std::vector<int>::iterator p_free;
            p_free = std::find(part.node_parts[0].begin(),
                               part.node_parts[0].end(), node);

            if (p_free != part.node_parts[0].end()) {
                std::cout << mesh->nodes[node*3]   << " "
                          << mesh->nodes[node*3+1] << " "
                          << sol_dense(p_free - part.node_parts[0].begin(), 0)
                          << std::endl;
            }
        }
    }

    delete mesh;
    delete meshfile;
}
