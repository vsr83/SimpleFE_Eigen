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

    std::map <int, int> partition;
    for (int ind_node=0; ind_node < mesh->num_nodes; ind_node++) {
        partition[ind_node] = 0;
    }
    for (int ind_line = 0; ind_line < mesh->num_lines; ind_line++) {
        int ind_elem = mesh->lines[ind_line];
        Mesh_Element *elem = mesh->elements[ind_elem];
        int node1 = elem->nodes[0];
        int node2 = elem->nodes[1];
        if (elem->physical == 206) {
            partition[node1-1] = 1;
            partition[node2-1] = 1;
        }
    }

    // Divide the nodes into free and bounary nodes.
    std::vector <int> nodes_bnd;
    std::vector <int> nodes_free;
    for (std::map<int, int>::iterator it = partition.begin(); it!= partition.end(); ++it) {
        int part = it->second;
        int node = it->first;
        if (part == 0) nodes_free.push_back(node);
        if (part == 1) nodes_bnd.push_back(node);
    }
    int num_free = nodes_free.size();
    int num_bnd  = nodes_bnd.size();

    // The partitions of the stiffness and mass matrices are computed by
    // multiplying the global matrices by from both sides with appropriate
    // sparse matrices, whichselect and sort the appropriate rows and columns.
    // As in assembly.cc, the sparse matrices are assembled using triplets,
    // which describe the rows, columns and values of non-zero elements.

    std::vector<Eigen::Triplet<double> > T_FF_L, T_FF_R,
                                         T_FI_L, T_FI_R,
                                         T_IF_L, T_IF_R,
                                         T_II_L, T_II_R;
    T_FF_L.reserve(num_free); T_FF_R.reserve(num_free);
    T_FI_L.reserve(num_bnd);  T_FI_R.reserve(num_bnd);
    T_IF_L.reserve(num_bnd);  T_IF_R.reserve(num_bnd);
    T_II_L.reserve(num_bnd);  T_II_R.reserve(num_bnd);

    for (int ind_free=0; ind_free<num_free; ind_free++) {
        Eigen::Triplet<double> triL(ind_free, nodes_free[ind_free], 1);
        Eigen::Triplet<double> triR(nodes_free[ind_free], ind_free, 1);
        T_FF_L.push_back(triL);
        T_FF_R.push_back(triR);
    }
    for (int ind_bnd=0; ind_bnd<num_bnd; ind_bnd++) {
        Eigen::Triplet<double> triL(ind_bnd, nodes_bnd[ind_bnd], 1);
        Eigen::Triplet<double> triR(nodes_bnd[ind_bnd], ind_bnd, 1);
        T_FI_L.push_back(triL);
        T_FI_R.push_back(triR);
    }
    Eigen::SparseMatrix<double> S_FF_L(num_free, mesh->num_nodes),
                                S_FF_R(mesh->num_nodes, num_free),
                                S_FI_L(num_free, mesh->num_nodes),
                                S_FI_R(mesh->num_nodes, num_bnd);
    S_FF_L.setFromTriplets(T_FF_L.begin(), T_FF_L.end());
    S_FF_R.setFromTriplets(T_FF_R.begin(), T_FF_R.end());
    S_FI_L.setFromTriplets(T_FI_L.begin(), T_FI_L.end());
    S_FI_R.setFromTriplets(T_FI_R.begin(), T_FI_R.end());

    // Compute the partitions of the stiffness matrix

    Eigen::SparseMatrix<double> S_FF(num_free, num_free),
                                S_FI(num_free, num_bnd),
                                F_F(num_free, 1);

    S_FF = S_FF_L * *ass.global_stiff * S_FF_R;
    S_FI = S_FI_L * *ass.global_stiff * S_FI_R;
    F_F  = S_FF_L * ass.excitation;

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

/*    for (int ind=0; ind < num_free; ind++) {
        int node = nodes_free[ind];
        std::cout << mesh->nodes[node*3] << " " << mesh->nodes[node*3+1]
                                         << " " << sol_dense(ind, 0)
                                         << std::endl;
    }
*/
    for (int ind_triangle=0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        Mesh_Element *elem = mesh->elements[ind_triangle];
        for (int ind_node=0; ind_node < 3; ind_node++) {
            int node = elem->nodes[ind_node]-1;

            std::vector<int>::iterator p_free;
            p_free = std::find(nodes_free.begin(), nodes_free.end(), node);

            if (p_free != nodes_free.end()) {
                std::cout << mesh->nodes[node*3]   << " "
                          << mesh->nodes[node*3+1] << " "
                          << sol_dense(p_free - nodes_free.begin(), 0) << std::endl;
            }
        }
    }


    delete mesh;
    delete meshfile;
}
