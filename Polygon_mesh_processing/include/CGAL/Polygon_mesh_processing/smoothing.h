#ifndef SMOOTHING_H
#define SMOOTHING_H

#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>
#include <Eigen/Dense>

namespace CGAL {

namespace Polygon_mesh_processing {



typedef CGAL::Eigen_solver_traits<> Eigen_solver;

using namespace Eigen;
using namespace std;

void solve_linear_system()
{


    std::size_t degree = 2;
    std::size_t nb_nonzero_coef = 3;

    Eigen_solver::Matrix A(degree);
    Eigen_solver::Vector B(degree);


    A.add_coef(0, 0, 2);
    A.add_coef(1, 1, 3);
    A.add_coef(0, 1, 1);

    B(0) = 1;
    B(1) = 1;

    A.assemble_matrix();

    Eigen_solver::Vector X(degree);


    Eigen_solver solver;

    Eigen_solver::NT d;

    if (!(solver.linear_solver (A, B, X, d)))
      {
        std::cerr << "Error: linear solver failed" << std::endl;
      }
    std::cerr << "Linear solve succeeded" << std::endl;


    std::cout<<"NT= "<<d<<std::endl;

    std::cout<<"X= "<<X<<std::endl;




}





} // PMP
} // CGAL


#endif // SMOOTHING_H
