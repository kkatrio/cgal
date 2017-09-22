#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/smoothing.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;



int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/shape.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }





  CGAL::Polygon_mesh_processing::solve_linear_system();


  std::cout<<"done"<<std::endl;


  return 0;
}
