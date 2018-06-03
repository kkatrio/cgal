#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/internal/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Triangle_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;


int main()
{

  Triangle_mesh tm;
  std::ifstream input("data-coref/elephant.off");
  input >> tm;
  input.close();


  CGAL::Bbox_3 bbox = ::CGAL::Polygon_mesh_processing::bbox(tm);
  Triangle_mesh tm_out;


  K::Point_3 p1(0,0,0);
  K::Point_3 p2(0,-0.5,0);
  K::Point_3 p3(0.3,0,0);
  K::Plane_3 plane(p1, p2, p3);

  PMP::clip_to_bbox(plane, bbox, tm_out, params::all_default());

  std::ofstream output("data-coref/clipped_result.off");
  output << tm_out;
  output.close();





  return 0;
}
