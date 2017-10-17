#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <fstream>

#include <CGAL/Polygon_mesh_processing/smoothing.h>

#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;



double measure_radius(Mesh& mesh)
{

  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;


  //choose randomly 4 points from the mesh
  CGAL::Random random;
  std::size_t numberOfPoints = mesh.number_of_vertices();


  std::vector<vertex_descriptor> verts;
  std::copy(vertices(mesh).begin(), vertices(mesh).end(), std::back_inserter(verts));

  K::Point_3 p[4];
  for(unsigned int i=0; i<4; ++i)
  {
    vertex_descriptor v = verts[random.get_int(0,numberOfPoints)];
    p[i] = get(vpmap, v);
  }

  K::Sphere_3 sphere(p[0], p[1], p[2], p[3]);

  return CGAL::sqrt(sphere.squared_radius());

}






int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/sphere.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }


  CGAL::Polygon_mesh_processing::smooth_shape(mesh);



  std::ofstream out("data/sphere_smoothed.off");
  out<<mesh;
  out.close();



  std::cout<<"done"<<std::endl;


  return 0;
}
