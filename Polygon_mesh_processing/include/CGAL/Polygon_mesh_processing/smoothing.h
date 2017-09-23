#ifndef SMOOTHING_H
#define SMOOTHING_H

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>
#include <Eigen/Sparse>

namespace CGAL {

namespace Polygon_mesh_processing {



template<typename PolygonMesh, typename VertexPointMap,
         typename CotangentValue = CGAL::internal::Cotangent_value_Meyer<PolygonMesh, VertexPointMap>>
struct Cotangent_weight : CotangentValue
{
    Cotangent_weight(PolygonMesh& pmesh_, VertexPointMap vpmap_)
      : CotangentValue(pmesh_, vpmap_)
    {}

    PolygonMesh& pmesh()
    {
      return CotangentValue::pmesh();
    }

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;

    double operator()(halfedge_descriptor he)
    {
      halfedge_descriptor h1 = next(he, pmesh());
      halfedge_descriptor h2 = prev(opposite(he, pmesh()), pmesh());

      vertex_descriptor vs = source(he, pmesh());
      vertex_descriptor vt = target(he, pmesh());
      vertex_descriptor v1 = target(h1, pmesh());
      vertex_descriptor v2 = source(h2, pmesh());

      return ( CotangentValue::operator()(vs, v1, vt) + CotangentValue::operator()(vs, v2, vt) ) / 2.0;
    }
};








template<typename PolygonMesh, typename VertexPointMap>
class Shape_smoother{

private:


    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    typedef CGAL::Eigen_solver_traits< Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > > Solver_traits;


    // vertex index map
    typedef typename boost::property_map<PolygonMesh, boost::vertex_index_t>::type IndexMap;
    IndexMap idxmap = get(boost::vertex_index, mesh_);





public:


    Shape_smoother(PolygonMesh& mesh, VertexPointMap& vpmap) : mesh_(mesh), vpmap_(vpmap), weight_calculator_(mesh, vpmap)
    {
        nb_vertices = vertices(mesh_).size();
    }




    void calculate_coeff_matrix()
    {


        // linear system data
        typename Solver_traits::Matrix A(nb_vertices, nb_vertices);



        for(vertex_descriptor vi : vertices(mesh_))
        {
            if(!is_border(vi, mesh_))
            {
                double sum_L = 0;
                for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
                {
                    // calculate weight
                    double Lij = weight_calculator_(h);
                    sum_L += Lij;

                    vertex_descriptor vj = target(h, mesh_);

                    A.add_coef(idxmap[vi], idxmap[vj], Lij);
                }

                A.add_coef(idxmap[vi], idxmap[vi], sum_L);

            }


        }



    }






    void solve_linear_system()
    {

        typedef CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> Eigen_solver;


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



private:

    // geometry data
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    Cotangent_weight<PolygonMesh, VertexPointMap> weight_calculator_;
    Solver_traits solver_;
    std::size_t nb_vertices;





};








template<typename PolygonMesh>
void smooth_shape(PolygonMesh& mesh)
{

    // VPmap type
    typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
    VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

    CGAL::Polygon_mesh_processing::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);
    smoother.calculate_coeff_matrix();


}











} // PMP
} // CGAL


#endif // SMOOTHING_H
