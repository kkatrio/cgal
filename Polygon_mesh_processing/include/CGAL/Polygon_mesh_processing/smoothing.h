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

#include <fstream>

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

// types
private:

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    typedef CGAL::Eigen_sparse_matrix<double>::EigenType EigenMatrix;
   // typedef CGAL::Eigen_solver_traits< Eigen::SimplicialLDLT< EigenMatrix > > Solver_traits;
     typedef CGAL::Eigen_solver_traits<> Solver_traits;

    typedef typename Solver_traits::Matrix Matrix;
    typedef typename Solver_traits::Vector Vector;


    // vertex index map
    typedef typename boost::property_map<PolygonMesh, boost::vertex_index_t>::type IndexMap;

    typedef typename GetGeomTraits<PolygonMesh>::type GeomTraits;

    typedef typename GeomTraits::Point_3 Point;

    typedef typename GeomTraits::FT NT;


// data
private:

    IndexMap vimap = get(boost::vertex_index, mesh_);

    // geometry data
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    Cotangent_weight<PolygonMesh, VertexPointMap> weight_calculator_;

    std::size_t nb_vert_;

    // linear solver
    Solver_traits solver_;


// operations
private:

    void compute_coeff_matrix(Matrix& A)
    {

        for(vertex_descriptor vi : vertices(mesh_))
        {
            //if(!is_border(vi, mesh_))
            //{
                double sum_Lik = 0;
                for(halfedge_descriptor h : halfedges_around_source(vi, mesh_))
                {
                    // calculate weight
                    double Lij = weight_calculator_(h);
                    sum_Lik -= Lij;

                    vertex_descriptor vj = target(h, mesh_);

                    A.add_coef(vimap[vi], vimap[vj], -Lij);
                }

                // diagonal
                A.add_coef(vimap[vi], vimap[vi], 1.0 - sum_Lik);
           // }
        }


       // A.assemble_matrix(); // remove with set
    }


    void compute_rhs(Vector& Bx, Vector& By, Vector& Bz)
    {

        for(vertex_descriptor vi : vertices(mesh_))
        {
            //if(!is_border(vi, mesh_)) // border ones do not move
           // {
                int index = vimap[vi];
                Point p = get(vpmap_, vi);
                Bx.set(index, p.x());
                By.set(index, p.y());
                Bz.set(index, p.z());
          //  }

        }


    }


/*
    void extract_matrix(Matrix& A)
    {


        std::ofstream out("matA.dat");
        for(int j=0; j < A.column_dimension(); ++j)
        {
            for(int i=0; i < A.row_dimension(); ++i)
            {
                NT val = A.
                out<<val;
            }
            out<<endl;
        }

        out.close();
    }
*/







public:

    Shape_smoother(PolygonMesh& mesh, VertexPointMap& vpmap) : mesh_(mesh), vpmap_(vpmap),
        weight_calculator_(mesh, vpmap),
        nb_vert_(static_cast<int>(vertices(mesh).size()))
    { }


    void solve_system()
    {
        Matrix A(nb_vert_, nb_vert_);
        Vector Bx(nb_vert_);
        Vector By(nb_vert_);
        Vector Bz(nb_vert_);
        Vector Xx(nb_vert_);
        Vector Xy(nb_vert_);
        Vector Xz(nb_vert_);

        compute_coeff_matrix(A);
        compute_rhs(Bx, By, Bz);

        //extract_matrix(A);

        NT Dx, Dy, Dz;
        if(!solver_.linear_solver(A, Bx, Xx, Dx) ||
           !solver_.linear_solver(A, By, Xy, Dy) ||
           !solver_.linear_solver(A, Bz, Xz, Dz) )
        {
            std::cerr<<"Could not solve linear system!"<<std::endl;
        }





        std::cout<<"Xx= "<<Xx<<std::endl;
        std::cout<<"Xy= "<<Xy<<std::endl;
        std::cout<<"Xz= "<<Xz<<std::endl;
    }











};








template<typename PolygonMesh>
void smooth_shape(PolygonMesh& mesh)
{

    // VPmap type
    typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
    VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

    CGAL::Polygon_mesh_processing::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);

    smoother.solve_system();


}











} // PMP
} // CGAL


#endif // SMOOTHING_H
