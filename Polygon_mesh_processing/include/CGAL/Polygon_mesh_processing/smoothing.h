#ifndef SMOOTHING_H
#define SMOOTHING_H

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>
#include <Eigen/Sparse>

#include <fstream>
#include <unordered_set>

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

template<typename PolygonMesh>
struct Incident_areas
{

  Incident_areas(PolygonMesh& mesh) : pmesh(mesh){}

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  double operator()(halfedge_descriptor he)
  {
    halfedge_descriptor hopp = opposite(he, pmesh);
    face_descriptor f1 = face(he, pmesh);
    face_descriptor f2 = face(hopp, pmesh);
    double A1 = face_area(f1, pmesh);
    double A2 = face_area(f2, pmesh);
    return A1 + A2;
  }

  // data
  PolygonMesh& pmesh;

};









template<typename PolygonMesh, typename VertexPointMap>
class Shape_smoother{

// types
private:

    typedef typename GetGeomTraits<PolygonMesh>::type GeomTraits;
    typedef typename GeomTraits::FT NT;
    typedef typename GeomTraits::Point_3 Point;

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    // vertex index map
    typedef typename boost::property_map<PolygonMesh, boost::vertex_index_t>::type IndexMap;

    // solver
    typedef typename CGAL::Eigen_sparse_matrix<double>::EigenType EigenMatrix; 

    //typedef CGAL::Eigen_solver_traits< Eigen::SimplicialLDLT< EigenMatrix > > Solver_traits;
    typedef CGAL::Eigen_solver_traits<> Solver_traits; //BicGSTAB

    typedef typename Solver_traits::Matrix Matrix;
    typedef typename Solver_traits::Vector Vector;



// data
private:

    IndexMap vimap_ = get(boost::vertex_index, mesh_);

    // geometry data
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    Cotangent_weight<PolygonMesh, VertexPointMap> weight_calculator_;
    Incident_areas<PolygonMesh> inc_areas_calculator_;

    std::size_t nb_vert_;

    // linear solver
    Solver_traits solver_;

    // constrained vertices
    std::unordered_set<vertex_descriptor> constrained_vertices_;


// operations
private:

    Matrix compute_stiffness_matrix()
    {

      Matrix L(nb_vert_, nb_vert_);
      for(vertex_descriptor vi : vertices(mesh_))
      {
        if(!is_border(vi, mesh_))
        {
          NT sum_Lik = 0;
          for(halfedge_descriptor h : halfedges_around_source(vi, mesh_))
          {
            // calculate L
            NT Lij = weight_calculator_(h);
            sum_Lik -= Lij;

            vertex_descriptor vj = target(h, mesh_);

            // A = (I - L) -> 0 - Lij
            //L.set_coef(vimap_[vi], vimap_[vj], -Lij, true);
            L.set_coef(vimap_[vi], vimap_[vj], Lij, true);
          }

          // diagonal, A = (I - L) -> 1 - Lii
          //L.set_coef(vimap_[vi], vimap_[vi], 1.0 - sum_Lik, true);
          L.set_coef(vimap_[vi], vimap_[vi], sum_Lik, true);
        }
      }

      return L;
    }


    Matrix compute_mass_matrix()
    {
      Matrix D(nb_vert_, nb_vert_);
      for(vertex_descriptor vi : vertices(mesh_))
      {
        if(!is_border(vi, mesh_))
        {

          NT sum_Dik = 0;
          for(halfedge_descriptor h : halfedges_around_source(vi, mesh_))
          {
            NT Tij = inc_areas_calculator_(h) / 12.0;
            sum_Dik += Tij;

            vertex_descriptor vj = target(h, mesh_);
            D.set_coef(vimap_[vi], vimap_[vj], Tij, true);
          }

          D.set_coef(vimap_[vi], vimap_[vi], sum_Dik, true);

        }
      }

      return D;
    }

    void compute_coeff_matrix(Matrix& A)
    {

      Matrix L = compute_stiffness_matrix();
      Matrix D = compute_mass_matrix();

      CGAL_assertion(L.row_dimension() == L.column_dimension());
      CGAL_assertion(D.row_dimension() == D.column_dimension());
      CGAL_assertion(L.row_dimension() == D.row_dimension() &&
                     L.column_dimension() == D.column_dimension());
      CGAL_assertion(L.row_dimension() == A.row_dimension() &&
                     L.column_dimension() == A.column_dimension());

      for(std::size_t i=0; i<A.row_dimension(); ++i)
      {
        for(std::size_t j=0; j<A.column_dimension(); ++j)
        {
          A.set_coef(i, j, D.get_coef(i,j) - L.get_coef(i,j), true);
        }
      }

      set_constraints(A);
      A.assemble_matrix();
    }


    void compute_rhs(Vector& Bx, Vector& By, Vector& Bz)
    {

        for(vertex_descriptor vi : vertices(mesh_))
        {
            int index = vimap_[vi];
            Point p = get(vpmap_, vi);
            Bx.set(index, p.x());
            By.set(index, p.y());
            Bz.set(index, p.z());
        }


        //temp solution: calc D again
        Matrix D = compute_mass_matrix();




    }

    // to be called after gather_constrained_vertices()
    void set_constraints(Matrix& A)
    {
        typename std::unordered_set<vertex_descriptor>::iterator it;
        for(it = constrained_vertices_.begin(); it != constrained_vertices_.end(); ++it)
        {
            int i = get(vimap_, *it);
            A.set_coef(i, i, 1.0, true);
        }
    }

    void extract_matrix(Matrix& A)
    {
        std::ofstream out("data/matA.dat");
        for(int j=0; j < A.column_dimension(); ++j)
        {
            for(int i=0; i < A.row_dimension(); ++i)
            {
                NT val = A.get_coef(i, j);
                out<<val<<" ";
            }
            out<<std::endl;
        }
        out.close();
    }

    void extract_solution(Vector& Xx, Vector& Xy, Vector& Xz)
    {
        std::ofstream out("data/solution.dat");
        for(int j=0; j < Xx.dimension(); ++j)
        {
            out<<Xx[j]<<"\t"<<Xy[j]<<"\t"<<Xx[j]<<"\t"<<std::endl;
        }
        out.close();
    }


    void update_mesh(Vector& Xx, Vector& Xy, Vector& Xz)
    {

        for (vertex_descriptor v : vertices(mesh_))
        {
            int index = get(vimap_, v);
            NT x_new = Xx[index];
            NT y_new = Xy[index];
            NT z_new = Xz[index];
            put(vpmap_, v, Point(x_new, y_new, z_new));
        }
    }

    void gather_constrained_vertices()
    {
        for(vertex_descriptor v : vertices(mesh_))
        {
            if(is_border(v, mesh_))
            {
                constrained_vertices_.insert(v);
            }
        }
    }

    // printing matrix
    void print(Matrix& Mat)
    {
      for(int i=0; i<Mat.row_dimension(); ++i)
      {
        for(int j=0; j<Mat.column_dimension(); ++j)
          std::cout<<Mat.get_coef(i,j)<<" ";
        std::cout<<std::endl;
      }
    }

public:
    // matrix multiplication
    Matrix multiply(Matrix& A, Matrix& B)
    {
      CGAL_assertion(A.row_dimension()>0 && A.column_dimension()>0);
      CGAL_assertion(B.row_dimension()>0 && B.column_dimension()>0);
      CGAL_assertion(A.row_dimension() == B.column_dimension());
      CGAL_assertion(B.row_dimension() == A.column_dimension());

      std::size_t degree = A.row_dimension();
      std::size_t sz = A.column_dimension();

      Matrix C(degree, degree);
      for(std::size_t i=0; i<degree; ++i)
      {
        for(std::size_t j=0; j<degree; ++j)
        {
          double v =0;
          for(std::size_t k=0; k<sz; ++k)
          {
            v += A.get_coef(i,k) * B.get_coef(k,i);
          }
          C.set_coef(i,j,v,true);
        }
      }

      return C;
    }







public:

    Shape_smoother(PolygonMesh& mesh, VertexPointMap& vpmap) : mesh_(mesh), vpmap_(vpmap),
        weight_calculator_(mesh, vpmap),
        inc_areas_calculator_(mesh),
        nb_vert_(static_cast<int>(vertices(mesh).size()))
    { }


    void solve_system()
    {

        gather_constrained_vertices();

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
            bool solved = false;
            CGAL_assertion(solved);
        }

        //extract_solution(Xx, Xy, Xz);

        update_mesh(Xx, Xy, Xz);

    }



    void test_matrix_product()
    {

      Matrix A(3,3);

   /*   double k = 0;
      for(int i=0; i<3; ++i)
      {
        for(int j=0; j<3; ++j)
        {
          A.set_coef(i,j,k,true);
          k++;
        }
      }

      Matrix Bh(3,1);
      for(int i=0; i<3; ++i)
      {
        Bh.set_coef(i,1, i*i, true);
      }

      Matrix Bg(1,3);
      for(int i=0; i<3; ++i)
      {
        Bg.set_coef(1,i, i*i, true);
      }
*/

      Matrix X = multiply(A, A);
      print(X);






    }






};








template<typename PolygonMesh>
void smooth_shape(PolygonMesh& mesh)
{

    // VPmap type
    typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
    VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

    CGAL::Polygon_mesh_processing::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);

    //smoother.solve_system();
    smoother.test_matrix_product();


}











} // PMP
} // CGAL


#endif // SMOOTHING_H
