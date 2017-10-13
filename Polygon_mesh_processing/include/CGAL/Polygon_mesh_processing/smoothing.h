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

    CGAL_assertion(A1>0 && A2>0);

    return A1 + A2;
  }

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



  //CGAL Sparse solver
  typedef typename CGAL::Eigen_sparse_matrix<double>::EigenType CEigenMatrix;

  typedef CGAL::Eigen_solver_traits< Eigen::SimplicialLDLT< CEigenMatrix > > Solver_traits;
  //typedef CGAL::Eigen_solver_traits<> Solver_traits; //BicGSTAB

  typedef typename Solver_traits::Matrix Matrix;
  typedef typename Solver_traits::Vector Vector;



  // Eigen sparse matrix & vector
  typedef typename Eigen::SparseMatrix<double> Eigen_matrix;
  typedef typename Eigen::SparseVector<double> Eigen_vector;



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


    Eigen_matrix get_stiffness_matrix()
    {
      Eigen_matrix mat(nb_vert_, nb_vert_);
      for(vertex_descriptor vi : vertices(mesh_))
      {
        if(!is_border(vi, mesh_))
        {
          NT sum_Lik = 0;
          for(halfedge_descriptor h : halfedges_around_source(vi, mesh_))
          {
            NT Lij = weight_calculator_(h);
            sum_Lik -= Lij;
            vertex_descriptor vj = target(h, mesh_);
            mat.coeffRef(vimap_[vi], vimap_[vj]) = Lij;
          }

          mat.coeffRef(vimap_[vi], vimap_[vi]) = sum_Lik;
        }
      }

      return mat;
    }

    Eigen_matrix get_mass_matrix()
    {
      Eigen_matrix mat(nb_vert_, nb_vert_);
      for(vertex_descriptor vi : vertices(mesh_))
      {
        if(!is_border(vi, mesh_))
        {
          NT sum_Dik = 0;
          for(halfedge_descriptor h : halfedges_around_source(vi, mesh_))
          {
            //vertex_descriptor vj = target(h, mesh_);
            NT Dij = inc_areas_calculator_(h) / 12.0;
            sum_Dik += Dij;
            vertex_descriptor vj = target(h, mesh_);
            mat.coeffRef(vimap_[vi], vimap_[vj]) = Dij;
          }

          mat.coeffRef(vimap_[vi], vimap_[vi]) = sum_Dik;
        }
      }

      return mat;
    }


    void compute_coeff_matrix(Matrix& A)
    {
      Eigen_matrix L = get_stiffness_matrix();
      Eigen_matrix D = get_mass_matrix();


      double delta = 0.01;


      //Eigen_matrix D(L.rows(), L.cols());
      //D.setIdentity();


      Eigen_matrix Ae = D - delta * L;


      fill_sparse_matrix(A, Ae);

    }


    void compute_rhs(Vector& Bx, Vector& By, Vector& Bz)
    {

      Eigen_vector Xx(nb_vert_);
      Eigen_vector Xy(nb_vert_);
      Eigen_vector Xz(nb_vert_);

      for(vertex_descriptor vi : vertices(mesh_))
      {
          int index = vimap_[vi];
          Point p = get(vpmap_, vi);
          Xx.coeffRef(index) = p.x();
          Xy.coeffRef(index) = p.y();
          Xz.coeffRef(index) = p.z();
      }


      Eigen_matrix D = get_mass_matrix();


      //Eigen_matrix D(Bx.rows(), Bx.rows());
      //D.setIdentity();




      Eigen_vector Bxe = D * Xx;
      Eigen_vector Bye = D * Xy;
      Eigen_vector Bze = D * Xz;


      fill_sparse_vector(Bx, Bxe);
      fill_sparse_vector(By, Bye);
      fill_sparse_vector(Bz, Bze);

    }


    void fill_sparse_matrix(Matrix&A, Eigen_matrix& Ae)
    {
      CGAL_assertion(A.row_dimension() == Ae.rows() &&
                     A.column_dimension() == Ae.cols());
      for(std::size_t i=0; i<Ae.rows(); ++i)
      {
        for(std::size_t j=0; j<Ae.cols(); ++j)
        {
          A.set_coef(i, j, Ae.coeff(i, j), true);
        }
      }
    }

    void fill_sparse_vector(Vector& V, Eigen_vector& Ve)
    {
      CGAL_assertion(V.rows() == Ve.rows());
      for(std::size_t i= 0; i<Ve.rows(); ++i)
      {
        V[i] = Ve.coeff(i);
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

    void apply_constraints(Matrix& A)
    { // to be called after gather_constrained_vertices()

      typename std::unordered_set<vertex_descriptor>::iterator it;
      for(it = constrained_vertices_.begin(); it != constrained_vertices_.end(); ++it)
      {
        int i = get(vimap_, *it);
        A.set_coef(i, i, 1.0);

        // also set all cols of the same row = 0 - not needed
        //for(std::size_t j = 0; j<A.column_dimension(); ++j)
        //  A.set_coef(i, j, 0);
      }
    }

    void extract_matrix(Matrix& A)
    {
        std::ofstream out("data/mat.dat");
        for(int i=0; i < A.row_dimension(); ++i)
        {
            for(int j=0; j < A.column_dimension(); ++j)
            {
                NT val = A.get_coef(i, j);
                out<<val<<" ";
            }
            out<<std::endl;
        }
        out.close();
    }

    void extract_vectors(Vector& Vx, Vector& Vy, Vector& Vz)
    {
      CGAL_assertion(Vx.rows() == Vy.rows());
      CGAL_assertion(Vx.rows() == Vz.rows());

      std::ofstream out("data/vecs.dat");
      for(int j=0; j < Vx.dimension(); ++j)
      {
          out<<Vx[j]<<"\t"<<Vy[j]<<"\t"<<Vz[j]<<"\t"<<std::endl;
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

    //-----------
    compute_coeff_matrix(A);
    apply_constraints(A);
    //----------

    //extract_matrix(A);


    //-------------
    compute_rhs(Bx, By, Bz);
    //------------

    //extract_vectors(Bx,By,Bz);


    NT dx, dy, dz;
    if(!solver_.linear_solver(A, Bx, Xx, dx) ||
       !solver_.linear_solver(A, By, Xy, dy) ||
       !solver_.linear_solver(A, Bz, Xz, dz) )
    {
        std::cerr<<"Could not solve linear system!"<<std::endl;
        bool solved = false;
        CGAL_assertion(solved);
    }


    update_mesh(Xx, Xy, Xz);

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
