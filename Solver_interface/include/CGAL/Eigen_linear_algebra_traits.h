
//license

#ifndef CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H
#define CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H

#include <CGAL/basic.h>
#include <Eigen/Dense>
#include <Eigen/QR>


namespace CGAL {


template<typename T, int D = Eigen::Dynamic>
class Eigen_dense_vector;

template<typename T, int D1 = Eigen::Dynamic, int D2 = Eigen::Dynamic>
class Eigen_dense_matrix
{
public:
  typedef Eigen::Matrix<T, D1, D2> EigenType;
  typedef Eigen::Matrix<T, 3, 3> EigenType3;
  typedef Eigen_dense_vector<T, 3> Eigen_dense_vector_type;

  Eigen_dense_matrix(std::size_t nrows, std::size_t ncols)
    : m_matrix(static_cast<int>(nrows), static_cast<int>(ncols))
  {
    CGAL_assertion(m_matrix.rows() > 0);
    CGAL_assertion(m_matrix.cols() > 0);
  }

  Eigen_dense_matrix(int nrows, int ncols)
    : m_matrix(nrows, ncols)
  {
    CGAL_assertion(m_matrix.rows() > 0);
    CGAL_assertion(m_matrix.cols() > 0);
  }

  Eigen_dense_matrix(const EigenType& eigen_mat)
    : m_matrix(eigen_mat) {}

  Eigen_dense_matrix() : m_matrix() {}

  ~Eigen_dense_matrix() {}

  T& operator() (int i_, int j_)
  {
    return m_matrix(i_, j_);
  }

  void set_coef(std::size_t i_, std::size_t j_, T val)
  {
    int i = static_cast<int>(i_);
    int j = static_cast<int>(j_);
    CGAL_precondition(i < m_matrix.rows());
    CGAL_precondition(j < m_matrix.cols());

    m_matrix.coeffRef(i,j) = val;
  }

  std::size_t rows() const {return m_matrix.rows();}
  std::size_t cols() const {return m_matrix.cols();}

  void resize(int i_, int j_) { m_matrix.resize(i_, j_);}

  const T& coeff(int i_) const
  {
    return m_matrix.coeff(i_);
  }

  const Eigen_dense_matrix transpose() const
  {
    Eigen_dense_matrix dm(m_matrix.transpose());
    return dm;
    //return Eigen_dense_matrix(m_matrix.transpose());
  }

  const double determinant() const
  {
    return m_matrix.determinant();
  }

  void qr_factorization(Eigen_dense_matrix& Q)
  {
    Q.resize(m_matrix.rows(), m_matrix.cols());
    Eigen::HouseholderQR<EigenType> qr(m_matrix);
    Eigen_dense_matrix qr_mat(qr.householderQ());
    Q = qr_mat;
  }

  friend void qr_factorization(std::vector<Eigen_dense_matrix>& simplex)
  {
    for(int i = 0; i < simplex.size(); ++i)
    {
      Eigen::HouseholderQR<EigenType> qr(simplex[i].m_matrix);
      Eigen_dense_matrix Q(qr.householderQ());
      simplex[i] = Q;
    }
  }

  friend const Eigen_dense_matrix operator* (const Eigen_dense_matrix& A,
                                             const Eigen_dense_matrix& B)
  {
    const Eigen_dense_matrix product(A.m_matrix * B.m_matrix);
    return product;
  }

  friend const Eigen_dense_matrix operator* (const T& scalar,
                                             const Eigen_dense_matrix& B)
  {
    const Eigen_dense_matrix product(scalar * B.m_matrix);
    return product;
  }

  friend const Eigen_dense_matrix operator+ (const Eigen_dense_matrix& A,
                                             const Eigen_dense_matrix& B)
  {
    const Eigen_dense_matrix sum_result(A.m_matrix + B.m_matrix);
    return sum_result;
  }

  friend const Eigen_dense_matrix operator/ (const Eigen_dense_matrix& A,
                                             const double& scalar)
  {
    const Eigen_dense_matrix product(A.m_matrix / scalar);
    return product;
  }

  friend const Eigen_dense_matrix operator/ (const double& scalar ,
                                             const Eigen_dense_matrix& A)
  {
    const Eigen_dense_matrix product(scalar / A.m_matrix);
    return product;
  }

  mutable EigenType m_matrix;
};


template <typename T, int D>
class Eigen_dense_vector
{
private:
  typedef Eigen::Matrix<T, D, 1> EigenType;

public:

  Eigen_dense_vector(const EigenType&  vec) : m_vector(vec) {}

  const T& coeff(std::size_t i)
  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < 3);
    return m_vector.coeff(i);
  }

  mutable EigenType m_vector;
};


class Eigen_linear_algebra_traits
{
public:
  typedef double NT;

  typedef CGAL::Eigen_dense_matrix<NT> MatrixXd; // dynamic size at run time
  typedef CGAL::Eigen_dense_matrix<NT, 3, 3> Matrix3d; // fixed at compile time
  typedef CGAL::Eigen_dense_vector<NT, 3> Vector3d; // fixed at compile time


  template <class NT, int D1, int D2>
  static CGAL::Eigen_dense_vector<NT, 3> row(const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                              int i)
  {
    return CGAL::Eigen_dense_vector<NT, 3>(A.m_matrix.row(i));
  }

};

typedef CGAL::Eigen_dense_matrix<double> MatrixXd; // dynamic size at run time
typedef CGAL::Eigen_dense_matrix<double, 3, 3> Matrix3d; // fixed at compile time
typedef CGAL::Eigen_dense_vector<double> Vector3d; // fixed at compile time


// matrix multiplication
template <class NT, int D1, int D2, int D3>
const CGAL::Eigen_dense_matrix<NT, D1, D3> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2 >& A,
                                                      const CGAL::Eigen_dense_matrix<NT, D2, D3 >& B)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D3>(A.m_matrix * B.m_matrix);
}

// scalar - matrix multiplication
template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator* (const NT& scalar,
                                                      const CGAL::Eigen_dense_matrix<NT, D1, D2>& B)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2>(scalar * B.m_matrix);
}

template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                                      const NT& scalar)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2>(A.m_matrix * scalar);
}

template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator/ (const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                                      const double& scalar)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2>(A.m_matrix / scalar);
}

template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator/ (const double& scalar,
                                                      const CGAL::Eigen_dense_matrix<NT, D1, D2> & A)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2> (scalar / A.m_matrix);
}

// addition
template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator+ (const CGAL::Eigen_dense_matrix<NT, D1, D2> & A,
                                                      const CGAL::Eigen_dense_matrix<NT, D1, D2> & B)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2> (A.m_matrix + B.m_matrix);
}


// vector - matrix multiplication
template <class NT, int D1, int D2>
const Eigen_dense_vector<NT, D1> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                            const CGAL::Eigen_dense_vector<NT, D2>& V)
{
  return Eigen_dense_vector<NT, D1>(A.m_matrix * V.m_vector);
}




}
// end namespace



#endif // CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H
