// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Konstantinos Katrioplas

#ifndef CGAL_OBB_H
#define CGAL_OBB_H

#include <CGAL/Optimal_bounding_box/optimization_algorithms.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <vector>
#include <fstream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Iso_cuboid_3.h>

#include <Eigen/Dense>


namespace CGAL {
namespace Optimal_bounding_box {





template <typename Matrix>
void evolution(Matrix& R, Matrix& points, std::size_t generations) // todo: points is const
{

  CGAL_assertion(points.rows() >= 3);
  CGAL_assertion(points.cols() == 3);
  CGAL_assertion(R.rows() == 3);
  CGAL_assertion(R.cols() == 3);


  std::size_t nelder_mead_iterations = 20;

  Population<Matrix> pop(50);

  //std::cout << "initial pop" << std::endl;
  //pop.show_population();
  //std::cout << std::endl;
  //std::cin.get();

  for(std::size_t t = 0; t < generations; ++t)
  {
    genetic_algorithm(pop, points);

    //std::cout << "pop after genetic" << std::endl;
   // pop.show_population();
   // std::cout << std::endl;
    //std::cin.get();


    for(std::size_t s = 0; s < pop.size(); ++s)
      nelder_mead(pop[s], points, nelder_mead_iterations);

    //std::cout << "pop after nelder mead: " << std::endl;
    //pop.show_population();
    //std::cout << std::endl;
    //std::cin.get();


    // debugging
    /*
    Fitness_map<Matrix> fitness_map(pop, points);
    Matrix R_now = fitness_map.get_best();
    std::cout << "det= " << R_now.determinant() << std::endl;
    */


  }

  // compute fitness of entire population
  Fitness_map<Matrix> fitness_map(pop, points);
  R = fitness_map.get_best();

}


// works on matrices only
template <typename Matrix>
void post_processing(Matrix& points, Matrix& R, Matrix& obb)
{

  // 1) rotate points with R

  Matrix rotated_points(points.rows(), points.cols());
  rotated_points = points * R.transpose();

  // 2) get AABB from rotated points

  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3 Point;
  std::vector<Point> v_points; // Simplex -> std::vector
  for(int i = 0; i < rotated_points.rows(); ++i)
  {
    Point p(rotated_points(i, 0), rotated_points(i, 1), rotated_points(i, 2));
    v_points.push_back(p);
  }
  CGAL::Bbox_3 bbox;
  bbox = bbox_3(v_points.begin(), v_points.end());
  K::Iso_cuboid_3 ic(bbox);

  Matrix aabb(8, 3);
  for(int i = 0; i < 8; ++i)
  {
    aabb(i, 0) = ic[i].x();
    aabb(i, 1) = ic[i].y();
    aabb(i, 2) = ic[i].z();
  }

  // 3) apply inverse rotation to rotated AABB

  obb = aabb * R;


}


/// @param points point coordinates of the input mesh
/// @param obb_points the 8 points of the obb.
template <typename Point>
void find_obb(std::vector<Point>& points, std::vector<Point>& obb_points)
{
  CGAL_assertion(points.size() >= 3);

  if(obb_points.size() != 8) // temp sanity until the API is decided
    obb_points.resize(8);
  CGAL_assertion(obb_points.size() == 8);

  // points: vector -> matrix
  typedef Eigen::MatrixXf MatrixXf;

  MatrixXf points_mat(points.size(), 3);
  for(std::size_t i = 0; i < points.size(); ++i)
  {
    Point p = points[i];
    points_mat(i, 0) = p.x();
    points_mat(i, 1) = p.y();
    points_mat(i, 2) = p.z();
  }


  MatrixXf R(3, 3);
  std::size_t generations = 10;
  CGAL::Optimal_bounding_box::evolution(R, points_mat, generations);

  MatrixXf obb(8, 3);
  CGAL::Optimal_bounding_box::post_processing(points_mat, R, obb);

  // matrix -> vector
  for(std::size_t i = 0; i < 8; ++i)
  {
    Point p(obb(i, 0), obb(i, 1), obb(i, 2));
    obb_points[i] = p;
  }

}

// it is called after post processing
template <typename Matrix>
void matrix_to_mesh_and_draw(Matrix& data_points, std::string filename)
{

  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> Mesh;

  // Simplex -> std::vector
  std::vector<Point> points;

  for(int i = 0; i < data_points.rows(); ++i)
  {
    Point p(data_points(i, 0), data_points(i, 1), data_points(i, 2));
    points.push_back(p);
  }


  Mesh mesh;
  CGAL::make_hexahedron(points[0], points[1], points[2], points[3], points[4], points[5],
      points[6], points[7], mesh);

  std::ofstream out(filename);
  out << mesh;
  out.close();

}




}} // end namespaces






#endif //CGAL_OBB_H



