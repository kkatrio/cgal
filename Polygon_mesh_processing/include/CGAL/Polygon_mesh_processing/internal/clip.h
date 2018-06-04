// Copyright (c) 2016 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_CLIP_H
#define CGAL_POLYGON_MESH_PROCESSING_CLIP_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/convex_hull_3.h>

namespace CGAL{
namespace Polygon_mesh_processing {

namespace internal
{
template <class TriangleMesh,
          class Ecm,
          class NamedParameters1,
          class NamedParameters2>
bool
clip_open_impl(      TriangleMesh& tm,
                     TriangleMesh& clipper,
               Ecm ecm,
               const NamedParameters1& np_tm,
               const NamedParameters2& np_c)
{
  // first corefine the meshes
  corefine(tm, clipper, np_tm, np_c);

  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::type Vpm;
  typedef typename GetFaceIndexMap<TriangleMesh,
                                   NamedParameters1>::type Fid_map;
  typedef typename GetVertexIndexMap<TriangleMesh,
                                     NamedParameters1>::type Vid_map;

  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters2>::type GeomTraits;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  Fid_map fid_map = boost::choose_param(get_param(np_tm, internal_np::face_index),
                                        get_property_map(boost::face_index, tm));
  Vid_map vid_map = boost::choose_param(get_param(np_tm, internal_np::vertex_index),
                                        get_property_map(boost::vertex_index, tm));
  Vpm vpm1 = boost::choose_param(get_param(np_tm, internal_np::vertex_point),
                                 get_property_map(vertex_point, tm));
  Vpm vpm_c = boost::choose_param(get_param(np_c, internal_np::vertex_point),
                                  get_property_map(vertex_point, clipper));

  // init indices if needed
  helpers::init_face_indices(tm, fid_map);
  helpers::init_vertex_indices(tm, vid_map);

  // set the connected component id of each face
  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));
  std::size_t nb_cc =
    connected_components(tm,
                         bind_property_maps(fid_map, make_property_map(face_cc)),
                         parameters::face_index_map(fid_map).
                         edge_is_constrained_map(ecm));


  boost::dynamic_bitset<> cc_not_handled(nb_cc);
  cc_not_handled.set();
  std::vector <std::size_t> ccs_to_remove;
  /// \todo clipper has been modified, this is not robust if inexact constructions are used
  CGAL::Side_of_triangle_mesh<TriangleMesh, GeomTraits, Vpm>
    side_of(clipper, vpm_c);
  BOOST_FOREACH(face_descriptor f, faces(tm))
  {
    std::size_t cc_id = face_cc[ get(fid_map, f) ];
    if ( !cc_not_handled.test(cc_id) ) continue;

    halfedge_descriptor h=halfedge(f, tm);
    for(int i=0;i<3;++i)
    {
      bool no_marked_edge=true;
      BOOST_FOREACH(halfedge_descriptor h2, halfedges_around_target(h, tm))
        if ( get(ecm, edge(h2, tm)) ){
          no_marked_edge=false;
          break;
        }
      if (no_marked_edge){
        if ( side_of( get(vpm1, target(h, tm) ) ) == ON_UNBOUNDED_SIDE )
          ccs_to_remove.push_back(cc_id);
        cc_not_handled.reset(cc_id);
        break;
      }
      h=next(h, tm);
    }
    if (!cc_not_handled.any()) break;
  }

  if (cc_not_handled.any())
  {
    ///\todo handle using barycenters? won't work for coplanar faces
  }
  //now remove the cc
  remove_connected_components(tm,
    ccs_to_remove,
    bind_property_maps(fid_map, make_property_map(face_cc)),
    np_tm);

  return true;
}

/// \TODO move this to property_map.h?
template <class Set>
struct Constrained_edge_map
{
  typedef boost::read_write_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef typename Set::value_type              key_type;

  Constrained_edge_map()
    : edge_set(NULL)
  {}

  Constrained_edge_map(Set& set)
    : edge_set(&set)
  {}

  friend bool get(const Constrained_edge_map<Set>& map, key_type k)
  {
    return map.edge_set->count(k)!=0;
  }

  friend void put(Constrained_edge_map<Set>& map, key_type k, bool b)
  {
    if (b)
      map.edge_set->insert(k);
    else
      map.edge_set->erase(k);
  }
private:
  Set* edge_set;
};

template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
clip_open_impl(      TriangleMesh& tm,
                     TriangleMesh& clipper,
               boost::param_not_found,
               const NamedParameters1& np_tm,
               const NamedParameters2& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>
    ::edge_descriptor edge_descriptor;
  boost::unordered_set<edge_descriptor> constrained_edges;
  Constrained_edge_map<boost::unordered_set<edge_descriptor> >
    cst_map(constrained_edges);

  return clip_open_impl(tm, clipper,
                        cst_map,
                        np_tm.edge_is_constrained_map(cst_map),
                        np_c);
}

} // end of internal namespace

#ifndef DOXYGEN_RUNNING

///\todo clipper const!
/// requires face_index_map, vertex_index_map for np_tm
/// requires face_index_map for np_c
/// if edge_is_constrained_map is not provided in np_tm a default one is
/// provided using boost::unordered_set<edge_descriptor>
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
clip(      TriangleMesh& tm,
     /*const*/ TriangleMesh& clipper,
           bool close,
     const NamedParameters1& np_tm,
     const NamedParameters2& np_c)
{
  if (close && is_closed(tm))
    return corefine_and_compute_intersection(tm, clipper, tm, np_tm, np_c);

  return internal::clip_open_impl(tm, clipper,
      get_param(np_tm, internal_np::edge_is_constrained), np_tm, np_c);
}

template <class Point_3, class TriangleMesh>
void construct_triangle_mesh(std::map<int, Point_3>& pid_map, TriangleMesh& tm,
                             cpp11::array<CGAL::Oriented_side,8>& orientations,
                             std::map<int, int>& fragments)
{

  cpp11::array<int, 8> edge_indices = {{ 0,1, 1,2, 2,3, 3,0 }};

  std::vector<std::vector<int>> face_indices;
  face_indices.push_back({{ 0,1,2,3 }});
  face_indices.push_back({{ 4,5,6,7 }});
  face_indices.push_back({{ 0,1,5,4 }});
  face_indices.push_back({{ 1,2,6,5 }});
  face_indices.push_back({{ 2,3,7,6 }});
  face_indices.push_back({{ 3,0,4,7 }});

  // collect intact faces
  std::vector<std::vector<int>> pfaces;
  for(int k = 0 ; k < 6; ++k)
  {
    bool is_on_negative_side = true;
    is_on_negative_side &= orientations[face_indices[k][0]] == ON_NEGATIVE_SIDE;
    is_on_negative_side &= orientations[face_indices[k][1]] == ON_NEGATIVE_SIDE;
    is_on_negative_side &= orientations[face_indices[k][2]] == ON_NEGATIVE_SIDE;
    is_on_negative_side &= orientations[face_indices[k][3]] == ON_NEGATIVE_SIDE;

    if(is_on_negative_side)
      pfaces.push_back(face_indices[k]);
  }
  std::cout << "faces_intact = " << pfaces.size() << std::endl;


  // collect cut faces
  std::vector<std::vector<int>> faces_cut;
  for(int k = 0 ; k < 6; ++k)
  {
    for(int j = 0; j < 4; ++j)
    {
      int s = face_indices[k][edge_indices[2*j]];
      int t = face_indices[k][edge_indices[2*j + 1]];
      if(orientations[s] != orientations[t])
      {
        faces_cut.push_back(face_indices[k]);
        break;
      }
    }
  }
  std::cout << "faces_cut = " << faces_cut.size() << std::endl;


  // repair cut faces
  for(std::vector<int>& f : faces_cut)
  {
    for(std::size_t i = 0; i < f.size(); ++i)
    {
      // remove bad vertices
      int vertex_id = f[i];
      //std::cout << "vertex_id = " << vertex_id << "\n";
      if(orientations[vertex_id] == ON_POSITIVE_SIDE)
      {
        assert(fragments.find(vertex_id) != fragments.end());
        // replace with the intersection vertex
        //std::cout << "fragments[vertex_id]= " << fragments[vertex_id] << "\n";
        f[i] = fragments[vertex_id];
      }
      
    }
  }
  
  // has all faces of the new mesh
  pfaces.insert(pfaces.end(), faces_cut.begin(), faces_cut.end());  
  
  for(auto f : pfaces)
  {
    for(int i : f)
    {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
  
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  
  std::cout << "#faces before= " << faces(tm).size() << std::endl;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::const_type Vpm;
  Vpm vpm = get(boost::vertex_point, tm);



  for(auto it = pid_map.begin(); it != pid_map.end(); ++it)
  {
    std::cout << it->first << " " << it->second << std::endl;
  }


  // connect in pmesh
  for(auto pf : pfaces)
  {

    std::vector<halfedge_descriptor> hedges;

    face_descriptor f = add_face(tm);

    // for every vertex - id  on the face
    for(int id : pf)
    {
      // add vertex to tm and point to vpm
      vertex_descriptor v = add_vertex(tm);
      Point_3 p = pid_map[id];
      put(vpm, v, p);

      // add halfedge and set its target
      halfedge_descriptor h = halfedge(add_edge(tm), tm);
      set_target(h, v, tm);
      hedges.push_back(h); // for last set_halfedge

      // set face and halfedge
      set_face(h, f, tm);
      set_halfedge(v, h, tm); // needed ?
    }

    set_halfedge(f, hedges[0], tm);


    // connect halfedges
    hedges.push_back(hedges.front());
    for(std::size_t i = 0; i < hedges.size() - 1; ++i)
    {
      set_next(hedges[i], hedges[i + 1], tm);
    }
  }
  std::cout << "connecting done.\n";


  std::cout << "#vertices= " << vertices(tm).size() << std::endl;
  std::cout << "#faces= " << faces(tm).size() << std::endl;




  std::cin.get();


  
  triangulate_faces(faces(tm), tm);
  std::cout << "#faces after= " << faces(tm).size() << std::endl;
  
  
  
}



/// \todo document me
template <class Plane_3,
          class TriangleMesh,
          class NamedParameters>
Oriented_side
clip_to_bbox(const Plane_3& plane,
             const Bbox_3& bbox,
                   TriangleMesh& tm_out,
             const NamedParameters& np )
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::type Vpm;

  Vpm vpm_out = boost::choose_param(get_param(np, internal_np::vertex_point),
                                    get_property_map(boost::vertex_point, tm_out));

  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  cpp11::array<Point_3,8> corners= {{
    Point_3(bbox.xmin(),bbox.ymin(),bbox.zmin()),
    Point_3(bbox.xmin(),bbox.ymax(),bbox.zmin()),
    Point_3(bbox.xmax(),bbox.ymax(),bbox.zmin()),
    Point_3(bbox.xmax(),bbox.ymin(),bbox.zmin()),
    Point_3(bbox.xmin(),bbox.ymin(),bbox.zmax()),
    Point_3(bbox.xmin(),bbox.ymax(),bbox.zmax()),
    Point_3(bbox.xmax(),bbox.ymax(),bbox.zmax()),
    Point_3(bbox.xmax(),bbox.ymin(),bbox.zmax())
  }};

  cpp11::array<CGAL::Oriented_side,8> orientations = {{
    plane.oriented_side(corners[0]),
    plane.oriented_side(corners[1]),
    plane.oriented_side(corners[2]),
    plane.oriented_side(corners[3]),
    plane.oriented_side(corners[4]),
    plane.oriented_side(corners[5]),
    plane.oriented_side(corners[6]),
    plane.oriented_side(corners[7])
  }};

  std::vector<Point_3> points;

  // look for intersections on edges
  cpp11::array<int,24> edge_indices = {{ // 2 *12 edges
    0,1, 1,2, 2,3, 3,0, // bottom face edges
    4,5, 5,6, 6,7, 7,4, // top face edges
    0,4, 1,5, 2,6, 3,7
  }};
  
  std::map<int, int> fragments; // <good vertex, intersection>
  std::vector<int> pids;
  std::map<int, Point_3> pid_map;

  int id = 8; // id of new intersetion vertices

  for (int i=0; i<12; ++i)
  {
    int i1=edge_indices[2*i], i2=edge_indices[2*i+1];
    if (orientations[i1]==ON_ORIENTED_BOUNDARY) continue;
    if (orientations[i2]==ON_ORIENTED_BOUNDARY) continue;
    if (orientations[i1]!=orientations[i2])
    {
      points.push_back(
        boost::get<Point_3>(
          *intersection(plane, Segment_3(corners[i1], corners[i2]) )
        )
      );

      pid_map[id] = points.back(); // todo: avoid points vec

      int bad_vertex = orientations[i1] == ON_POSITIVE_SIDE ? i1 : i2;
      // fragmented edges indices
      fragments[bad_vertex] = id; // not enough: 0,1 and 0,4

      id++;


      int good_vertex = orientations[i1] == ON_NEGATIVE_SIDE ? i1 : i2;
      pid_map[good_vertex] = corners[good_vertex];


      /*

      // update the vpm_out & tm_out with the intersection points
      vertex_descriptor v = add_vertex(tm_out);
      Point_3 p = boost::get<Point_3>(
          *intersection(plane, Segment_3(corners[i1], corners[i2]) ) );
      put(vpm_out, v, p);

      // put in vpm the coors of the good corners as well.
      int good_vertex = orientations[i1] == ON_NEGATIVE_SIDE ? i1 : i2;
      vertex_descriptor vg = add_vertex(tm_out);
      Point_3 pg = corners[good_vertex];
      std::cout << pg.x() << " " << pg.y() << " " << pg.z() << std::endl;
      put(vpm_out, vg, pg);

      */

    }
  }


  /*
  std::cout << "fragments:\n";
  for(auto i = fragments.begin(); i != fragments.end(); ++i)
  {
    std::cout << i->first << " " << i->second << "\n";
  }
  */

  Oriented_side last_os = ON_ORIENTED_BOUNDARY;
  for (int i=0; i<8; ++i)
    if (orientations[i]!=ON_ORIENTED_BOUNDARY)
    {
      if (last_os==ON_ORIENTED_BOUNDARY)
        last_os=orientations[i];
      else
      {
        if(last_os!=orientations[i])
        {
          last_os=ON_ORIENTED_BOUNDARY;
          break;
        }
      }
    }

  // the intersection is the full bbox
  if (last_os!=ON_ORIENTED_BOUNDARY)
    return last_os;

  //add points on negative side and on the plane
  for (int i=0; i<8; ++i)
    if (orientations[i]!=ON_POSITIVE_SIDE)
      points.push_back(corners[i]);

  // take the convex hull of the points on the negative side+intersection points
  // overkill...
  //Polyhedron_3<Geom_traits> P;
  //CGAL::convex_hull_3(points.begin(), points.end(), P);



  construct_triangle_mesh(pid_map, tm_out, orientations, fragments);



  //face_descriptor f = Euler::add_face(vertices_out, tm_out);

  /*
  std::cout << "tm_out: \n";
  std::cout << "vertices.size()= " << vertices(tm_out).size() << std::endl;
  for(auto v : vertices(tm_out))
  {
    Point_3 p = get(vpm_out, v);
    std::cout << p.x() << " " << p.y() << " " << p.z() << std::endl;
  }
  */

  //triangulate_face(f, tm_out);





  /*
  copy_face_graph(P, tm_out,
                  Emptyset_iterator(), Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, P), vpm_out);
*/




  return ON_ORIENTED_BOUNDARY;
}


// convenience overload
template <class TriangleMesh,
          class NamedParameters1>
bool
clip(      TriangleMesh& tm,
     /*const*/ TriangleMesh& clipper,
           bool close,
     const NamedParameters1& np_tm)
{
  return clip(tm, clipper, close, np_tm, parameters::all_default());
}

// convenience overload
template <class TriangleMesh>
bool
clip(      TriangleMesh& tm,
     /*const*/ TriangleMesh& clipper,
           bool close)
{
  return clip(tm, clipper, close, parameters::all_default());
}

// works only with the default point map, for more complex use cases, use
// clip_to_bbox first and the other overload of clip with two meshes
/// \todo document me
template <class TriangleMesh,
          class Plane_3>
void clip(      TriangleMesh& tm,
          const Plane_3& plane,
          bool close)
{
  if( boost::begin(faces(tm))==boost::end(faces(tm)) ) return;
  CGAL::Bbox_3 bbox = ::CGAL::Polygon_mesh_processing::bbox(tm);
  //extend the bbox a bit to avoid border cases
  double xd=(bbox.xmax()-bbox.xmin())/100;
  double yd=(bbox.ymax()-bbox.ymin())/100;
  double zd=(bbox.zmax()-bbox.zmin())/100;
  bbox=CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                    bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);
  TriangleMesh clipper;
  Oriented_side os = clip_to_bbox(plane, bbox, clipper, parameters::all_default());

  switch(os)
  {
    case ON_NEGATIVE_SIDE:
      return; // nothing to clip, the full mesh is on the negative side
    case ON_POSITIVE_SIDE:
      clear(tm); // clear the mesh that is fully on the positive side
      return;
    default:
      clip(tm, clipper, close);
  }
}

#endif // !DOXYGEN_RUNNING

} } //end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_CLIP_H
