#ifndef SHRINK_H
#define SHRINK_H

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>


namespace CGAL {

namespace Polygon_mesh_processing {



template<typename PolygonMesh, typename NamedParameters>
void shrink(PolygonMesh& pmesh, const double& shrink_factor, const NamedParameters& np)
{

    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    //vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    //geom_traits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

    typedef typename GeomTraits::Point_3 Point;


    std::map<face_descriptor, Point> face_centroid_map;

    // find each face's centroid
    BOOST_FOREACH(face_descriptor f, faces(pmesh))
    {

        Point points[3];
        halfedge_descriptor he = halfedge(f, pmesh);
        for(std::size_t i=0; i<3; ++i)
        {
            points[i] = get(vpmap, target(he, pmesh));
            he = next(he, pmesh);
        }

        Point centroid = CGAL::centroid(points[0], points[1], points[2]);

        face_centroid_map[f] = centroid;

    }

    // perform movements
    typename std::map<face_descriptor, Point>::iterator it;
    for(it = face_centroid_map.begin(); it != face_centroid_map.end(); ++it)
    {
        face_descriptor f = it->first;
        Point centroid = it->second;
        BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f, pmesh), pmesh))
        {
            Point p = get(vpmap, v);
            typename GeomTraits::Vector_3 vec(p, centroid);
            Point new_p = p + shrink_factor * vec;


            // create new vertex




            //put(vpmap, v, new_p);
        }

    }

}






template<typename PolygonMesh, typename NamedParameters>
PolygonMesh real_shrink(PolygonMesh& pmesh, const double& shrink_factor, const NamedParameters& np)
{
    // create a new graph
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point3;
    typedef CGAL::Surface_mesh<Point3> Surface_mesh;

    //Surface_mesh new_pmesh;
    PolygonMesh new_pmesh;

    // create new_vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap new_vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, new_pmesh));

    // create existing graph's vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    //-------------------------------------------------------------------------------

    // geom_traits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;
    typedef typename GeomTraits::Point_3 Point;

    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    //-------------------------------------------------------------------------------

    // find each face's centroid
    std::map<face_descriptor, Point> face_centroid_map;
    BOOST_FOREACH(face_descriptor f, faces(pmesh))
    {

        Point points[3];
        halfedge_descriptor he = halfedge(f, pmesh);
        for(std::size_t i=0; i<3; ++i)
        {
            points[i] = get(vpmap, target(he, pmesh));
            he = next(he, pmesh);
        }

        Point centroid = CGAL::centroid(points[0], points[1], points[2]);

        face_centroid_map[f] = centroid;

    }
    unsigned int i=0;
    // create new shrunk vertices towards the centroid and add them to the new graph
    typename std::map<face_descriptor, Point>::iterator it;
    for(it = face_centroid_map.begin(); it != face_centroid_map.end(); ++it)
    {
        std::cout<<"face"<<i<<std::endl;
        face_descriptor f = it->first;
        Point centroid = it->second;
        std::vector<vertex_descriptor> vertices;
        BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f, pmesh), pmesh))
        {
            // calc new location
            Point p = get(vpmap, v);
            typename GeomTraits::Vector_3 vec(p, centroid);
            Point new_p = p + shrink_factor * vec;

            // add new vertex in the new pmesh
            vertex_descriptor new_v = add_vertex(new_pmesh);
            put(new_vpmap, new_v, new_p);

            vertices.push_back(new_v);
        }

        // add new face
        CGAL_assertion_code(face_descriptor fd = )
        CGAL::Euler::add_face(vertices, new_pmesh);
        CGAL_assertion(fd != boost::graph_traits<PolygonMesh>::null_face());

        std::cout<<"done with face"<<i<<std::endl;
        i++;


    }

    return new_pmesh;

}








} // namespace Polygon_mesh_processing
} // namespace CGAL














#endif // SHRINK_H
