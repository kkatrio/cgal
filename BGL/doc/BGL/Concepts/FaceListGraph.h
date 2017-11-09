/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `FaceListGraph` refines the concept `FaceGraph` and adds
the requirement for traversal of all faces in a graph.

\cgalRefines `FaceGraph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`
\cgalHasModel `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
class FaceListGraph{};


/*! \relates FaceListGraph
 * returns an iterator range over all faces.
 */
template <typename FaceListGraph>
std::pair<boost::graph_traits<FaceListGraph>::face_iterator,
          boost::graph_traits<FaceListGraph>::face_iterator>
faces(const FaceListGraph& g);


/*! \relates FaceListGraph
  returns an upper bound of the number of faces of the graph.
  \attention `num_faces()` may return a number larger than `std::distance(faces(g).first, faces(g).second)`.
  This is the case for implementations only marking faces deleted in the face container.
 */
template <typename FaceListGraph>
boost::graph_traits<FaceListGraph>::face_size_type
num_faces(const FaceListGraph& g);

