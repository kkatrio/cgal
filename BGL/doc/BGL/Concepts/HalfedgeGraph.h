/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `HalfedgeGraph` is a refinement of the \sc{Bgl} concept
<A HREF="http://www.boost.org/libs/graph/doc/Graph.html">`Graph`</A> and adds the notion of a *halfedge*: Each edge is
associated with two *opposite* halfedges with source and target vertices swapped.
Furthermore, halfedges have a *successor* and *predecessor*,
and form cycles we call *faces*. However, this concept 
does not introduce a face type. 
A `HalfedgeGraph` is undirected and does not allow parallel edges.

Using the composition of the *successor* and *opposite* functions results 
in another cycle, namely the cycle of halfedges which are incident to
the same vertex. We refer to \ref PkgBGLIterators for a description of
iterators and circulators for these halfedge cycles.


\cgalRefines <A HREF="http://www.boost.org/libs/graph/doc/IncidenceGraph.html">`IncidenceGraph`</A>
\cgalRefines <A HREF="http://www.boost.org/libs/graph/doc/PropertyGraph.html">`PropertyGraph`</A>

A model of `HalfedgeGraph` must have the interior property `vertex_point` attached to its vertices.

\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`
*/
class HalfedgeGraph {};

/*! \relates HalfedgeGraph
returns the edge corresponding to halfedges `h` and `opposite(h,g)`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::edge_descriptor
edge(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns one of the halfedges corresponding to `e`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
halfedge(boost::graph_traits<HalfedgeGraph>::edge_descriptor f, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns a halfedge with target `v`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
halfedge(boost::graph_traits<HalfedgeGraph>::vertex_descriptor v, const HalfedgeGraph& g);


/*! \relates HalfedgeGraph
returns the halfedge with source `u` and target `v`. The Boolean is `true`, iff this halfedge exists.
 */
template <typename HalfedgeGraph>
std::pair<boost::graph_traits<HalfedgeGraph>::halfedge_descriptor,bool>
halfedge(boost::graph_traits<HalfedgeGraph>::vertex_descriptor u,
         boost::graph_traits<HalfedgeGraph>::vertex_descriptor v,
         const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the halfedge with source and target swapped.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
opposite(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the source vertex of `h`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::vertex_descriptor
source(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the target vertex of `h`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::vertex_descriptor
target(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the next halfedge around its face.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
next(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the previous halfedge around its face.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
prev(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns a special halfedge that is not equal to any other halfedge.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
null_halfedge();
