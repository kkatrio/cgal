#include <QtCore/qglobal.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include "Scene_surface_mesh_item.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "Scene_polyhedron_selection_item.h"

#include <boost/graph/graph_traits.hpp>
#include <CGAL/Optimal_bounding_box/obb.h>

//typedef Scene_surface_mesh_item Scene_facegraph_item;
//typedef Scene_facegraph_item Scene_facegraph_item;

typedef Scene_facegraph_item::Face_graph FaceGraph;
typedef Polyhedron::Point_3 Point_3;
using namespace CGAL::Three;

class Create_obb_mesh_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const;

  bool applicable(QAction*) const {

    if (scene->selectionIndices().size() == 1)
    {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()))
    || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
    }

    Q_FOREACH(int index, scene->selectionIndices())
    {
      if (qobject_cast<Scene_facegraph_item*>(scene->item(index)))
        return true;
    }
    return false;

    /*
    if(scene->mainSelectionIndex() != -1
       && scene->item(scene->mainSelectionIndex())->isFinite())
      return true;
  return false;
  */

}

protected:
  void gather_mesh_points(std::vector<Point_3>& points);
  void obb();

public Q_SLOTS:
  void createObb() {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    obb();
    QApplication::restoreOverrideCursor();
  }

private:
  Scene_interface* scene;
  QMainWindow* mw;
  QAction* actionObb;

}; // end Create_obb_mesh_plugin class



void Create_obb_mesh_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionObb = new QAction(tr("Create &Optimal Bbox Mesh"), mainWindow);
  actionObb->setObjectName("createObbMeshAction");
  connect(actionObb, SIGNAL(triggered()), this, SLOT(createObb()));
}

QList<QAction*> Create_obb_mesh_plugin::actions() const {
  return QList<QAction*>() << actionObb;
}


void Create_obb_mesh_plugin::gather_mesh_points(std::vector<Point_3>& points)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_facegraph_item* poly_item =
    qobject_cast<Scene_facegraph_item*>(scene->item(index));

  Scene_polyhedron_selection_item* selection_item =
    qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

  if(poly_item || selection_item)
  {
    typedef typename boost::property_map<FaceGraph, boost::vertex_point_t>::type PointPMap;
    typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

    std::vector<vertex_descriptor> selected_vertices;

    if(poly_item != NULL)
    {
      FaceGraph& pmesh = *poly_item->polyhedron();
      selected_vertices.assign(vertices(pmesh).begin(), vertices(pmesh).end());
      PointPMap pmap = get(CGAL::vertex_point, pmesh);
      BOOST_FOREACH(vertex_descriptor v, selected_vertices)
        points.push_back(get(pmap, v ));

    }
    else if(selection_item != NULL) // using selection of faces
    {
      FaceGraph& pmesh = *selection_item->polyhedron();
      BOOST_FOREACH(face_descriptor f, selection_item->selected_facets)
      {
        BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f, pmesh), pmesh))
        {
          selected_vertices.push_back(v);
        }
      }

      PointPMap pmap = get(CGAL::vertex_point, pmesh);
      BOOST_FOREACH(vertex_descriptor v, selected_vertices)
        points.push_back(get(pmap, v ));
    }
    CGAL_assertion(points.size() >= 3);
  }
}

void Create_obb_mesh_plugin::obb()
{
  // gather point coordinates
  std::vector<Point_3> points;
  gather_mesh_points(points);

  // find obb
  std::vector<Point_3> obb_points(8);
  CGAL::Optimal_bounding_box::find_obb(points, obb_points, true);

  Scene_item* item;
  if(mw->property("is_polyhedorn_mode").toBool())
  {
    Polyhedron* p = new Polyhedron;
    CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                          obb_points[4], obb_points[5], obb_points[6], obb_points[7], *p);
    item = new Scene_polyhedron_item(p);
  }
  else {
    SMesh* p = new SMesh;
    CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                          obb_points[4], obb_points[5], obb_points[6], obb_points[7], *p);
    item = new Scene_surface_mesh_item(p);
  }

  item->setName("Optimal bbox mesh");
  item->setRenderingMode(Wireframe);
  scene->addItem(item);
}

#include "Create_obb_mesh_plugin.moc"
