
#include <QtCore/qglobal.h>
#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/utility.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>

#include "ui_Smoothing_plugin.h"


//using namespace CGAL::Polygon_mesh_processing;
using namespace CGAL::Three;
class Polyhedron_demo_smothing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")



public:
    void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
    {
        scene = scene_interface;
        mw = mainWindow;

        actionSmoothing_ = new QAction(tr("Smoothing"), mw);
        actionSmoothing_->setProperty("subMenuName", "Polygon Mesh Processing");

        connect(actionSmoothing_, SIGNAL(triggered()), this, SLOT(smoothing_action()));

        dock_widget = new QDockWidget("Smoothing", mw);
        dock_widget->setVisible(false);

        ui_widget.setupUi(dock_widget);
        addDockWidget(dock_widget);

        //connect(ui_widget.Apply_button,  SIGNAL(clicked()), this, SLOT(on_Apply_by_type_clicked()));
        connect(ui_widget.shape_smoothing_button,  SIGNAL(clicked()), this, SLOT(on_Apply_smoothing_clicked()));
        connect(ui_widget.setup_system_button,  SIGNAL(clicked()), this, SLOT(on_Setup_system_clicked()));
        //connect(ui_widget.Run_convergence_button,  SIGNAL(clicked()), this, SLOT(on_Run_convergence_clicked()));

    }

    QList<QAction*> actions() const
    {
        return QList<QAction*>() << actionSmoothing_;
    }

    bool applicable(QAction*) const
    {
      const Scene_interface::Item_id index = scene->mainSelectionIndex();
      if (qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
        return true;
      else if (qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)))
        return true;
      else
        return false;
    }

    virtual void closure()
    {
      dock_widget->hide();
    }

private:
    void init_ui()
    {
        ui_widget.smothing_spinBox->setValue(1);
        ui_widget.smothing_spinBox->setSingleStep(10);
        ui_widget.smothing_spinBox->setMinimum(1);

        /*
        ui_widget.gd_dSpinBox->setSingleStep(0.0001);
        ui_widget.gd_dSpinBox->setDecimals(4);
        ui_widget.gd_dSpinBox->setMinimum(0.0001);
        ui_widget.gd_dSpinBox->setValue(0.001);

        ui_widget.use_weights_checkBox->setChecked(true);

        ui_widget.dist_dSpinBox->setValue(0.01);
        ui_widget.dist_dSpinBox->setSingleStep(0.0001);
        ui_widget.dist_dSpinBox->setDecimals(4);
        ui_widget.dist_dSpinBox->setMinimum(0.0001);

        ui_widget.gd_precision_label->setToolTip("Tradeoff between precision and speed. Less is more precise.");
        ui_widget.distance_label->setToolTip("Tradeoff between precision and speed. Less is more precise.");

        ui_widget.iterations_spinBox->setValue(20);
        ui_widget.iterations_spinBox->setSingleStep(1);
        ui_widget.iterations_spinBox->setMinimum(1);

        ui_widget.curv_iterations_spinBox->setValue(1);
        ui_widget.curv_iterations_spinBox->setSingleStep(1);
        ui_widget.curv_iterations_spinBox->setMinimum(1);
        */
    }

    void init_parameters()
    {

    }


public Q_SLOTS:
    void smoothing_action()
    {
        dock_widget->show();
        dock_widget->raise();

        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));

        if(poly_item)
        {
            init_ui();
        }
    }

    void on_Setup_system_clicked()
    {
       const Scene_interface::Item_id index = scene->mainSelectionIndex();
       Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
       Polyhedron& pmesh = *poly_item->polyhedron();

       QApplication::setOverrideCursor(Qt::WaitCursor);

       CGAL::Polygon_mesh_processing::setup_system(pmesh, stiffness_matrix);

       QApplication::restoreOverrideCursor();

    }


    void on_Apply_smoothing_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Polyhedron& pmesh = *poly_item->polyhedron();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        unsigned int nb_iter = ui_widget.smothing_spinBox->value();
        std::cout << "nb_iter= " << nb_iter << std::endl;
        CGAL::Polygon_mesh_processing::smooth_shape(pmesh, nb_iter, stiffness_matrix);

        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();

        QApplication::restoreOverrideCursor();
    }

    /*
    void on_Run_convergence_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Polyhedron& pmesh = *poly_item->polyhedron();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        unsigned int nb_iter = ui_widget.iterations_spinBox->value();
        double dist = ui_widget.dist_dSpinBox->value();
        double gd_precision = ui_widget.gd_dSpinBox->value();
        bool use_weights = ui_widget.use_weights_checkBox->isChecked();

        compatible_smoothing(pmesh,
                             parameters::number_of_iterations(nb_iter).
                             distance_precision(dist).
                             gradient_descent_precision(gd_precision).
                             use_weights(use_weights));

        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();

        QApplication::restoreOverrideCursor();
    }
    */




private:
    QAction* actionSmoothing_;
    QDockWidget* dock_widget;
    Ui::Smoothing ui_widget;

    typedef typename Eigen::SparseMatrix<double> Eigen_matrix;
    Eigen_matrix stiffness_matrix;



};



#include "Smoothing_plugin.moc"


