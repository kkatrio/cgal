// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_widget_toolbar_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_TOOLBAR_LAYERS_H
#define CGAL_QT_WIDGET_TOOLBAR_LAYERS_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/IO/Qt_widget.h>

#include <CGAL/IO/Qt_layer_show_mouse_coordinates.h>

#include "Qt_layer_show_triangulation.h"
#include "Qt_layer_show_voronoi.h"
#include "Qt_layer_show_points.h"
#include "Qt_layer_show_nearest_vertex.h"
#include "Qt_layer_show_circum_circle.h"

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qbuttongroup.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rp;
typedef CGAL::Delaunay_triangulation_2<Rp>  Delaunay;

namespace CGAL {

class Layers_toolbar : public QObject
{
	Q_OBJECT
public:
  Layers_toolbar(Qt_widget *w, QMainWindow *mw, Delaunay *t);
  ~Layers_toolbar()
  {
    delete showT;
    delete showV;
    delete showP;
    delete showNV;
    delete showMC;
  };
  QToolBar*	toolbar(){return maintoolbar;};

signals:
  void new_object(CGAL::Object);
	
private:
  QToolBar	*maintoolbar;
  QToolButton	*but[10];
  Qt_widget	*widget;
  QMainWindow	*window;
  Delaunay	*dt;
  QButtonGroup  *button_group;

  int	  nr_of_buttons;
	
  CGAL::Qt_layer_show_triangulation < Delaunay >  *showT;
  CGAL::Qt_layer_show_voronoi < Delaunay >	  *showV;
  CGAL::Qt_layer_show_points < Delaunay >	  *showP;
  CGAL::Qt_layer_nearest_vertex < Delaunay >	  *showNV;
  CGAL::Qt_layer_mouse_coordinates		  *showMC;
  CGAL::Qt_layer_circum_circle < Delaunay>        *showCC;
};//end class

};//end namespace

#endif
