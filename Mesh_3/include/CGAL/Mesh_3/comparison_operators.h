// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://jtournoi@scm.gforge.inria.fr/svnroot/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/Mesh_complex_3_in_triangulation_3.h $
// $Id: comparison_operators.h  $
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef COMPARISON_OPERATORS_H
#define COMPARISON_OPERATORS_H

namespace CGAL {
  namespace Mesh_3 {

    template<typename Vertex_handle>
    struct Vertex_handle_comparator 
      : public std::binary_function<Vertex_handle, Vertex_handle, bool> 
    {
      bool operator()(const Vertex_handle& v1, const Vertex_handle& v2) const 
      { 
        return v1->point() < v2->point(); 
      }
    };

    template<typename Tr>
    struct Cell_handle_comparator
      : public std::binary_function<typename Tr::Cell_handle,
                                    typename Tr::Cell_handle,
                                    bool>
    {
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Vertex_handle Vertex_handle;

      bool operator()(const Cell_handle& c1, const Cell_handle& c2) const
      {
        std::vector<Vertex_handle> v1;
        std::vector<Vertex_handle> v2;
        for(int i = 0; i < 4; ++i)
        {
          v1.push_back(c1->vertex(i));
          v2.push_back(c2->vertex(i));
        }
        Vertex_handle_comparator<Vertex_handle> vcomp;
        std::sort(v1.begin(), v1.end(), vcomp);
        std::sort(v2.begin(), v2.end(), vcomp);
        for(std::size_t i = 0; i < 4; ++i)
        {
          if(v1[i] == v2[i])
            continue;
          else return vcomp(v1[i], v2[i]);
        }
        return false;
      }
    };

    template<typename Tr>
    struct Triangulation_canonical_facets_comparator
      : public std::binary_function<typename Tr::Facet,
                                    typename Tr::Facet,
                                    bool>
    {
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Vertex_handle Vertex_handle;

      bool operator()(const Facet& f1, const Facet& f2) const
      {
        Vertex_handle v1 = (f1.first)->vertex(f1.second);
        Vertex_handle v2 = (f2.first)->vertex(f2.second);
        Vertex_handle_comparator<Vertex_handle> vcomp;
        if(vcomp(v1,v2))
          return true;
        else if(v1 == v2)
        {
          std::vector<Vertex_handle> vf1;
          std::vector<Vertex_handle> vf2;
          for(int i = 0; i < 4; ++i)
          {
            if(i != f1.second)
              vf1.push_back(f1.first->vertex(i));
            if(i != f2.second)
              vf2.push_back(f2.first->vertex(i));
          }
          std::sort(vf1.begin(), vf1.end(), vcomp);
          std::sort(vf2.begin(), vf2.end(), vcomp);
          for(std::size_t i = 0; i < 3; ++i)
          {
            if(vf1[i] == vf2[i])
              continue;
            else return vcomp(vf1[i], vf2[i]);
          }
          return false;
        }
        else return false;        
      }
    };

    template<typename Tr>
    struct Triangulation_Finite_facets_comparator
      : public std::binary_function<typename Tr::Facet,
                                    typename Tr::Facet,
                                    bool>
    {
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Vertex_handle Vertex_handle;

      Triangulation_Finite_facets_comparator(const Tr& tr)
        : tr_(tr) {}

      bool operator()(const Facet& f1, const Facet& f2) const
      {
        return Triangulation_canonical_facets_comparator<Tr>(
          canonical(f1),
          canonical(f2));
      }

      Facet canonical(const Facet& f) const
      {
        const Facet mirror = tr_.mirror_facet(f);
        Vertex_handle vf = f.first->vertex(f.second);
        Vertex_handle vm = mirror.first->vertex(mirror.second);
        if(tr_.is_infinite(vf))
          return mirror;
        else if(tr_.is_infinite(vm))
          return f;

        Vertex_handle_comparator<Vertex_handle> vcomp;
        return vcomp(vf,vm) ? f : mirror;
      }
    
    private:
      const Tr& tr_;
    };

    template<typename Facet>
    struct Polyhedron_Facet_handle_comparator  
      : public std::binary_function<typename Facet::Facet_handle, 
                                    typename Facet::Facet_handle, 
                                    bool> 
    {
      typedef typename Facet::Facet_handle  Facet_handle;
      typedef typename Facet::Vertex_handle Vertex_handle;
      typedef typename Facet::Halfedge_around_facet_circulator Facet_he_circ;

      bool operator()(const Facet_handle& pf1, const Facet_handle& pf2) const
      { 
        //collect vertices of both facets
        std::vector<Vertex_handle> vertices_f1;
        std::vector<Vertex_handle> vertices_f2;
        Facet_he_circ begin =	pf1->facet_begin();
        Facet_he_circ end = begin;
        do
        {
          vertices_f1.push_back(begin->vertex()); 
          ++begin;
        }while(begin != end);

        begin =	pf2->facet_begin();
        end = begin;
        do
        {
          vertices_f2.push_back(begin->vertex()); 
          ++begin;
        }while(begin != end);

        //compare vertices
        Vertex_handle_comparator<Vertex_handle> comparator;
        std::sort(vertices_f1.begin(), vertices_f1.end(), comparator);
        std::sort(vertices_f2.begin(), vertices_f2.end(), comparator);
        std::size_t nmax = std::min(vertices_f1.size(), vertices_f2.size());
        for(std::size_t i = 0; i < nmax; i++)
        {
          if(vertices_f1[i] == vertices_f2[i])
            continue;
          else //either < or >
            return comparator(vertices_f1[i], vertices_f2[i]);
        }
        return false; //it is the same facet
      }
    };

    template<typename Halfedge>
    struct Halfedge_handle_comparator 
      : public std::binary_function<typename Halfedge::Halfedge_handle,
                                    typename Halfedge::Halfedge_handle,
                                    bool> 
    {
      typedef typename Halfedge::Halfedge_handle Halfedge_handle;
      typedef typename Halfedge::Vertex_handle Vertex_handle;

      bool operator()(const Halfedge_handle& he1, const Halfedge_handle& he2) const
      {     
        //collect vertices of both facets
        std::vector<Vertex_handle> vertices_he1;
        vertices_he1.push_back(he1->vertex());
        vertices_he1.push_back(he1->opposite()->vertex());
      
        std::vector<Vertex_handle> vertices_he2;
        vertices_he2.push_back(he2->vertex());
        vertices_he2.push_back(he2->opposite()->vertex());

        //compare vertices
        Vertex_handle_comparator<Vertex_handle> comparator;
        std::sort(vertices_he1.begin(), vertices_he1.end(), comparator);
        std::sort(vertices_he2.begin(), vertices_he2.end(), comparator);

        //we want he and he->opposite() to be both in the set, for flooding
        if(he1 == he2->opposite())
          return comparator(he1->vertex(), he2->vertex());
      
        if(vertices_he1[0] == vertices_he2[0])
          return comparator(vertices_he1[1], vertices_he2[1]);
        else //either < or >
          return comparator(vertices_he1[0], vertices_he2[0]);
      }
    };


  } //  namespace Mesh_3 {
} // namespace CGAL

#endif COMPARISON_OPERATORS_H
