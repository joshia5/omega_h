#ifndef OMEGA_H_CURVE_COARSEN_HPP
#define OMEGA_H_CURVE_COARSEN_HPP

#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_atomics.hpp"

#include <iostream>
#include <fstream>
namespace Omega_h {

static Few<Real, 10> const BlendedTriangleGetValues(
    Vector<3> const xi, LO b) {
  Few<Real, 10> values;
  Real const blendingTol = 1.e-12;
  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

  for (LO i = 0; i < 3; ++i)
    values[i] = -std::pow(xii[i], b);
  LO const n = 10;//10 nodes per tri
  for (LO i = 3; i < n; ++i)
    values[i] = 0.0;


  const LO nE = 2;//2 nodes per edge
  Few<LO, 6> tev;
  tev[0] = 0;
  tev[1] = 1;
  tev[2] = 1;
  tev[3] = 2;
  tev[4] = 2;
  tev[5] = 0;

  for (LO i = 0; i < 3; ++i) {
    Real x, xiix, xv;
    Few<Real, 4> v;
    x = xii[tev[i*2 + 0]] + xii[tev[i*2 + 1]];

    if (x < blendingTol)
      xiix = 0.5;
    else
      xiix = xii[tev[i*2 +1]]/x;

    xv = xiix;
    v[0] = Bi(3, 0, xv);
    v[1] = Bi(3, 3, xv);
    v[2] = Bi(3, 1, xv);
    v[3] = Bi(3, 2, xv);

    for (LO j = 0; j < 2; ++j)
      values[tev[i*2 + j]] += v[j]*std::pow(x, b);
    for (LO j = 0; j < nE; ++j)
      values[3+i*nE+j] = v[2+j]*std::pow(x, b);
  }
  return values;
}

template <Int dim>
void swap_curved_verts_and_edges(Mesh *mesh, Mesh *new_mesh, const LOs old2new,
    const LOs prods2new, const LOs keys2prods, const LOs keys2edges) {
  auto const nold_verts = mesh->nverts();
  auto const nold_edges = mesh->nedges();
  auto const old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto const old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto const old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto const old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto const old_fv2v = mesh->ask_down(2, 0).ab2b;
  auto const old_v2vf = mesh->ask_up(0,2).a2ab;
  auto const old_vf2f = mesh->ask_up(0,2).ab2b;
  auto const nkeys = keys2prods.size()-1; 
  OMEGA_H_CHECK(nkeys == keys2edges.size());
  // keys2prods is the offset array from the CSR storing keys2new (prods2new is the values array)
  fprintf(stderr, "nedges %d, nkeys %d\n", nold_edges, nkeys);
  if (!mesh->has_tag(0, "bezier_pts")) {
    mesh->add_tag<Real>(0, "bezier_pts", dim, mesh->coords());
  }
  auto const old_vertCtrlPts = mesh->get_ctrlPts(0);
  new_mesh->add_tag<Real>(0, "bezier_pts", dim);
  new_mesh->set_tag_for_ctrlPts(0, old_vertCtrlPts);

  auto const old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto const old_coords = mesh->coords();

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_edge = new_mesh->nedges();
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto const nnew_verts = new_mesh->nverts();

  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim, INT8_MAX);
  Write<I8> edge_crv2bdry_dim(nnew_edge, -1);
  Write<I8> edge_dualCone(nnew_edge, -1);

  //copy ctrl pts for edges
  auto copy_sameedgePts = OMEGA_H_LAMBDA(LO i) {
    if (old2new[i] != -1) {
      LO new_edge = old2new[i];
      for (I8 d = 0; d < dim; ++d) {
        edge_ctrlPts[new_edge*n_edge_pts*dim + d] = 
          old_edgeCtrlPts[i*n_edge_pts*dim + d];
        edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = 
          old_edgeCtrlPts[i*n_edge_pts*dim + dim + d];
      }
    }
  };
  parallel_for(nold_edges, std::move(copy_sameedgePts),
      "copy_same_edgectrlPts");
  auto prod_edge_points = OMEGA_H_LAMBDA(LO i) {
    LO e = prods2new[i];
    auto const v0 = new_ev2v[e*2 + 0];
    auto const v1 = new_ev2v[e*2 + 1];
    for (LO j=0; j<dim; ++j) {
      edge_ctrlPts[e*n_edge_pts*dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])/3.0;
      edge_ctrlPts[e*n_edge_pts*dim + dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*(2.0/3.0);
    }
  };
  parallel_for(prods2new.size(), std::move(prod_edge_points),
      "prod_edge_points");

  auto const newedge_gdim = new_mesh->get_array<I8>(1, "class_dim");
  auto const newedge_gid = new_mesh->get_array<LO>(1, "class_id");
  auto const oldvert_gdim = mesh->get_array<I8>(0, "class_dim");
  auto const oldvert_gid = mesh->get_array<LO>(0, "class_id");
  auto const oldface_gdim = mesh->get_array<I8>(2, "class_dim");
  auto const oldedge_gdim = mesh->get_array<I8>(1, "class_dim");
  auto const oldface_gid = mesh->get_array<LO>(2, "class_id");

  auto const v2v_old = mesh->ask_star(0);
  auto const v2vv_old = v2v_old.a2ab;
  auto const vv2v_old = v2v_old.ab2b;
  auto const v2v_new = new_mesh->ask_star(0);
  auto const v2vv_new = v2v_new.a2ab;
  auto const vv2v_new = v2v_new.ab2b;
  auto const old_v2e = mesh->ask_up(0, 1);
  auto const old_v2ve = old_v2e.a2ab;
  auto const old_ve2e = old_v2e.ab2b;

  Write<LO> nedge_shared_gface_w(nnew_edge, 0);
  Write<Real> cav_edge_len_w(nnew_edge, 0.0);

  // for every old edge, calc and store 2 the unit tangent vectors from each vertex
  //fprintf(stderr, "line 184\n");
  Write<Real> tangents(nold_edges*2*dim, 0.0);
  auto calc_tangents = OMEGA_H_LAMBDA (LO e) {
    LO const e_v0 = old_ev2v[e*2 + 0];
    LO const e_v1 = old_ev2v[e*2 + 1];
    Real length1 = 0.0;
    Real length2 = 0.0;
    for (LO d = 0; d < dim; ++d) {
      tangents[e*2*dim + d] = old_edgeCtrlPts[e*2*dim + d] - 
        old_vertCtrlPts[e_v0*dim + d];
      tangents[e*2*dim + dim + d] = old_edgeCtrlPts[e*2*dim + dim + d] - 
        old_vertCtrlPts[e_v1*dim + d];
      length1 += tangents[e*2*dim + d] * tangents[e*2*dim + d];
      length2 += tangents[e*2*dim + dim + d] * tangents[e*2*dim + dim + d];
    }
    for (LO d = 0; d < dim; ++d) {
      tangents[e*2*dim + d] = tangents[e*2*dim + d]/std::sqrt(length1);
      tangents[e*2*dim + dim + d] = tangents[e*2*dim + dim + d]/std::sqrt(length2);
    }
  };
  parallel_for(nold_edges, std::move(calc_tangents));
  //fprintf(stderr, "line 206\n");

  Write<LO> count_dualCone_cavities(1, 0);
  //Write<LO> count_interior_dualCone_cavities(1, 0); // dont need
  // all keys are in interior
  Write<LO> count_interior_cavities(1, nkeys);
  auto count_dualCone_cav = OMEGA_H_LAMBDA(LO i) {
    /*
    if () {
      atomic_increment(&count_interior_cavities[0]);
    }
    */
    if ((keys2prods[i+1] - keys2prods[i]) == 1) {
      atomic_increment(&count_dualCone_cavities[0]);
      /*
      if () {
        atomic_increment(&count_interior_dualCone_cavities[0]);
      }
      */
    }
  };
  parallel_for(nkeys, std::move(count_dualCone_cav));
  printf("nkeys with dual cone cavities %d\n", count_dualCone_cavities[0]);

  if (dim == 2) {
    auto v2t = mesh->ask_up(0, dim);
    auto e2t = mesh->ask_up(1, dim);
    // here 't' is for triangle 
    LO const t2e_degree = element_degree(OMEGA_H_SIMPLEX, dim, 1); // 3 edges that bound a tri
    auto v2vt = v2t.a2ab;
    auto vt2t = v2t.ab2b;
    auto e2et = e2t.a2ab;
    auto et2t = e2t.ab2b;
    auto v2e = mesh->ask_up(0, 1);
    auto v2ve = v2e.a2ab;
    auto ve2e = v2e.ab2b;
    auto te2e = mesh->ask_down(dim, 1).ab2b;
 
    //fprintf(stderr, "ok 1120\n");
    auto curve_dualCone_cav = OMEGA_H_LAMBDA (LO i) {
      LO const nprods = keys2prods[i+1] - keys2prods[i];
      OMEGA_H_CHECK(nprods == 1);
      LO const e_key = keys2edges[i];
      for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
        LO const new_edge = prods2new[prod];
        //interior (always for 2d swap)
        //if (checks interior swap) {
          LO const new_edge_v0 = new_ev2v[new_edge*2 + 0];
          LO const new_edge_v1 = new_ev2v[new_edge*2 + 1];
          LO const new_edge_v0_old = new_edge_v0;
          LO const new_edge_v1_old = new_edge_v1;
          LO const v_onto = new_edge_v0_old; // just pick one and assign this for now
          auto new_edge_v0_c = get_vector<dim>(new_coords, new_edge_v0);
          auto new_edge_v1_c = get_vector<dim>(new_coords, new_edge_v1);

          Real new_length = 0.;
          for (LO d = 0; d < dim; ++d) {
            new_length += std::pow((new_edge_v1_c[d] - new_edge_v0_c[d]), 2);
          } 
          new_length = std::sqrt(new_length);

          //find lower vtx & check first vertex is vlower or vupper
          LO v_onto_is_first = -1;
          LO v_lower = -1;
          if (v_onto == new_edge_v0_old) {
            v_onto_is_first = 1;
            v_lower = new_edge_v1_old;
          }
          else {
            OMEGA_H_CHECK (v_onto == new_edge_v1_old);
            v_lower = new_edge_v0_old;
          }
          LO e_lower = -1;
          //e_lower is lower half of collapsing edge i.e. edge connecting vkey
          //to vlower;
          /* not needing for 2d so commenting out routine to fine v_lower and
           * t_e_lower */

          //fprintf(stderr, "ok 1173\n");
          Few<LO, 128> lower_edges; // 128 is considered max number of edges defining cones
          Few<I8, 128> from_first_vtx_low;
          LO count_lower_edge = 0;
          for (LO et = e2et[e_key]; et < e2et[e_key + 1]; ++et) {
            //adj tris of ekey
            LO adj_t = et2t[et];
            for (LO te = 0; te < t2e_degree; ++te) {
              LO adj_t_e = te2e[adj_t*t2e_degree + te];
              //adj edges of tri

              LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
              LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
              //adj verts of edge

              if ((adj_t_e_v0 == v_lower) && (adj_t_e != e_key)) {
                OMEGA_H_CHECK(count_lower_edge < 128);
                LO is_duplicate = -1;
                for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                  if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_low[count_lower_edge] = 1;
                  ++count_lower_edge;
                }
              }
              if ((adj_t_e_v1 == v_lower) && (adj_t_e != e_key)) {
                LO is_duplicate = -1;
                for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                  if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_low[count_lower_edge] = -1;
                  ++count_lower_edge;
                }
              }
            }
          }
          //fprintf(stderr, "ok 1213\n");

          Few<Real, dim> t_avg_l;
          for (LO d = 0; d < dim; ++d) t_avg_l[d] = 0.0; 
          for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
            if (from_first_vtx_low[lower_e] == 1) {
              for (LO d = 0; d < dim; ++d) {
                t_avg_l[d] += tangents[lower_edges[lower_e]*2*dim + d];
              }
            }
            if (from_first_vtx_low[lower_e] == -1) {
              for (LO d = 0; d < dim; ++d) {
                t_avg_l[d] += tangents[lower_edges[lower_e]*2*dim + dim + d];
              }
            }
          }
          for (LO d = 0; d < dim; ++d) t_avg_l[d]=t_avg_l[d]/count_lower_edge;
          Real length_t_l = 0.0;
          for (LO d = 0; d < dim; ++d) length_t_l += t_avg_l[d]*t_avg_l[d]; 
          for (LO d = 0; d < dim; ++d) 
            t_avg_l[d] = t_avg_l[d]/std::sqrt(length_t_l);
          Vector<dim> c_lower;
          for (LO d = 0; d < dim; ++d) {
            c_lower[d] = old_coords[v_lower*dim + d] + 
              t_avg_l[d]*new_length/3.0;
          }

          // only using this makes very curved and skinny tets so
          // average t_avg_l and t_e_lower

          Vector<dim> c_upper;
          //init c_upper for inner edges as mid pt of clower and vonto making
          //that 2/3 of the bezier curve straight sided
          for (LO d = 0; d < dim; ++d) {
            //c_upper[d] = (old_coords[v_onto*dim+d] + old_coords[v_key*dim+d])*0.5;
            c_upper[d] = old_coords[v_onto*dim + d] +
              //(c_lower[d] - old_coords[v_onto*dim + d])*0.5;
              (old_coords[v_lower*dim + d] - old_coords[v_onto*dim + d])/3.0;
          }

          //fprintf(stderr, "ok 1268\n");
          //###dual cone
          if (nprods == 1) {
            edge_dualCone[new_edge] = 1;
            //count upper edges
            Few<LO, 128> upper_edges;
            Few<LO, 128> vtx_ring;
            Few<I8, 128> from_first_vtx;
            LO count_upper_edge = 0;
            for (LO et = e2et[e_key]; et < e2et[e_key + 1]; ++et) {
              //adj tris of ekey
              LO adj_t = et2t[et];
              for (LO te = 0; te < t2e_degree; ++te) {
                LO adj_t_e = te2e[adj_t*t2e_degree + te];
                //adj edges of tri

                LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
                LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
                //adj verts of edge

                if ((adj_t_e_v0 == v_onto) && (adj_t_e != e_key)) {
                  OMEGA_H_CHECK(count_upper_edge < 128);
                  LO is_duplicate = -1;
                  for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                    if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                  }
                  if (is_duplicate == -1) {
                    upper_edges[count_upper_edge] = adj_t_e;
                    from_first_vtx[count_upper_edge] = 1;
                    vtx_ring[count_upper_edge] = adj_t_e_v1;
                    ++count_upper_edge;
                  }
                }
                if ((adj_t_e_v1 == v_onto) && (adj_t_e != e_key)) {
                  LO is_duplicate = -1;
                  for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                    if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                  }
                  if (is_duplicate == -1) {
                    upper_edges[count_upper_edge] = adj_t_e;
                    from_first_vtx[count_upper_edge] = -1;
                    vtx_ring[count_upper_edge] = adj_t_e_v0;
                    ++count_upper_edge;
                  }
                }
              }
            }

            Few<Real, dim> t_avg;
            for (LO d = 0; d < dim; ++d) t_avg[d] = 0.0; 
            for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
              if (from_first_vtx[upper_e] == 1) {
                for (LO d = 0; d < dim; ++d) {
                  t_avg[d] += tangents[upper_edges[upper_e]*2*dim + d];
                }
              }
              if (from_first_vtx[upper_e] == -1) {
                for (LO d = 0; d < dim; ++d) {
                  t_avg[d] += tangents[upper_edges[upper_e]*2*dim + dim + d];
                }
              }
            }
            for (LO d = 0; d < dim; ++d) t_avg[d] = t_avg[d]/count_upper_edge;
            Real length_t = 0.0;
            for (LO d = 0; d < dim; ++d) length_t += t_avg[d]*t_avg[d]; 
            for (LO d = 0; d < dim; ++d) t_avg[d] = t_avg[d]/std::sqrt(length_t);

            for (LO d = 0; d < dim; ++d) {
              c_upper[d] = old_coords[v_onto*dim + d] + t_avg[d]*new_length/3.0;
            }
          }
          for (LO d = 0; d < dim; ++d) {
            if (v_onto_is_first == 1) {
              edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_upper[d];
              edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_lower[d];
            }
            else {
              OMEGA_H_CHECK (v_onto_is_first == -1);
              edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_lower[d];
              edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_upper[d];
            }
          }
        //}//interior if
      }
    };
    parallel_for(nkeys, std::move(curve_dualCone_cav), "curve_dualCone_cav");
  }
 
  if (dim == 3) {
    auto v2t = mesh->ask_up(0, dim);
    auto e2t = mesh->ask_up(1, dim);
    // here 't' is for tet 
    LO const t2e_degree = element_degree(OMEGA_H_SIMPLEX, dim, 1); // 6 edges that bound a tet
    auto v2vt = v2t.a2ab;
    auto vt2t = v2t.ab2b;
    auto e2et = e2t.a2ab;
    auto et2t = e2t.ab2b;
    auto v2e = mesh->ask_up(0, 1);
    auto v2ve = v2e.a2ab;
    auto ve2e = v2e.ab2b;
    auto te2e = mesh->ask_down(dim, 1).ab2b;
 
    //fprintf(stderr, "ok 1120\n");
    auto curve_dualCone_cav = OMEGA_H_LAMBDA (LO i) {
      LO const nprods = keys2prods[i+1] - keys2prods[i]; // can be >1
      LO const e_key = keys2edges[i];
      for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
        LO const new_edge = prods2new[prod];
        //interior (only allowing interior swaps right now)
        //if (checks interior swap) {
          LO const new_edge_v0 = new_ev2v[new_edge*2 + 0];
          LO const new_edge_v1 = new_ev2v[new_edge*2 + 1];
          LO const new_edge_v0_old = new_edge_v0;
          LO const new_edge_v1_old = new_edge_v1;
          LO const v_onto = new_edge_v0_old; // just pick one and assign this for now
          auto new_edge_v0_c = get_vector<dim>(new_coords, new_edge_v0);
          auto new_edge_v1_c = get_vector<dim>(new_coords, new_edge_v1);
          
          Real new_length = 0.;
          for (LO d = 0; d < dim; ++d) {
            new_length += std::pow((new_edge_v1_c[d] - new_edge_v0_c[d]), 2);
          } 
          new_length = std::sqrt(new_length);

          //find lower vtx & check first vertex is vlower or vupper
          LO v_onto_is_first = -1;
          LO v_lower = -1;
          if (v_onto == new_edge_v0_old) {
            v_onto_is_first = 1;
            v_lower = new_edge_v1_old;
          }
          else {
            OMEGA_H_CHECK (v_onto == new_edge_v1_old);
            v_lower = new_edge_v0_old;
          }
          LO e_lower = -1;
          //e_lower is lower half of collapsing edge i.e. edge connecting vkey
          //to vlower;
          /* not needing for 2d so commenting out
          for (LO ve = old_v2ve[v_key]; ve < old_v2ve[v_key + 1]; ++ve) {
            LO const e = old_ve2e[ve];
            if ((v_lower == old_ev2v[e*2+0]) || 
                (v_lower == old_ev2v[e*2+1])) e_lower = e;
          }
          Vector<dim> t_e_lower;//tangent unit vec from v_lower
          if (v_lower == old_ev2v[e_lower*2+0]) {
            for (LO d=0; d<dim; ++d) {
              t_e_lower[d] = tangents[e_lower*2*dim + d];
            }
          }
          else {
            OMEGA_H_CHECK (v_lower == old_ev2v[e_lower*2+1]);
            for (LO d=0; d<dim; ++d) {
              t_e_lower[d] = tangents[e_lower*2*dim + dim + d];
            }
          }
          */

          //fprintf(stderr, "ok 1173\n");
          Few<LO, 128> lower_edges; // 128 is considered max number of edges defining cones
          Few<I8, 128> from_first_vtx_low;
          LO count_lower_edge = 0;
          for (LO et = e2et[e_key]; et < e2et[e_key + 1]; ++et) {
            //adj tris of ekey
            LO adj_t = et2t[et];
            for (LO te = 0; te < t2e_degree; ++te) {
              LO adj_t_e = te2e[adj_t*t2e_degree + te];
              //adj edges of tri

              LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
              LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
              //adj verts of edge

              if ((adj_t_e_v0 == v_lower) && (adj_t_e != e_key)) {
                OMEGA_H_CHECK(count_lower_edge < 128);
                LO is_duplicate = -1;
                for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                  if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_low[count_lower_edge] = 1;
                  ++count_lower_edge;
                }
              }
              if ((adj_t_e_v1 == v_lower) && (adj_t_e != e_key)) {
                LO is_duplicate = -1;
                for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                  if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_low[count_lower_edge] = -1;
                  ++count_lower_edge;
                }
              }
            }
          }
          //fprintf(stderr, "ok 1213\n");

          Few<Real, dim> t_avg_l;
          for (LO d = 0; d < dim; ++d) t_avg_l[d] = 0.0; 
          for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
            if (from_first_vtx_low[lower_e] == 1) {
              for (LO d = 0; d < dim; ++d) {
                t_avg_l[d] += tangents[lower_edges[lower_e]*2*dim + d];
              }
            }
            if (from_first_vtx_low[lower_e] == -1) {
              for (LO d = 0; d < dim; ++d) {
                t_avg_l[d] += tangents[lower_edges[lower_e]*2*dim + dim + d];
              }
            }
          }
          for (LO d = 0; d < dim; ++d) t_avg_l[d]=t_avg_l[d]/count_lower_edge;
          Real length_t_l = 0.0;
          for (LO d = 0; d < dim; ++d) length_t_l += t_avg_l[d]*t_avg_l[d]; 
          for (LO d = 0; d < dim; ++d) 
            t_avg_l[d] = t_avg_l[d]/std::sqrt(length_t_l);
          Vector<dim> c_lower;
          for (LO d = 0; d < dim; ++d) {
            c_lower[d] = old_coords[v_lower*dim + d] + 
              t_avg_l[d]*new_length/3.0;
          }

          // only using this makes very curved and skinny tets so
          // average t_avg_l and t_e_lower
          /*
          // wont need this for swap, also theres no elower for swap
          if (nprods > 1) {
            for (LO d = 0; d < dim; ++d) t_avg_l[d] = (t_avg_l[d]+t_e_lower[d])*0.5;
            length_t_l = 0.0;
            for (LO d = 0; d < dim; ++d) length_t_l += t_avg_l[d]*t_avg_l[d]; 
            for (LO d = 0; d < dim; ++d)
              t_avg_l[d] = t_avg_l[d]/std::sqrt(length_t_l);

            for (LO d = 0; d < dim; ++d) {
              //c_lower[d] = (old_coords[v_lower*dim+d] + old_coords[v_key*dim+d])*0.5;
              c_lower[d] = old_coords[v_lower*dim + d] + 
              //t_avg_l[d]*new_length/3.0;
              (old_coords[v_onto*dim + d] - old_coords[v_lower*dim + d])/3.0;
            }
          }
          */

          Vector<dim> c_upper;
          //init c_upper for inner edges as mid pt of clower and vonto making
          //that 2/3 of the bezier curve straight sided
          for (LO d = 0; d < dim; ++d) {
            //c_upper[d] = (old_coords[v_onto*dim+d] + old_coords[v_key*dim+d])*0.5;
            c_upper[d] = old_coords[v_onto*dim + d] +
              //(c_lower[d] - old_coords[v_onto*dim + d])*0.5;
              (old_coords[v_lower*dim + d] - old_coords[v_onto*dim + d])/3.0;
          }

          //fprintf(stderr, "ok 1268\n");
          //###dual cone
          if (nprods == 1) {
            edge_dualCone[new_edge] = 1;
            //count upper edges
            Few<LO, 128> upper_edges;
            Few<LO, 128> vtx_ring;
            Few<I8, 128> from_first_vtx;
            LO count_upper_edge = 0;
            for (LO et = e2et[e_key]; et < e2et[e_key + 1]; ++et) {
              //adj tris of ekey
              LO adj_t = et2t[et];
              for (LO te = 0; te < t2e_degree; ++te) {
                LO adj_t_e = te2e[adj_t*t2e_degree + te];
                //adj edges of tri

                LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
                LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
                //adj verts of edge

                if ((adj_t_e_v0 == v_onto) && (adj_t_e != e_key)) {
                  OMEGA_H_CHECK(count_upper_edge < 128);
                  LO is_duplicate = -1;
                  for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                    if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                  }
                  if (is_duplicate == -1) {
                    upper_edges[count_upper_edge] = adj_t_e;
                    from_first_vtx[count_upper_edge] = 1;
                    vtx_ring[count_upper_edge] = adj_t_e_v1;
                    ++count_upper_edge;
                  }
                }
                if ((adj_t_e_v1 == v_onto) && (adj_t_e != e_key)) {
                  LO is_duplicate = -1;
                  for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                    if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                  }
                  if (is_duplicate == -1) {
                    upper_edges[count_upper_edge] = adj_t_e;
                    from_first_vtx[count_upper_edge] = -1;
                    vtx_ring[count_upper_edge] = adj_t_e_v0;
                    ++count_upper_edge;
                  }
                }
              }
            }

            Few<Real, dim> t_avg;
            for (LO d = 0; d < dim; ++d) t_avg[d] = 0.0; 
            for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
              if (from_first_vtx[upper_e] == 1) {
                for (LO d = 0; d < dim; ++d) {
                  t_avg[d] += tangents[upper_edges[upper_e]*2*dim + d];
                }
              }
              if (from_first_vtx[upper_e] == -1) {
                for (LO d = 0; d < dim; ++d) {
                  t_avg[d] += tangents[upper_edges[upper_e]*2*dim + dim + d];
                }
              }
            }
            for (LO d = 0; d < dim; ++d) t_avg[d] = t_avg[d]/count_upper_edge;
            Real length_t = 0.0;
            for (LO d = 0; d < dim; ++d) length_t += t_avg[d]*t_avg[d]; 
            for (LO d = 0; d < dim; ++d) t_avg[d] = t_avg[d]/std::sqrt(length_t);

            for (LO d = 0; d < dim; ++d) {
              c_upper[d] = old_coords[v_onto*dim + d] + t_avg[d]*new_length/3.0;
            }
          }
          for (LO d = 0; d < dim; ++d) {
            if (v_onto_is_first == 1) {
              edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_upper[d];
              edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_lower[d];
            }
            else {
              OMEGA_H_CHECK (v_onto_is_first == -1);
              edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_lower[d];
              edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_upper[d];
            }
          }
        //}//interior if
      }
    };
    parallel_for(nkeys, std::move(curve_dualCone_cav), "curve_dualCone_cav");
  }

  /*
   * commenting this out for now because we have interior cavities as
   * per Omega_h_swap.cpp line 22
  Write<LO> count_bdry_cavities(1, 0);
  if (dim == 3) {
    auto calc_gface_prods = OMEGA_H_LAMBDA(LO i) {
      for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
        LO const new_edge = prods2new[prod];
        // new edge on model face
        if (newedge_gdim[new_edge] == 2) {
          nedge_shared_gface_w[new_edge] += 1;
          for (LO prod2 = keys2prods[i]; prod2 < keys2prods[i+1]; ++prod2) {
            LO const other_edge = prods2new[prod2];
            if ((other_edge != new_edge) && (newedge_gdim[other_edge] == 2) &&
                (newedge_gid[new_edge] == newedge_gid[other_edge])) {
              //calc new edge length
              LO const new_edge_v0 = new_ev2v[new_edge*2 + 0];
              LO const new_edge_v1 = new_ev2v[new_edge*2 + 1];
              auto const v0 = get_vector<dim>(new_coords, new_edge_v0);
              auto const v1 = get_vector<dim>(new_coords, new_edge_v1);
              Real length = 0.0;
              for (LO d=0; d<dim; ++d) {
                length += (v1[d] - v0[d])*(v1[d] - v0[d]);
              }
              length = std::sqrt(length);
              cav_edge_len_w[new_edge] += length;
              nedge_shared_gface_w[new_edge] += 1;
            }
          }
          cav_edge_len_w[new_edge] = cav_edge_len_w[new_edge]/
            (nedge_shared_gface_w[new_edge]*1.0);
        }
      }
    };
    parallel_for(nkeys, std::move(calc_gface_prods));
  }
  auto nedge_shared_gface = Read<LO>(nedge_shared_gface_w);
  auto cav_edge_len = Read<Real>(cav_edge_len_w);

  auto count_bdry_cavs = OMEGA_H_LAMBDA(LO i) {
    LO const v_key = keys2verts[i];
    LO const v_onto = keys2verts_onto[i];
    LO has_bdry_cav = -1;
    for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
      LO new_edge = prods2new[prod];
      //edges classified on g_edges or g_faces
      // oldvert classified on g_verts or g-edges
      if (((oldvert_gdim[v_key] <= 1) && (oldvert_gdim[v_onto] <= 1) &&
          (newedge_gdim[new_edge] == 1)) || 
         ((dim == 3) && (newedge_gdim[new_edge] == 2))) {
          has_bdry_cav = 1;
        }
    }
    if (has_bdry_cav == 1) 
      atomic_increment(&count_bdry_cavities[0]);
  };
  parallel_for(nkeys, std::move(count_bdry_cavs));

  auto curve_bdry_edges = OMEGA_H_LAMBDA(LO i) {
    LO const v_key = keys2verts[i];
    LO const v_onto = keys2verts_onto[i];
    Few<I8, 256> bdry_cands_avail;
    Few<Real, 256> all_cands_dist;//max 16 cand and 16 prods
    for (LO k=0; k<256; ++k) all_cands_dist[k] = DBL_MAX;
    for (LO a=0; a<256; ++a) bdry_cands_avail[a] = 1;

    for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
      LO const new_edge = prods2new[prod];
      LO const nedge_shared_gface_i = nedge_shared_gface[new_edge];
      LO const new_edge_v0 = new_ev2v[new_edge*2 + 0];
      LO const new_edge_v1 = new_ev2v[new_edge*2 + 1];
      LO const new_edge_v0_old = same_verts2old_verts[ab2b[a2ab[new_edge_v0]]];
      LO const new_edge_v1_old = same_verts2old_verts[ab2b[a2ab[new_edge_v1]]];
      auto const c0 = get_vector<dim>(old_vertCtrlPts, new_edge_v0);
      auto const c0_coord = get_vector<dim>(old_coords, new_edge_v0_old);
      auto const c3 = get_vector<dim>(old_vertCtrlPts, new_edge_v1);
      auto const c3_coord = get_vector<dim>(old_coords, new_edge_v1_old);
      auto new_length = 
        (c3_coord[0] - c0_coord[0])*(c3_coord[0] - c0_coord[0]) + 
        (c3_coord[1] - c0_coord[1])*(c3_coord[1] - c0_coord[1]) + 
        (c3_coord[2] - c0_coord[2])*(c3_coord[2] - c0_coord[2]);
      new_length = std::sqrt(new_length);
      Vector<dim> c1,c2;

      //check for new edge, first vertex is vlower or vupper
      LO v_onto_is_first = -1;
      LO v_lower = -1;
      if (v_onto == new_edge_v0_old) {
        v_onto_is_first = 1;
        v_lower = new_edge_v1_old;
      }
      else {
        OMEGA_H_CHECK (v_onto == new_edge_v1_old);
        v_lower = new_edge_v0_old;
      }

      //edges classified on g_edges
      if ((oldvert_gdim[v_key] <= 1) && (oldvert_gdim[v_onto] <= 1) &&
          (newedge_gdim[new_edge] == 1)) {
	edge_crv2bdry_dim[new_edge] = 1;
	auto old_p1 = get_vector<dim>(old_vertCtrlPts, v_key);
	Real xi_1 = 0.5;
        
	Real sum_dist1 = 0.0;
	Real sum_dist2 = 0.0;
	for (Int j = 0; j < dim; ++j) {
	  sum_dist1 += (old_p1[j] - c0[j])*(old_p1[j] - c0[j]);
	  sum_dist2 += (old_p1[j] - c3[j])*(old_p1[j] - c3[j]);
	}
	sum_dist1 = std::sqrt(sum_dist1/(dim*1.0));
	sum_dist2 = std::sqrt(sum_dist2/(dim*1.0));
	xi_1 = sum_dist1/(sum_dist1 + sum_dist2);
        
	Vector<dim> old_c1;
	for (Int j = 0; j < dim; ++j) {
	  old_c1[j] = (old_p1[j] - B0_quad(xi_1)*c0[j] - 
                       B2_quad(xi_1)*c3[j])/B1_quad(xi_1);
	}
	for (LO d = 0; d < dim; ++d) {
	  c1[d] = (1.0/3.0)*c0[d] + (2.0/3.0)*old_c1[d];
	  c2[d] = (2.0/3.0)*old_c1[d] + (1.0/3.0)*c3[d];
          edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c1[d];
          edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c2[d];
	}
      }
      //end:edges classified on g_edges

      //edges classified on g_faces
      if ((dim == 3) && (newedge_gdim[new_edge] == 2)) {
	edge_crv2bdry_dim[new_edge] = 2;
        //init with straight sided
        for (LO d = 0; d < dim; ++d) {
          c1[d] = edge_ctrlPts[new_edge*n_edge_pts*dim + d];
          c2[d] = edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d];
	}

        //find lower vtx, here upper is v_onto and lower is other end of edge
        Vector<dim> c_lower;
        LO e_lower = -1;
        LO e_upper = -1;
        LO e_upper_from_v_onto = -1;
        LO e_lower_from_v_lower= -1;
        //e_lower is lower half of collapsing edge i.e. edge connecting vkey
        //to vlower;
        //e_upper is upper half of collapsing edge i.e. edge connecting vkey
        //to vonto;
        for (LO ve = old_v2ve[v_key]; ve < old_v2ve[v_key + 1]; ++ve) {
          LO const e = old_ve2e[ve];
          if ((v_lower == old_ev2v[e*2+0]) || 
              (v_lower == old_ev2v[e*2+1])) e_lower = e;
          if ((v_onto == old_ev2v[e*2+0]) || 
              (v_onto == old_ev2v[e*2+1])) e_upper = e;
          if (v_onto == old_ev2v[e*2+0]) e_upper_from_v_onto = 1;
          if (v_lower== old_ev2v[e*2+0]) e_lower_from_v_lower= 1;
        }
        Vector<dim> t_e_lower;//tangent unit vec from either one end-vtx
        Vector<dim> t_e_upper;//tangent unit vec from vonto
        for (LO d=0; d<dim; ++d) {
          if (e_lower_from_v_lower> 0) {
            t_e_lower[d] = tangents[e_lower*2*dim + d];
          }
          else {
            t_e_lower[d] = tangents[e_lower*2*dim + dim + d];
          }
          if (e_upper_from_v_onto > 0) {
            t_e_upper[d] = tangents[e_upper*2*dim + d];
          }
          else {
            t_e_upper[d] = tangents[e_upper*2*dim + dim + d];
          }
        }

        //count upper edges
        Few<LO, 2> upper_edges;
        Few<I8, 2> from_first_vtx;
        LO count_upper_edge = 0;
        for (LO vf = old_v2vf[v_key]; vf < old_v2vf[v_key + 1]; ++vf) {
          //adj tris of vkey
          LO const adj_t = old_vf2f[vf];
          for (LO te = 0; te < 3; ++te) {
            LO adj_t_e = old_fe2e[adj_t*3 + te];
            //adj edges of tri
            LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
            LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
            //adj verts of edge
            if ((adj_t_e_v0 == v_onto) && (adj_t_e_v1 != v_key)) {
              LO is_duplicate = -1;
              for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
              }
              if (is_duplicate == -1) {
                if ((oldface_gid[adj_t] == newedge_gid[new_edge]) && (oldface_gdim[adj_t] == 2)) {
                  OMEGA_H_CHECK(count_upper_edge < 2);
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = 1;
                  ++count_upper_edge;
                }
              }
            }
            if ((adj_t_e_v1 == v_onto) && (adj_t_e_v0 != v_key)) {
              LO is_duplicate = -1;
              for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
              }
              if (is_duplicate == -1) {
                if ((oldface_gid[adj_t] == newedge_gid[new_edge]) && (oldface_gdim[adj_t] == 2)) {
                  OMEGA_H_CHECK(count_upper_edge < 2);
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = -1;
                  ++count_upper_edge;
                }
              }
            }
            if (oldvert_gdim[v_key] == 1) {
              if ((adj_t_e_v0 == v_onto) && (adj_t_e_v1 == v_key)) {
                LO is_duplicate = -1;
                for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                  if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  OMEGA_H_CHECK(count_upper_edge < 2);
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = 1;
                  ++count_upper_edge;
                }
              }
              if ((adj_t_e_v1 == v_onto) && (adj_t_e_v0 == v_key)) {
                LO is_duplicate = -1;
                for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                  if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  OMEGA_H_CHECK(count_upper_edge < 2);
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = -1;
                  ++count_upper_edge;
                }
              }
            }
          }
        }

        Vector<dim> c_upper;
        Vector<dim> c_avg_upper;
        Few<Real, 2*dim> upper_tangents;
        Vector<dim> t_avg_upper;
        for (LO d = 0; d < dim; ++d) t_avg_upper[d] = 0.0; 
        for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
          if (from_first_vtx[upper_e] == 1) {
            for (LO d = 0; d < dim; ++d) {
              t_avg_upper[d] += tangents[upper_edges[upper_e]*2*dim + d];
              upper_tangents[upper_e*dim + d] = tangents[upper_edges[upper_e]*2*dim + d];
            }
          }
          if (from_first_vtx[upper_e] == -1) {
            for (LO d = 0; d < dim; ++d) {
              t_avg_upper[d] += tangents[upper_edges[upper_e]*2*dim + dim + d];
              upper_tangents[upper_e*dim + d] = tangents[upper_edges[upper_e]*2*dim + dim + d];
            }
          }
        }

        //calc c using upper avg dir to check concavity
        Real length_t = 0.0;
        LO concave_upper = -1;
        for (LO d = 0; d < dim; ++d) 
          t_avg_upper[d] = t_avg_upper[d]/count_upper_edge;
        for (LO d = 0; d < dim; ++d) 
          length_t += t_avg_upper[d]*t_avg_upper[d]; 
        for (LO d = 0; d < dim; ++d) 
          t_avg_upper[d] = t_avg_upper[d]/std::sqrt(length_t);
        for (LO d = 0; d < dim; ++d) {
          c_avg_upper[d] = old_coords[v_onto*dim + d] + 
                           t_avg_upper[d]*new_length/3.0;
        }
        //calc upper concavity
        Real dist_to_lower = 0.0;
        for (LO d = 0; d < dim; ++d) {
          dist_to_lower += std::pow(
              (c_avg_upper[d] - old_coords[v_lower*dim + d]), 2);
        }
        dist_to_lower = std::sqrt(dist_to_lower);
        if (dist_to_lower > new_length) {
          concave_upper = 1;
          for (LO d = 0; d < dim; ++d) {
            c_avg_upper[d] = 
              old_coords[v_onto*dim + d] - t_avg_upper[d]*new_length/3.0;
          }
        }
        //

        if (nedge_shared_gface_i == 1) {
          for (LO d = 0; d < dim; ++d) c_upper[d] = c_avg_upper[d];
        }
        else {
          assert(nedge_shared_gface_i > 1);
          assert(nedge_shared_gface_i <= 32);
          Few<Real, 3*32> cand_tangents;
          Few<Real, 3*32> cand_c;

          Few<Real, 32> cand_angle_to_uppere0;
          Few<LO, 32> sorted_cands;
          Few<LO, 32> cand_ids;
          for (LO cand = 0; cand < 32; ++cand) {
            sorted_cands[cand] = -1;
            cand_ids[cand] = cand;
            cand_angle_to_uppere0[cand] = DBL_MAX;
          }
          LO const uppere0 = upper_edges[0];
          Vector<dim> uppere0_pt;
          LO const uppere0_v0 = old_ev2v[uppere0*2 + 0];
          LO const uppere0_v1 = old_ev2v[uppere0*2 + 1];
          LO uppere0_vlower = -1;
          if (uppere0_v0 == v_onto) {
            uppere0_vlower = uppere0_v1;
          }
          if (uppere0_v1 == v_onto) {
            uppere0_vlower = uppere0_v0;
          }
          Real uppere0_len = 0.0;
          for (LO d=0; d<dim; ++d) {
            uppere0_len += std::pow((old_coords[uppere0_v0*dim + d] -
                                     old_coords[uppere0_v1*dim + d]), 2);
          }
          uppere0_len = std::sqrt(uppere0_len);
          for (LO d=0; d<dim; ++d) {
            uppere0_pt[d] = old_coords[v_onto*dim + d] +
              (old_coords[uppere0_vlower*dim + d] - old_coords[v_onto*dim + d])/2.0;
          }

          //calc n locations of new tangent pts
          Real upper_theta = acos(
              upper_tangents[0]*upper_tangents[dim + 0] +
              upper_tangents[1]*upper_tangents[dim + 1] +
              upper_tangents[2]*upper_tangents[dim + 2]);
          Vector<dim> n;
          //cross
          n[0] = (upper_tangents[1]*upper_tangents[dim + 2]) -
                 (upper_tangents[2]*upper_tangents[dim + 1]);
          n[1] =-(upper_tangents[0]*upper_tangents[dim + 2]) +
                 (upper_tangents[2]*upper_tangents[dim + 0]);
          n[2] = (upper_tangents[0]*upper_tangents[dim + 1]) -
                 (upper_tangents[1]*upper_tangents[dim + 0]);
          Real length_n = 0.0;
          for (LO d = 0; d < dim; ++d) length_n += n[d]*n[d];
          for (LO d = 0; d < dim; ++d) n[d] = n[d]/std::sqrt(length_n);
          for (LO cand = 0; cand < nedge_shared_gface_i; ++cand) {
            Real theta_c = (cand+1)*upper_theta/(nedge_shared_gface_i + 1);
            Vector<3> b;
            b[0] = 0.0; b[1] = cos(theta_c); b[2] = cos(upper_theta-theta_c);
            auto A = tensor_3(
                n[0], n[1], n[2],
                upper_tangents[0], upper_tangents[1], upper_tangents[2],
                upper_tangents[dim+0], upper_tangents[dim+1], upper_tangents[dim+2]
                );
            auto A_inv = invert(A);
            auto X = A_inv*b;

            if (std::abs(upper_theta - PI) < 1e-3) {

              n[0] = (upper_tangents[1]*t_e_upper[2]) - 
                (upper_tangents[2]*t_e_upper[1]);
              n[1] =-(upper_tangents[0]*t_e_upper[2]) +
                (upper_tangents[2]*t_e_upper[0]);
              n[2] = (upper_tangents[0]*t_e_upper[1]) - 
                (upper_tangents[1]*t_e_upper[0]);
              Real len_n =0.0;
              for (LO d=0;d<3;++d) len_n+= std::pow(n[d],2);
              len_n = std::sqrt(len_n);
              for (LO d=0;d<3;++d) n[d] = n[d]/len_n;
              if (n[0] > 0.99) {
                n[0] = 1.0000;
                n[1] = 0.0000;
                n[2] = 0.0000;
              }
              if (n[0] < -0.99) {
                n[0] = -1.0000;
                n[1] = 0.0000;
                n[2] = 0.0000;
              }
              if (n[1] > 0.99) {
                n[0] = 0.0000;
                n[1] = 1.0000;
                n[2] = 0.0000;
              }
              if (n[1] < -0.99) {
                n[0] = 0.0000;
                n[1] = -1.0000;
                n[2] = 0.0000;
              }
              if (n[2] > 0.99) {
                n[0] = 0.0000;
                n[1] = 0.0000;
                n[2] = 1.0000;
              }
              if (n[2] < -0.99) {
                n[0] = 0.0000;
                n[1] = 0.0000;
                n[2] = -1.0000;
              }

              //hand-calc +/- y axis as tang if A is ill-cond
              if ((std::abs(n[0])<1e-3) && ((std::abs(n[1])-1.0)<1e-3)
                  && (std::abs(n[2])<1e-3)) {
                X[0] = upper_tangents[0]*cos(theta_c) +
                       upper_tangents[2]*sin(theta_c);
                X[1] = upper_tangents[1];
                X[2] =-upper_tangents[0]*sin(theta_c) + 
                       upper_tangents[2]*cos(theta_c);

                //calc temp candidate to check if this solved value of tangent
                //vector is pointing in right direction or should be flipped
                Vector<dim> temp_c, temp_t;
                length_t = 0.0;
                for (LO d=0; d<dim; ++d) length_t += std::pow(X[d], 2);
                for (LO d=0; d<dim; ++d) temp_t[d] = X[d]/std::sqrt(length_t);
                for (LO d=0; d<dim; ++d) {
                  temp_c[d] = old_coords[v_onto*dim + d] +
                    temp_t[d]*new_length/3.0;
                }
                Real temp_dist = 0.0;//dist from tempc to vlower
                for (LO d = 0; d < dim; ++d) {
                  temp_dist+=std::pow((temp_c[d]-old_coords[v_lower*dim+d]),2);
                }
                temp_dist = std::sqrt(temp_dist);
                if (temp_dist > new_length) {
                  X[0] = upper_tangents[0]*cos(-theta_c) +
                    upper_tangents[2]*sin(-theta_c);
                  X[2] =-upper_tangents[0]*sin(-theta_c) + 
                    upper_tangents[2]*cos(-theta_c);
                }
              }
              //hand-calc +/- x axis as tang if A is ill-cond
              if ((std::abs(n[1])<1e-3) && ((std::abs(n[0])-1.0)<1e-3)
                  && (std::abs(n[2])<1e-3)) {
                //made by straight edge with uppere0
                X[0] = upper_tangents[0];
                X[1] = upper_tangents[1]*cos(theta_c) -
                       upper_tangents[2]*sin(theta_c);
                X[2] = upper_tangents[1]*sin(theta_c) + 
                       upper_tangents[2]*cos(theta_c);

                //calc temp candidate to check if this solved value of tangent
                //vector is pointing in right direction or should be flipped
                Vector<dim> temp_c, temp_t;
                length_t = 0.0;
                for (LO d=0; d<dim; ++d) length_t += std::pow(X[d], 2);
                for (LO d=0; d<dim; ++d) temp_t[d] = X[d]/std::sqrt(length_t);
                for (LO d=0; d<dim; ++d) {
                  temp_c[d] = old_coords[v_onto*dim + d] +
                    temp_t[d]*new_length/3.0;
                }
                Real temp_dist = 0.0;//dist from tempc to vlower
                for (LO d = 0; d < dim; ++d) {
                  temp_dist+=std::pow((temp_c[d]-old_coords[v_lower*dim+d]),2);
                }
                temp_dist = std::sqrt(temp_dist);
                if (temp_dist > new_length) {
                  X[1] = upper_tangents[1]*cos(-theta_c) -
                    upper_tangents[2]*sin(-theta_c);
                  X[2] = upper_tangents[1]*sin(-theta_c) + 
                    upper_tangents[2]*cos(-theta_c);
                }
              }
              //hand-calc +/- z axis as tang if A is ill-cond
              if ((std::abs(n[0])<1e-3) && ((std::abs(n[2])-1.0)<1e-3)
                  && (std::abs(n[1])<1e-3)) {
                //made by straight edge with uppere0
                X[0] = upper_tangents[0]*cos(theta_c) -
                       upper_tangents[1]*sin(theta_c);
                X[1] = upper_tangents[0]*sin(theta_c) +
                       upper_tangents[1]*cos(theta_c);
                X[2] = upper_tangents[2];

                //calc temp candidate to check if this solved value of tangent
                //vector is pointing in right direction or should be flipped
                Vector<dim> temp_c, temp_t;
                length_t = 0.0;
                for (LO d=0; d<dim; ++d) length_t += std::pow(X[d], 2);
                for (LO d=0; d<dim; ++d) temp_t[d] = X[d]/std::sqrt(length_t);
                for (LO d=0; d<dim; ++d) {
                  temp_c[d] = old_coords[v_onto*dim + d] +
                    temp_t[d]*new_length/3.0;
                }
                Real temp_dist = 0.0;//dist from tempc to vlower
                for (LO d = 0; d < dim; ++d) {
                  temp_dist+=std::pow((temp_c[d]-old_coords[v_lower*dim+d]),2);
                }
                temp_dist = std::sqrt(temp_dist);
                if (temp_dist > new_length) {
                  X[0] = upper_tangents[0]*cos(-theta_c) -
                         upper_tangents[1]*sin(-theta_c);
                  X[1] = upper_tangents[0]*sin(-theta_c) +
                         upper_tangents[1]*cos(-theta_c);
                }
              }
            }

            length_t = 0.0;
            for (LO d=0; d<dim; ++d) {
              cand_tangents[cand*dim + d] = X[d];
              length_t += std::pow(cand_tangents[cand*dim + d], 2);
            }
            for (LO d = 0; d < dim; ++d) {
              cand_tangents[cand*dim + d] = cand_tangents[cand*dim + d]/
                std::sqrt(length_t);
            }
            for (LO d = 0; d < dim; ++d) {
              cand_c[cand*dim + d] = old_coords[v_onto*dim + d] +
                cand_tangents[cand*dim + d]*new_length/3.0;
                //v_onto is always upper
            }
            if (concave_upper == 1) {
              Vector<dim> cand_vec;
              Real cand_len = 0.0;
              for (LO d = 0; d < dim; ++d) {
                //flip upper tangent
                cand_c[cand*dim + d] = old_coords[v_onto*dim + d] -
                  cand_tangents[cand*dim + d]*new_length/3.0;
                cand_vec[d] = cand_c[cand*dim+d] - old_coords[v_onto*dim+d];
                cand_len += std::pow(cand_vec[d], 2);
              }
              cand_len = std::sqrt(cand_len);
              for (LO d = 0; d < dim; ++d) {
                cand_vec[d] = cand_vec[d]/cand_len;
              }
              cand_angle_to_uppere0[cand] = acos(
                  upper_tangents[0]*cand_vec[0] +
                  upper_tangents[1]*cand_vec[1] +
                  upper_tangents[2]*cand_vec[2]);
            }
          }
          
          //unless concave sorting cands is not needed for cands
          if (concave_upper == 1) {
            //sort cands
            for (LO count_c2 = 0; count_c2 < nedge_shared_gface_i; ++count_c2) {
              for (LO count_cand2_2 = count_c2; count_cand2_2 < nedge_shared_gface_i; ++count_cand2_2) {
                if (cand_angle_to_uppere0[count_c2] < cand_angle_to_uppere0[count_cand2_2]) {
                  sorted_cands[count_c2] = cand_ids[count_c2];
                }
                else {
                  sorted_cands[count_c2] = cand_ids[count_cand2_2];
                  swap2(cand_angle_to_uppere0[count_c2], cand_angle_to_uppere0[count_cand2_2]);
                  swap2(cand_ids[count_c2], cand_ids[count_cand2_2]);
                }
              }
            }
          }

          //###FOR PRODS
          //find dist of all relevant prods to uppere0
          Few<Real, 32> prod_angle_to_uppere1;
          Few<LO, 32> sorted_prods;
          Few<LO, 32> prod_ids;
          for (LO count_p2 = 0; count_p2 < 32; ++count_p2) {
            sorted_prods[count_p2] = -1;
            prod_ids[count_p2] = -1;
            prod_angle_to_uppere1[count_p2] = DBL_MAX;
          }
          LO count_prod2 = 0;
          //0. find unit vec for e0 as stored earlir in tags
          //1. we are finding otherp so find unit vec in dir of other p
          //2. take dot product 
          //3. store acos angle value as angle_touppere0
          //4. sort
          for (LO prod2 = keys2prods[i]; prod2 < keys2prods[i+1]; ++prod2) {
            LO const other_edge = prods2new[prod2];
            if ((newedge_gdim[other_edge] == 2) &&
                (newedge_gid[new_edge] == newedge_gid[other_edge])) {
              LO const oth_edge_v0 = new_ev2v[other_edge*2 + 0];
              LO const oth_edge_v1 = new_ev2v[other_edge*2 + 1];
              LO const oth_edge_v0_old = same_verts2old_verts[ab2b[a2ab[oth_edge_v0]]];
              LO const oth_edge_v1_old = same_verts2old_verts[ab2b[a2ab[oth_edge_v1]]];
              auto const c0_coord2= get_vector<dim>(old_coords, oth_edge_v0_old);
              auto const c3_coord2= get_vector<dim>(old_coords, oth_edge_v1_old);
              Vector<dim> other_p;
              Vector<dim> other_vec;
              Real other_len = 0.0;

              //check for new edge, first vertex is vlower or vupper
              if (v_onto == oth_edge_v0_old) {
                for (LO d = 0; d < dim; ++d) {
                  other_p[d] = c0_coord2[d] + (1.0/3.0)*(c3_coord2[d] - c0_coord2[d]);
                  other_vec[d] = other_p[d] - c0_coord2[d];
                  other_len += std::pow(other_vec[d], 2);
                }
              }
              else {
                OMEGA_H_CHECK(v_onto == oth_edge_v1_old);
                for (LO d = 0; d < dim; ++d) {
                  other_p[d] = c3_coord2[d] + (1.0/3.0)*(c0_coord2[d] - c3_coord2[d]);
                  other_vec[d] = other_p[d] - c3_coord2[d];
                  other_len += std::pow(other_vec[d], 2);
                }
              }
              for (LO d = 0; d < dim; ++d) {
                other_vec[d] = other_vec[d]/std::sqrt(other_len);
              }
              //found unit vec in dir of other edge
              //now find angle

              prod_angle_to_uppere1[count_prod2] = acos(
                 other_vec[0]*upper_tangents[dim+0] + 
                 other_vec[1]*upper_tangents[dim+1] + 
                 other_vec[2]*upper_tangents[dim+2]); 
              prod_ids[count_prod2] = other_edge;

              ++count_prod2;
            }
          }
          OMEGA_H_CHECK(count_prod2 == nedge_shared_gface_i);
          //sort prods by angle to upper e0
          for (LO count_p2 = 0; count_p2 < nedge_shared_gface_i; ++count_p2) {
            for (LO count_prod2_2 = count_p2; count_prod2_2 < nedge_shared_gface_i; ++count_prod2_2) {
              if (prod_angle_to_uppere1[count_p2] > prod_angle_to_uppere1[count_prod2_2]) {
                sorted_prods[count_p2] = prod_ids[count_p2];
              }
              else {
                sorted_prods[count_p2] = prod_ids[count_prod2_2];
                swap2(prod_angle_to_uppere1[count_p2], prod_angle_to_uppere1[count_prod2_2]);
                swap2(prod_ids[count_p2], prod_ids[count_prod2_2]);
              }
            }
          }
          //assign optimal cand with optimal prod
          LO opt_cand_id = -1;
          for (LO count_p2 = 0; count_p2 < nedge_shared_gface_i; ++count_p2) {
            if (prod_ids[count_p2] == new_edge) opt_cand_id = count_p2;
          }
          if (concave_upper == 1) {
            opt_cand_id = sorted_cands[opt_cand_id];
          }
          for (LO d = 0; d < dim; ++d) {
            c_upper[d] = cand_c[opt_cand_id*dim + d];
          }
        }

        Few<LO, 2> lower_edges;
        Few<I8, 2> from_first_vtx_l;
        LO count_lower_edge = 0;
        for (LO vf = old_v2vf[v_key]; vf < old_v2vf[v_key + 1]; ++vf) {
          //adj tris of vkey
          LO const adj_t = old_vf2f[vf];
          for (LO te = 0; te < 3; ++te) {
            LO adj_t_e = old_fe2e[adj_t*3 + te];
            //adj edges of tri
            LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
            LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
            //adj verts of edge
            if ((adj_t_e_v0 == v_lower) && (adj_t_e_v1 != v_key)) {
              LO is_duplicate = -1;
              for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
              }
              if (is_duplicate == -1) {
                if ((oldface_gid[adj_t] == newedge_gid[new_edge]) && (oldface_gdim[adj_t] == 2)) {
                  OMEGA_H_CHECK(count_lower_edge < 2);
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_l[count_lower_edge] = 1;
                  ++count_lower_edge;
                }
              }
            }
            if ((adj_t_e_v1 == v_lower) && (adj_t_e_v0 != v_key)) {
              LO is_duplicate = -1;
              for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
              }
              if (is_duplicate == -1) {
                if ((oldface_gid[adj_t] == newedge_gid[new_edge]) && (oldface_gdim[adj_t] == 2)) {
                  OMEGA_H_CHECK(count_lower_edge < 2);
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_l[count_lower_edge] = -1;
                  ++count_lower_edge;
                }
              }
            }
          }
        }

        Vector<dim> t_lower;
        Few<Real, 2*dim> lower_tangents;
        for (LO d = 0; d < dim; ++d) t_lower[d] = 0.0; 
        for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
          if (from_first_vtx_l[lower_e] == 1) {
            for (LO d = 0; d < dim; ++d) {
              t_lower[d] += tangents[lower_edges[lower_e]*2*dim + d];
              lower_tangents[lower_e*dim + d] = tangents[lower_edges[lower_e]*2*dim + d];
            }
          }
          if (from_first_vtx_l[lower_e] == -1) {
            for (LO d = 0; d < dim; ++d) {
              t_lower[d] += tangents[lower_edges[lower_e]*2*dim + dim + d];
              lower_tangents[lower_e*dim + d] = tangents[lower_edges[lower_e]*2*dim + dim + d];
            }
          }
        }

        OMEGA_H_CHECK(count_lower_edge == 2);
        for (LO d = 0; d < dim; ++d) t_lower[d] = t_lower[d]/count_lower_edge;
        length_t = 0.0;

        for (LO d = 0; d < dim; ++d) length_t += t_lower[d]*t_lower[d]; 

        if (std::abs(length_t - 0.0) < 1e-3) {
          Real lower_theta = acos(
            lower_tangents[0]*lower_tangents[dim + 0] +
            lower_tangents[1]*lower_tangents[dim + 1] +
            lower_tangents[2]*lower_tangents[dim + 2]);
          Vector<dim> n;
          //cross
          n[0] = (lower_tangents[1]*lower_tangents[dim + 2]) - 
                 (lower_tangents[2]*lower_tangents[dim + 1]);
          n[1] =-(lower_tangents[0]*lower_tangents[dim + 2]) +
                 (lower_tangents[2]*lower_tangents[dim + 0]);
          n[2] = (lower_tangents[0]*lower_tangents[dim + 1]) - 
                 (lower_tangents[1]*lower_tangents[dim + 0]);

          Real theta_c = lower_theta/2.0;
          Vector<3> b;
          b[0] = 0.0; b[1] = cos(theta_c); b[2] = cos(lower_theta-theta_c);
          auto A = tensor_3(
              n[0], n[1], n[2],
              lower_tangents[0], lower_tangents[1], lower_tangents[2],
              lower_tangents[dim+0], lower_tangents[dim+1], lower_tangents[dim+2]
              );

          auto A_inv = invert(A);
          auto X = A_inv*b;
          //hand-calc +/- x axis as tang if A is ill-cond
          if (std::abs(lower_theta - PI) < 1e-3) {
            n[0] = (lower_tangents[1]*t_e_lower[2]) - 
              (lower_tangents[2]*t_e_lower[1]);
            n[1] =-(lower_tangents[0]*t_e_lower[2]) +
              (lower_tangents[2]*t_e_lower[0]);
            n[2] = (lower_tangents[0]*t_e_lower[1]) - 
              (lower_tangents[1]*t_e_lower[0]);
            Real len_n =0.0;
            for (LO d=0;d<3;++d) len_n+= std::pow(n[d],2);
            for (LO d=0;d<3;++d) n[d] = n[d]/std::sqrt(len_n);//unit normal
            if (n[0] > 0.99) {
              n[0] = 1.0000;
              n[1] = 0.0000;
              n[2] = 0.0000;
            }
            if (n[0] < -0.99) {
              n[0] = -1.0000;
              n[1] = 0.0000;
              n[2] = 0.0000;
            }
            if (n[1] > 0.99) {
              n[0] = 0.0000;
              n[1] = 1.0000;
              n[2] = 0.0000;
            }
            if (n[1] < -0.99) {
              n[0] = 0.0000;
              n[1] = -1.0000;
              n[2] = 0.0000;
            }
            if (n[2] > 0.99) {
              n[0] = 0.0000;
              n[1] = 0.0000;
              n[2] = 1.0000;
            }
            if (n[2] < -0.99) {
              n[0] = 0.0000;
              n[1] = 0.0000;
              n[2] = -1.0000;
            }
            //+x
            if (((std::abs(n[0])-1.0)<1e-3) && (std::abs(n[1])<1e-3)
                && (std::abs(n[2])<1e-3)) {
              X[0] = lower_tangents[0];
              X[1] = lower_tangents[1]*cos(theta_c) -
                     lower_tangents[2]*sin(theta_c);
              X[2] = lower_tangents[1]*sin(theta_c) +
                     lower_tangents[2]*cos(theta_c);
  
              //calc temp candidate to check if this solved value of tangent
              //vector is pointing in right direction or should be flipped
              Vector<dim> temp_c, temp_t;
              Real length_t2 = 0.0;
              for (LO d=0; d<dim; ++d) length_t2 += std::pow(X[d], 2);
              for (LO d=0; d<dim; ++d) temp_t[d] = X[d]/std::sqrt(length_t2);
              for (LO d=0; d<dim; ++d) {
                temp_c[d] = old_coords[v_lower*dim + d] +
                  temp_t[d]*new_length/3.0;
              }
              Real temp_dist = 0.0;//dist from tempc to vupper
              for (LO d = 0; d < dim; ++d) {
                temp_dist+=std::pow((temp_c[d]-old_coords[v_onto*dim+d]),2);
              }
              temp_dist = std::sqrt(temp_dist);
              if (temp_dist > new_length) {
              X[1] = lower_tangents[1]*cos(-theta_c) -
                     lower_tangents[2]*sin(-theta_c);
              X[2] = lower_tangents[1]*sin(-theta_c) +
                     lower_tangents[2]*cos(-theta_c);
              }
            }

            //hand-calc +/- y axis as tang if A is ill-cond
            if ((std::abs(n[0])<1e-3) && ((std::abs(n[1])-1.0)<1e-3)
             && (std::abs(n[2])<1e-3)) {
              X[0] = lower_tangents[0]*cos(theta_c) +
                     lower_tangents[2]*sin(theta_c);
              X[1] = lower_tangents[1];
              X[2] =-lower_tangents[0]*sin(theta_c) +
                     lower_tangents[2]*cos(theta_c);

              //calc temp candidate to check if this solved value of tangent
              //vector is pointing in right direction or should be flipped
              Vector<dim> temp_c, temp_t;
              Real length_t2 = 0.0;
              for (LO d=0; d<dim; ++d) length_t2 += std::pow(X[d], 2);
              for (LO d=0; d<dim; ++d) temp_t[d] = X[d]/std::sqrt(length_t2);
              for (LO d=0; d<dim; ++d) {
                temp_c[d] = old_coords[v_lower*dim + d] +
                  temp_t[d]*new_length/3.0;
              }
              Real temp_dist = 0.0;//dist from tempc to vupper
              for (LO d = 0; d < dim; ++d) {
                temp_dist+=std::pow((temp_c[d]-old_coords[v_onto*dim+d]),2);
              }
              temp_dist = std::sqrt(temp_dist);
              if (temp_dist > new_length) {
              X[0] = lower_tangents[0]*cos(-theta_c) +
                     lower_tangents[2]*sin(-theta_c);
              X[2] =-lower_tangents[0]*sin(-theta_c) +
                     lower_tangents[2]*cos(-theta_c);
              }
            }

            //hand-calc +/- z axis as tang if A is ill-cond
            if ((std::abs(n[0])<1e-3) && ((std::abs(n[2])-1.0)<1e-3)
             && (std::abs(n[1])<1e-3)) {
              X[0] = lower_tangents[0]*cos(theta_c) -
                     lower_tangents[1]*sin(theta_c);
              X[1] = lower_tangents[0]*sin(theta_c) +
                     lower_tangents[1]*cos(theta_c);
              X[2] = lower_tangents[2];

              //calc temp candidate to check if this solved value of tangent
              //vector is pointing in right direction or should be flipped
              Vector<dim> temp_c, temp_t;
              Real length_t2 = 0.0;
              for (LO d=0; d<dim; ++d) length_t2 += std::pow(X[d], 2);
              for (LO d=0; d<dim; ++d) temp_t[d] = X[d]/std::sqrt(length_t2);
              for (LO d=0; d<dim; ++d) {
                temp_c[d] = old_coords[v_lower*dim + d] +
                  temp_t[d]*new_length/3.0;
              }
              Real temp_dist = 0.0;//dist from tempc to vupper
              for (LO d = 0; d < dim; ++d) {
                temp_dist+=std::pow((temp_c[d]-old_coords[v_onto*dim+d]),2);
              }
              temp_dist = std::sqrt(temp_dist);
              if (temp_dist > new_length) {
                X[0] = lower_tangents[0]*cos(-theta_c) -
                       lower_tangents[1]*sin(-theta_c);
                X[1] = lower_tangents[0]*sin(-theta_c) +
                       lower_tangents[1]*cos(-theta_c);
              }
            }
          }
          for (LO d=0; d<dim; ++d) {
            t_lower[d] = X[d];
            length_t += std::pow(X[d], 2);
          }
        }

        for (LO d = 0; d < dim; ++d) t_lower[d] = t_lower[d]/std::sqrt(length_t);
        for (LO d = 0; d < dim; ++d) {
          c_lower[d] = old_coords[v_lower*dim + d] + t_lower[d]*new_length/3.0;
        }

        //For concave lower angle//
        {
          Real dist_to_upper = 0.0;
          for (LO d = 0; d < dim; ++d) {
            dist_to_upper += (c_lower[d] - old_coords[v_onto*dim + d])*
              (c_lower[d] - old_coords[v_onto*dim + d]);
          }
          dist_to_upper = std::sqrt(dist_to_upper);
          if (dist_to_upper > new_length) {
            for (LO d = 0; d < dim; ++d) {
              c_lower[d] = old_coords[v_lower*dim + d] - t_lower[d]*new_length/3.0;
            }
          }
        }
        //

        for (LO d = 0; d < dim; ++d) {
          if (v_onto_is_first == 1) {
            edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_lower[d];
            edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_upper[d];
          }
          else {
            OMEGA_H_CHECK (v_onto_is_first == -1);
            edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_lower[d];
            edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_upper[d];
          }
        }
      } // dim=3 conditional and for edges classified on g_faces
    } // prods loop
  };
  parallel_for(nkeys, curve_bdry_edges);
  */

  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim, Reals(edge_ctrlPts));
  new_mesh->add_tag<I8>(1, "edge_crv2bdry_dim", 1, Read<I8>(edge_crv2bdry_dim));
  new_mesh->add_tag<I8>(1, "edge_dualCone", 1, Read<I8>(edge_dualCone));

  return;
}

template <Int dim>
void swap_curved_faces(Mesh *mesh, Mesh *new_mesh, const LOs old2new,
    const LOs prods2new) {
  auto const new_fv2v = new_mesh->ask_down(2, 0).ab2b;
  auto const new_fe2e = new_mesh->get_adj(2, 1).ab2b;
  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_faces = new_mesh->nfaces();
  auto const nnew_edges = new_mesh->nedges();
  auto const nold_faces = mesh->nfaces();
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto const vertCtrlPts = new_mesh->get_ctrlPts(0);
  auto const edgeCtrlPts = new_mesh->get_ctrlPts(1);
  Write<Real> face_ctrlPts(nnew_faces*1*dim, INT8_MAX);

  fprintf(stderr, "nfaces %d\n", nold_faces);
  //copy ctrl pts for faces
  auto copy_samefacePts = OMEGA_H_LAMBDA(LO i) {
    if (old2new[i] != -1) {
      LO new_face = old2new[i];
      for (I8 d = 0; d < dim; ++d) {
        face_ctrlPts[new_face*dim + d] = old_faceCtrlPts[i*dim + d];
      }
    }
  };
  parallel_for(nold_faces, std::move(copy_samefacePts),
      "copy_same_facectrlPts");

  auto face_centroids = OMEGA_H_LAMBDA(LO i) {
    LO const tri = prods2new[i];
    LO const tri2v_degree = element_degree(OMEGA_H_SIMPLEX, FACE, 1);
    LO const v0 = new_fv2v[tri*tri2v_degree + 0];
    LO const v1 = new_fv2v[tri*tri2v_degree + 1];
    LO const v2 = new_fv2v[tri*tri2v_degree + 2];
    for (LO j = 0; j < dim; ++j) {
      face_ctrlPts[tri*dim + j] = (new_coords[v0*dim + j] +
         new_coords[v1*dim + j] + new_coords[v2*dim + j])/3.0;
    }
  };
  parallel_for(prods2new.size(), std::move(face_centroids), "face_centroids");

  //loop over new edges
  //if new edge dual cone
  //  find adjacent faces in a for loop
  //  pass face ids as input to blended tri fn 
  //  then take first 9 values and mult them with all ctrl pts of face
  //  that gives the interp pt
  //  face interp to ctrl pt
  if (dim == 3) {
    Vector<3> face_xi;
    face_xi[0] = 1.0/3.0;
    face_xi[1] = 1.0/3.0;
    face_xi[2] = 1.0/3.0;
    auto const weights = BlendedTriangleGetValues(face_xi, 2);
    auto edge_dualCone = new_mesh->get_array<I8>(1, "edge_dualCone");
    auto edge_crv2bdry_dim = new_mesh->get_array<I8>(1, "edge_crv2bdry_dim");
    auto const e2et = new_mesh->ask_up(1, 2).a2ab;
    auto const et2t = new_mesh->ask_up(1, 2).ab2b;
    auto const newedge_gdim = new_mesh->get_array<I8>(1, "class_dim");
    auto const newedge_gid = new_mesh->get_array<LO>(1, "class_id");
    auto const newface_gdim = new_mesh->get_array<I8>(2, "class_dim");
    auto const newface_gid = new_mesh->get_array<LO>(2, "class_id");

    auto face_blends = OMEGA_H_LAMBDA(LO e) {
      if (edge_dualCone[e] == 1) {
        for (LO et = e2et[e]; et < e2et[e + 1]; ++et) {
          LO const tri = et2t[et];

          auto p11 = face_blend_interp_3d(3, tri, new_ev2v, new_fe2e,
              vertCtrlPts, edgeCtrlPts, new_fv2v, weights);
          auto newface_c11 = face_interpToCtrlPt_3d(3, tri, new_ev2v, new_fe2e,
              vertCtrlPts, edgeCtrlPts, p11, new_fv2v);
          for (LO k = 0; k < 3; ++k) {
            face_ctrlPts[tri*dim + k] = newface_c11[k];
          }

        }
      }
      if ((dim == 3) && (newedge_gdim[e] == 2)) {
      //if (edge_crv2bdry_dim[e] == 2) 
        for (LO et = e2et[e]; et < e2et[e + 1]; ++et) {
          LO const tri = et2t[et];
          if ((newedge_gid[e] == newface_gid[tri]) && (newface_gdim[tri == 2])) {
            auto p11 = face_blend_interp_3d(3, tri, new_ev2v, new_fe2e,
                vertCtrlPts, edgeCtrlPts, new_fv2v, weights);
            auto newface_c11 = face_interpToCtrlPt_3d(3, tri, new_ev2v, new_fe2e,
                vertCtrlPts, edgeCtrlPts, p11, new_fv2v);
            for (LO k = 0; k < 3; ++k) {
              face_ctrlPts[tri*dim + k] = newface_c11[k];
            }
          }
        }
      }
    };
    parallel_for(nnew_edges, std::move(face_blends), "face_blends");
  }

  new_mesh->add_tag<Real>(FACE, "bezier_pts", dim);
  new_mesh->set_tag_for_ctrlPts(FACE, Reals(face_ctrlPts));

  return;
}

#define OMEGA_H_INST(T)                                                       
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

} // namespace Omega_h

#endif
