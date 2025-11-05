#include "Omega_h_refine.hpp"

#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_indset.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_modify.hpp"
#include "Omega_h_profile.hpp"
#include "Omega_h_refine_qualities.hpp"
#include "Omega_h_refine_topology.hpp"
#include "Omega_h_transfer.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_curve_coarsen.hpp"
#include "Omega_h_curve_validity_3d.hpp"

#include "Omega_h_file.hpp"
#include "Omega_h_build.hpp"

namespace Omega_h {

static bool refine_ghosted(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto edges_are_cands = mesh->get_array<I8>(EDGE, "candidate");
  mesh->remove_tag(EDGE, "candidate");
  auto cands2edges = collect_marked(edges_are_cands);
  auto cand_quals = refine_qualities(mesh, cands2edges);
  auto cands_are_good = each_geq_to(cand_quals, opts.min_quality_allowed);
  if (get_max(comm, cands_are_good) != 1) return false;
  auto nedges = mesh->nedges();
  auto edges_are_initial =
      map_onto(cands_are_good, cands2edges, nedges, I8(0), 1);
  auto edge_quals = map_onto(cand_quals, cands2edges, nedges, 0.0, 1);
  auto edges_are_keys = find_indset(mesh, EDGE, edge_quals, edges_are_initial);
  mesh->add_tag(EDGE, "key", 1, edges_are_keys);
  mesh->add_tag(EDGE, "rep_vertex2md_order", 1,
      get_rep2md_order_adapt(mesh, EDGE, VERT, edges_are_keys));
  auto keys2edges = collect_marked(edges_are_keys);
  Graph edges2elems;
  if (mesh->dim() == 1)
    edges2elems = identity_graph(mesh->nedges());
  else
    edges2elems = mesh->ask_up(EDGE, mesh->dim());
  set_owners_by_indset(mesh, EDGE, keys2edges, edges2elems);
  return true;
}

static bool refine_ghosted_crv(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto edges_are_cands = mesh->get_array<I8>(EDGE, "candidate");
  mesh->remove_tag(EDGE, "candidate");
  auto cands2edges = collect_marked(edges_are_cands);
  auto cand_quals = refine_qualities(mesh, cands2edges);
  auto cands_are_good = each_geq_to(cand_quals, INT64_MIN*1.);
  if (get_max(comm, cands_are_good) != 1) return false;
  auto nedges = mesh->nedges();
  auto edges_are_initial =
      map_onto(cands_are_good, cands2edges, nedges, I8(0), 1);
  auto edge_quals = map_onto(cand_quals, cands2edges, nedges, 0.0, 1);
  auto edges_are_keys = find_indset(mesh, EDGE, edge_quals, edges_are_initial);
  mesh->add_tag(EDGE, "key", 1, edges_are_keys);
  mesh->add_tag(EDGE, "rep_vertex2md_order", 1,
      get_rep2md_order_adapt(mesh, EDGE, VERT, edges_are_keys));
  auto keys2edges = collect_marked(edges_are_keys);
  Graph edges2elems;
  if (mesh->dim() == 1)
    edges2elems = identity_graph(mesh->nedges());
  else
    edges2elems = mesh->ask_up(EDGE, mesh->dim());
  set_owners_by_indset(mesh, EDGE, keys2edges, edges2elems);
  return true;
}

static void refine_element_based(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto edges_are_keys = mesh->get_array<I8>(EDGE, "key");
  auto keys2edges = collect_marked(edges_are_keys);
  auto nkeys = keys2edges.size();
  auto ntotal_keys = comm->allreduce(GO(nkeys), OMEGA_H_SUM);
  if (opts.verbosity >= EACH_REBUILD && comm->rank() == 0) {
    std::cout << "refining " << ntotal_keys << " edges\n";
  }
  auto new_mesh = mesh->copy_meta();
  auto keys2midverts = LOs();
  auto old_verts2new_verts = LOs();
  auto old_lows2new_lows = LOs();

  auto keys2old_faces = LOs();

  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    auto keys2prods = LOs();
    auto prod_verts2verts = LOs();
    if (ent_dim == VERT) {
      keys2prods = LOs(nkeys + 1, 0, 1);
    } else {
      refine_products(mesh, ent_dim, keys2edges, keys2midverts,
          old_verts2new_verts, keys2prods, prod_verts2verts);
    }
    auto prods2new_ents = LOs();
    auto same_ents2old_ents = LOs();
    auto same_ents2new_ents = LOs();
    auto old_ents2new_ents = LOs();
    modify_ents_adapt(mesh, &new_mesh, ent_dim, EDGE, keys2edges, keys2prods,
        prod_verts2verts, old_lows2new_lows, &prods2new_ents,
        &same_ents2old_ents, &same_ents2new_ents, &old_ents2new_ents);
    if (ent_dim == VERT) {
      keys2midverts = prods2new_ents;
      old_verts2new_verts = old_ents2new_ents;
    }
    transfer_refine(mesh, opts.xfer_opts, &new_mesh, keys2edges, keys2midverts,
        ent_dim, keys2prods, prods2new_ents, same_ents2old_ents,
        same_ents2new_ents);

    old_lows2new_lows = old_ents2new_ents;
  }
  *mesh = new_mesh;
}

static void refine_element_based_crv(Mesh* mesh, AdaptOpts const& opts,
    I8 const should_modify_mesh, Read<I8> edges_are_cands) {
  auto comm = mesh->comm();
  auto edges_are_keys = mesh->get_array<I8>(EDGE, "key");
  auto keys2edges = collect_marked(edges_are_keys);
  auto nkeys = keys2edges.size();
  auto ntotal_keys = comm->allreduce(GO(nkeys), OMEGA_H_SUM);
  if (opts.verbosity >= EACH_REBUILD && comm->rank() == 0) {
    std::cout << "refining " << ntotal_keys << " edges\n";
  }
  auto new_mesh = mesh->copy_meta();
  auto keys2midverts = LOs();
  auto old_verts2new_verts = LOs();
  auto old_lows2new_lows = LOs();

  auto keys2old_faces = LOs();

  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    auto keys2prods = LOs();
    auto prod_verts2verts = LOs();
    if (ent_dim == VERT) {
      keys2prods = LOs(nkeys + 1, 0, 1);
    } else {
      refine_products(mesh, ent_dim, keys2edges, keys2midverts,
          old_verts2new_verts, keys2prods, prod_verts2verts);
    }
    auto prods2new_ents = LOs();
    auto same_ents2old_ents = LOs();
    auto same_ents2new_ents = LOs();
    auto old_ents2new_ents = LOs();
    modify_ents_adapt(mesh, &new_mesh, ent_dim, EDGE, keys2edges, keys2prods,
        prod_verts2verts, old_lows2new_lows, &prods2new_ents,
        &same_ents2old_ents, &same_ents2new_ents, &old_ents2new_ents);
    if (ent_dim == VERT) {
      keys2midverts = prods2new_ents;
      old_verts2new_verts = old_ents2new_ents;
    }
    transfer_refine(mesh, opts.xfer_opts, &new_mesh, keys2edges, keys2midverts,
        ent_dim, keys2prods, prods2new_ents, same_ents2old_ents,
        same_ents2new_ents);

    if (ent_dim == EDGE) {
      if (mesh->is_curved() > 0) {
        printf("refining mesh of order %d\n", mesh->get_max_order());
        new_mesh.set_curved(1);
        new_mesh.set_max_order(mesh->get_max_order());
        new_mesh.add_tag<I8>
          (1, "n_bezier_pts", 1, Bytes(new_mesh.nents(1), mesh->get_max_order()-1,
                                       "numBezierPts"));
        if (mesh->get_max_order() == 3) {
          if (mesh->dim() == 2) {
            keys2old_faces = create_curved_verts_and_edges_2d
              (mesh, &new_mesh, old_ents2new_ents, prods2new_ents, keys2prods,
               keys2midverts, old_verts2new_verts, keys2edges);
          }
          else {
            OMEGA_H_CHECK(mesh->dim() == 3);
            keys2old_faces = create_curved_verts_and_edges_3d
              (mesh, &new_mesh, old_ents2new_ents, prods2new_ents, keys2prods,
               keys2midverts, old_verts2new_verts, keys2edges);
          }
        }
        else {
          OMEGA_H_CHECK(mesh->get_max_order() == 2);
          if (mesh->dim() == 2) {
            Omega_h_fail("error aborting\n");
          }
          else {
            OMEGA_H_CHECK(mesh->dim() == 3);
            create_curved_verts_and_edges_3d_p2
              (mesh, &new_mesh, old_ents2new_ents, prods2new_ents, keys2prods,
               keys2midverts, old_verts2new_verts, keys2edges,
               opts.transfer_cands, edges_are_cands);
          }
        }
      }
    }
    if (ent_dim == FACE) {
      if (mesh->is_curved() > 0){
        if (mesh->get_max_order() == 3) {
          if (mesh->dim() == 2) {
            create_curved_faces_2d(mesh, &new_mesh, old_ents2new_ents, prods2new_ents,
                keys2prods, keys2edges, keys2old_faces,
                old_verts2new_verts);
          }
          else {
            OMEGA_H_CHECK (mesh->dim() == 3);
            create_curved_faces_3d(mesh, &new_mesh, old_ents2new_ents, prods2new_ents,
                keys2prods, keys2edges, keys2old_faces,
                old_verts2new_verts);
          }
        }
      }
    }
    if (ent_dim == REGION) {
      if (should_modify_mesh < 0) {
        /* from coarsen
            auto edge_cand_codes = get_edge_codes(mesh);
            auto edges_are_cands = each_neq_to(edge_cand_codes, I8(DONT_COLLAPSE));
            auto cands2edges = collect_marked(edges_are_cands);
            auto cand_edge_codes = read(unmap(cands2edges, edge_cand_codes, 1));
            // this does not require ghosting as 
            // ghosting called later for ind. set selection
            auto cand_edge_invalidities = coarsen_invalidities_new_mesh
              (mesh, cands2edges, cand_edge_codes, &new_mesh,
               old_verts2new_verts, verts_are_keys, keys2verts_onto,
               prods2new_ents);
            cand_edge_codes = filter_coarsen_invalids(
                cand_edge_codes, cand_edge_invalidities, -1);
            filter_coarsen_candidates(&cands2edges, &cand_edge_codes);
            if (mesh->is_curved() > 0) {
              put_edge_codes(mesh, cands2edges, cand_edge_codes);
            }
            */
      }
      else {
        if (opts.check_crv_qual > 0) {
          printf("checking validity after refine\n");
          check_validity_all_tet(&new_mesh);
          auto quals = calc_crvQuality_3d(&new_mesh);
        }
      }
    }

    old_lows2new_lows = old_ents2new_ents;
  }

  /*
  if (mesh->is_curved() > 0) {
    if (opts.check_crv_qual > 0) {
      printf("checking validity after refine\n");
      check_validity_all_tet(&new_mesh);
      printf("calc quality after refine\n");
      auto qual = calc_crvQuality_3d(&new_mesh);
    }
  }
  */
  if (should_modify_mesh > 0) {
    *mesh = new_mesh;
    if (opts.verbosity >= WRITE_FILE && comm->rank() == 0) {
      printf("writing refine mesh\n");
      auto wireframe_mesh = Mesh(comm->library());
      wireframe_mesh.set_comm(comm);
      build_cubic_wireframe_3d(mesh, &wireframe_mesh, 4);
      std::string vtuPath =
        "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/refine_wire.vtu";
      vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
      auto cubic_curveVtk_mesh = Mesh(comm->library());
      cubic_curveVtk_mesh.set_comm(comm);
      build_cubic_curveVtk_3d(mesh, &cubic_curveVtk_mesh, 4);
      vtuPath =
        "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/refine.vtu";
      vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
    }
  }
  return;
}

bool refine(Mesh* mesh, AdaptOpts const& opts) {

  auto edges_are_cands = mesh->get_array<I8>(EDGE, "candidate");
  mesh->set_parting(OMEGA_H_GHOSTED);

  if (mesh->is_curved() < 0) { // linear mesh
    if (!refine_ghosted(mesh, opts)) return false;
    mesh->set_parting(OMEGA_H_ELEM_BASED);
    refine_element_based(mesh, opts);
  }
  if (mesh->is_curved() > 0) {
    if (!refine_ghosted_crv(mesh, opts)) return false;
    mesh->set_parting(OMEGA_H_ELEM_BASED);
    if (opts.min_crv_qual_allowed > 0.) {
      refine_element_based_crv(mesh, opts, 1, edges_are_cands);
    }
  }

  return true;
}

bool refine_by_size(Mesh* mesh, AdaptOpts const& opts) {
  OMEGA_H_TIME_FUNCTION;
  auto comm = mesh->comm();
  auto lengths = mesh->ask_lengths();
  auto edge_is_cand = each_gt(lengths, opts.max_length_desired);
  if (get_max(comm, edge_is_cand) != 1) return false;
  mesh->add_tag(EDGE, "candidate", 1, edge_is_cand);
  return refine(mesh, opts);
}

}  // end namespace Omega_h
