#include "Omega_h_build.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_box.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_inertia.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_migrate.hpp"
#include "Omega_h_owners.hpp"
#include "Omega_h_simplify.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_collapse.hpp"

#include <iostream>
#include <fstream>

namespace Omega_h {

void add_ents2verts(
    Mesh* mesh, Int ent_dim, LOs ev2v, GOs vert_globals, GOs elem_globals) {
  auto comm = mesh->comm();
  auto nverts_per_ent = element_degree(mesh->family(), ent_dim, VERT);
  auto ne = divide_no_remainder(ev2v.size(), nverts_per_ent);
  Remotes owners;
  if (comm->size() > 1) {
    if (mesh->could_be_shared(ent_dim)) {
      if (ent_dim == mesh->dim()) {
        owners = owners_from_globals(comm, elem_globals, Read<I32>());
      } else {
        resolve_derived_copies(
            comm, vert_globals, nverts_per_ent, &ev2v, &owners);
      }
    } else {
      owners = identity_remotes(comm, ne);
    }
  }
  if (ent_dim == 1) {
    mesh->set_ents(ent_dim, Adj(ev2v));
  } else {
    auto ldim = ent_dim - 1;
    auto lv2v = mesh->ask_verts_of(ldim);
    auto v2l = mesh->ask_up(VERT, ldim);
    auto down = reflect_down(ev2v, lv2v, v2l, mesh->family(), ent_dim, ldim);
    mesh->set_ents(ent_dim, down);
  }
  if (comm->size() > 1) {
    mesh->set_owners(ent_dim, owners);
    if (ent_dim == mesh->dim() && elem_globals.exists()) {
      mesh->add_tag(ent_dim, "global", 1, elem_globals);
    } else {
      mesh->add_tag(ent_dim, "global", 1, globals_from_owners(mesh, ent_dim));
    }
  } else {
    mesh->add_tag(ent_dim, "global", 1, GOs(ne, 0, 1));
  }
}

void build_verts_from_globals(Mesh* mesh, GOs vert_globals) {
  auto comm = mesh->comm();
  auto nverts = vert_globals.size();
  mesh->set_verts(nverts);
  mesh->add_tag(VERT, "global", 1, vert_globals);
  if (comm->size() > 1) {
    mesh->set_owners(
        VERT, owners_from_globals(comm, vert_globals, Read<I32>()));
  }
}

void build_ents_from_elems2verts(
    Mesh* mesh, LOs ev2v, GOs vert_globals, GOs elem_globals) {
  auto comm = mesh->comm();
  auto elem_dim = mesh->dim();
  for (Int mdim = 1; mdim < elem_dim; ++mdim) {
    auto mv2v = find_unique(ev2v, mesh->family(), elem_dim, mdim);
    add_ents2verts(mesh, mdim, mv2v, vert_globals, elem_globals);
  }
  add_ents2verts(mesh, elem_dim, ev2v, vert_globals, elem_globals);
  if (!comm->reduce_and(is_sorted(vert_globals))) {
    reorder_by_globals(mesh);
  }
}

void build_from_elems2verts(Mesh* mesh, CommPtr comm, Omega_h_Family family,
    Int edim, LOs ev2v, Read<GO> vert_globals) {
  mesh->set_comm(comm);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->set_family(family);
  mesh->set_dim(edim);
  build_verts_from_globals(mesh, vert_globals);
  build_ents_from_elems2verts(mesh, ev2v, vert_globals);
}

void build_from_elems2verts(
    Mesh* mesh, Omega_h_Family family, Int edim, LOs ev2v, LO nverts) {
  auto vert_globals = Read<GO>(nverts, 0, 1);
  build_from_elems2verts(
      mesh, mesh->library()->self(), family, edim, ev2v, vert_globals);
}

void build_from_elems_and_coords(
    Mesh* mesh, Omega_h_Family family, Int edim, LOs ev2v, Reals coords) {
  auto nverts = coords.size() / edim;
  build_from_elems2verts(mesh, family, edim, ev2v, nverts);
  mesh->add_coords(coords);
}

void build_box_internal(Mesh* mesh, Omega_h_Family family, Real x, Real y,
    Real z, LO nx, LO ny, LO nz, bool symmetric) {
  OMEGA_H_CHECK(nx > 0);
  OMEGA_H_CHECK(ny >= 0);
  OMEGA_H_CHECK(nz >= 0);
  if (ny == 0) {
    LOs ev2v;
    Reals coords;
    make_1d_box(x, nx, &ev2v, &coords);
    build_from_elems_and_coords(mesh, family, EDGE, ev2v, coords);
  } else if (nz == 0) {
    LOs fv2v;
    Reals coords;
    make_2d_box(x, y, nx, ny, &fv2v, &coords);
    if (family == OMEGA_H_SIMPLEX && (!symmetric)) fv2v = tris_from_quads(fv2v);
    auto fam2 = symmetric ? OMEGA_H_HYPERCUBE : family;
    build_from_elems_and_coords(mesh, fam2, FACE, fv2v, coords);
    if (family == OMEGA_H_SIMPLEX && symmetric) tris_from_quads_symmetric(mesh);
  } else {
    LOs rv2v;
    Reals coords;
    make_3d_box(x, y, z, nx, ny, nz, &rv2v, &coords);
    if (family == OMEGA_H_SIMPLEX && (!symmetric)) rv2v = tets_from_hexes(rv2v);
    auto fam2 = symmetric ? OMEGA_H_HYPERCUBE : family;
    build_from_elems_and_coords(mesh, fam2, REGION, rv2v, coords);
    if (family == OMEGA_H_SIMPLEX && symmetric) tets_from_hexes_symmetric(mesh);
  }
}

Mesh build_box(CommPtr comm, Omega_h_Family family, Real x, Real y, Real z,
    LO nx, LO ny, LO nz, bool symmetric) {
  auto lib = comm->library();
  auto mesh = Mesh(lib);
  if (comm->rank() == 0) {
    build_box_internal(&mesh, family, x, y, z, nx, ny, nz, symmetric);
    reorder_by_hilbert(&mesh);
    classify_box(&mesh, x, y, z, nx, ny, nz);
    mesh.class_sets = get_box_class_sets(mesh.dim());
  }
  mesh.set_comm(comm);
  mesh.balance();
  return mesh;
}

/* When we try to build a mesh from _partitioned_
   element-to-vertex connectivity only, we have to derive
   consistent edges and faces in parallel.
   We'll start by using the usual local derivation on each
   part, which creates the right entities but with inconsistent
   alignment and no ownership information.
   This function establishes the ownership of derived entities
   and canonicalizes their connectivity.
   It does this by expressing connectivity in terms of vertex
   global numbers, locally sorting, and sending each entity
   to its lowest-global-number vertex.
   It uses a linear partitioning of the vertices by global number,
   and each "server" vertex handles the entity copies for which
   it is the lowest-global-number vertex.
   We reuse the matching code that powers reflect_down() to identify
   copies which have the exact same connectivity and establish ownership.
*/

LOs sort_locally_based_on_rank(
    LOs servers_to_served, Read<I32> served_to_rank) {
  OMEGA_H_TIME_FUNCTION;
  auto const served_order = Write<LO>(served_to_rank.size());
  auto functor = OMEGA_H_LAMBDA(LO const server) {
    auto const begin = servers_to_served[server];
    auto const end = servers_to_served[server + 1];
    I32 last_smallest_rank = -1;
    for (LO i = begin; i < end;) {
      I32 next_smallest_rank = ArithTraits<I32>::max();
      for (LO j = begin; j < end; ++j) {
        auto const rank = served_to_rank[j];
        if (rank > last_smallest_rank && rank < next_smallest_rank) {
          next_smallest_rank = rank;
        }
      }
      OMEGA_H_CHECK(next_smallest_rank > last_smallest_rank);
      OMEGA_H_CHECK(next_smallest_rank < ArithTraits<I32>::max());
      for (LO j = begin; j < end; ++j) {
        if (served_to_rank[j] == next_smallest_rank) {
          served_order[i++] = j;
        }
      }
      last_smallest_rank = next_smallest_rank;
    }
  };
  parallel_for(servers_to_served.size() - 1, std::move(functor));
  return served_order;
}

void resolve_derived_copies(CommPtr comm, Read<GO> verts2globs, Int deg,
    LOs* p_ent_verts2verts, Remotes* p_ents2owners) {
  // entity vertices to vertices
  auto const ev2v = *p_ent_verts2verts;
  // entity vertices to vertex globals
  auto const ev2vg = read(unmap(ev2v, verts2globs, 1));
  auto const canon_codes = get_codes_to_canonical(deg, ev2vg);
  // entity vertices to vertices, flipped so smallest global goes in front
  auto const ev2v_canon = align_ev2v(deg, ev2v, canon_codes);
  *p_ent_verts2verts = ev2v_canon;
  // entity vertices to vertex globals, with smallest global first
  auto const ev2vg_canon = align_ev2v(deg, ev2vg, canon_codes);
  // entity to adj vertex with smallest global
  auto const e2fv = get_component(ev2v_canon, deg, 0);
  auto const total_verts = find_total_globals(comm, verts2globs);
  // vertices to owners in linear partitioning of globals
  auto const v2ov = globals_to_linear_owners(comm, verts2globs, total_verts);
  // entities to owning vertices in linear partitioning of vertex globals
  auto const e2ov = unmap(e2fv, v2ov);
  auto const linsize = linear_partition_size(comm, total_verts);
  // dist from entities to linear partitioning of vertex globals
  auto const in_dist = Dist(comm, e2ov, linsize);
  // each MPI rank in the linear partitioning of vertex globals
  // is "serving" a contiguous subset of vertex globals, and hence
  // a subset of the vertices
  // Each of these vertices is "serving" all entities for which they
  // are the smallest-global adjacent vertex
  //
  // "served" entity vertices to vertex globals
  auto const sev2vg = in_dist.exch(ev2vg_canon, deg);
  auto out_dist = in_dist.invert();
  // "serving" vertices to served entities
  auto const sv2svse = out_dist.roots2items();
  // number of served entities
  auto const nse = out_dist.nitems();
  // items2dests() returns, for each served entity, the (rank, local index) pair
  // of where it originally came from
  //
  // owning served entity to original entity
  auto const se2orig = out_dist.items2dests();
  auto const se2orig_rank = se2orig.ranks;
  // the ordering of this array is critical for the resulting output to be
  // both deterministic and serial-parallel consistent.
  // it ensures that entities will be explored for matches in an order that
  // depends only on which rank they came from (which is unique), therefore
  // basing ownership on a minimum-rank rule.
  // note that otherwise the ordering of entities in these arrays is not even
  // deterministic, let alone partition-independent
  auto const svse2se = sort_locally_based_on_rank(sv2svse, se2orig_rank);
  // helper to make "serving vertices to served entities" look like an Adj
  auto const sv2se_codes = Read<I8>(nse, make_code(false, 0, 0));
  auto const sv2se = Adj(sv2svse, svse2se, sv2se_codes);
  // served entities to serving vertex
  auto const se2fsv = invert_fan(sv2svse);
  // find_matches_ex will match up all duplicate entities by assigning a
  // temporary owner to each set of duplicates. note that this owner is entirely
  // ordering-dependent, i.e. it will be just the first duplicate in the input
  // arrays, and the ordering of those arrays is not even deterministic!
  //
  // served entity to owning served entity
  Write<LO> se2ose;
  Write<I8> ignored_codes;
  constexpr bool allow_duplicates = true;
  find_matches_ex(deg, se2fsv, sev2vg, sev2vg, sv2se, &se2ose, &ignored_codes,
      allow_duplicates);
  // the problem with the temporary owner is that it is not deterministic, and
  // definitely not serial-parallel consistent. we need to choose an owner based
  // on a deterministic and serial-parallel consistent rule.
  // Two commonly used rules are:
  //   1) smallest rank
  //   2) rank who owns the fewest entities (tiebreaker is rank).
  //      The problem with this rule is that we are inside the code that
  //      determines ownership, so there is a bit of a chicken-and-egg obstacle
  //      to applying this rule.
  //
  //      TODO: fix non-determinism and partition-dependence of owner
  //
  // served entity to owning original entity
  auto const se2owner_orig = unmap(read(se2ose), se2orig);
  // remove roots2items so that the next exch() call will accept items, not
  // roots
  out_dist.set_roots2items(LOs());
  // send back, to each served entity, the (rank, local index) pair of its new
  // owner
  //
  // entities to owning entities
  auto const e2owner = out_dist.exch(se2owner_orig, 1);
  *p_ents2owners = e2owner;
}

void suggest_slices(
    GO total, I32 comm_size, I32 comm_rank, GO* p_begin, GO* p_end) {
  auto comm_size_gt = GO(comm_size);
  auto quot = total / comm_size_gt;
  auto rem = total % comm_size_gt;
  if (comm_rank < rem) {
    *p_begin = quot * comm_rank + comm_rank;
    *p_end = *p_begin + quot + 1;
  } else {
    *p_begin = quot * comm_rank + rem;
    *p_end = *p_begin + quot;
  }
}

void assemble_slices(CommPtr comm, Omega_h_Family family, Int dim,
    GO global_nelems, GO elem_offset, GOs conn_in, GO global_nverts,
    GO vert_offset, Reals vert_coords, Dist* p_slice_elems2elems,
    LOs* p_conn_out, Dist* p_slice_verts2verts) {
  auto verts_per_elem = element_degree(family, dim, 0);
  auto nslice_elems = divide_no_remainder(conn_in.size(), verts_per_elem);
  auto nslice_verts = divide_no_remainder(vert_coords.size(), dim);
  {  // require that slicing was done as suggested
    GO suggested_elem_offset, suggested_elem_end;
    suggest_slices(global_nelems, comm->size(), comm->rank(),
        &suggested_elem_offset, &suggested_elem_end);
    OMEGA_H_CHECK(suggested_elem_offset == elem_offset);
    OMEGA_H_CHECK(suggested_elem_end == elem_offset + nslice_elems);
    GO suggested_vert_offset, suggested_vert_end;
    suggest_slices(global_nverts, comm->size(), comm->rank(),
        &suggested_vert_offset, &suggested_vert_end);
    OMEGA_H_CHECK(suggested_vert_offset == vert_offset);
    OMEGA_H_CHECK(suggested_vert_end == vert_offset + nslice_verts);
  }
  // generate communication pattern from sliced vertices to their sliced
  // elements
  auto slice_elems2slice_verts = copies_to_linear_owners(comm, conn_in);
  auto slice_verts2slice_elems = slice_elems2slice_verts.invert();
  // compute sliced element coordinates
  auto slice_elem_vert_coords = slice_verts2slice_elems.exch(vert_coords, dim);
  auto slice_elem_coords_w = Write<Real>(nslice_elems * dim);
  auto f = OMEGA_H_LAMBDA(LO slice_elem) {
    for (Int i = 0; i < dim; ++i) {
      slice_elem_coords_w[slice_elem * dim + i] = 0.;
      for (Int j = 0; j < verts_per_elem; ++j) {
        slice_elem_coords_w[slice_elem * dim + i] +=
            slice_elem_vert_coords[(slice_elem * verts_per_elem + j) * dim + i];
      }
      slice_elem_coords_w[slice_elem * dim + i] /= verts_per_elem;
    }
  };
  parallel_for(nslice_elems, f, "sliced elem coords");
  auto slice_elem_coords = Reals(slice_elem_coords_w);
  // geometrically partition the elements (RIB)
  auto rib_masses = Reals(nslice_elems, 1.0);
  auto rib_tol = 2.0;
  auto elem_owners = identity_remotes(comm, nslice_elems);
  auto rib_coords = resize_vectors(slice_elem_coords, dim, 3);
  inertia::Rib rib_hints;
  inertia::recursively_bisect(
      comm, rib_tol, &rib_coords, &rib_masses, &elem_owners, &rib_hints);
  // communication pattern from sliced elements to partitioned elements
  auto elems2slice_elems = Dist{comm, elem_owners, nslice_elems};
  auto slice_elems2elems = elems2slice_elems.invert();
  // communication pattern from sliced vertices to partitioned elements
  auto slice_elem_vert_owners = slice_elems2slice_verts.items2dests();
  auto elem_vert_owners =
      slice_elems2elems.exch(slice_elem_vert_owners, verts_per_elem);
  auto elems2slice_verts = Dist{comm, elem_vert_owners, nslice_verts};
  // unique set of vertices needed for partitioned elements
  auto slice_vert_globals =
      GOs{nslice_verts, vert_offset, 1, "slice vert globals"};
  auto verts2slice_verts =
      get_new_copies2old_owners(elems2slice_verts, slice_vert_globals);
  // new (local) connectivity
  auto slice_verts2elems = elems2slice_verts.invert();
  auto new_conn = form_new_conn(verts2slice_verts, slice_verts2elems);
  *p_slice_elems2elems = slice_elems2elems;
  *p_conn_out = new_conn;
  *p_slice_verts2verts = verts2slice_verts.invert();
}

void build_quadratic_wireframe_3d(Mesh* mesh, Mesh* wireframe_mesh,
                               LO n_sample_pts) {
  auto dim = mesh->dim();
  auto nedge = mesh->nedges();
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto vert_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(0));
  auto vert_ctrlPts = mesh->get_ctrlPts(0);
  if (!vert_ctrlPts.exists()) {
    vert_ctrlPts_h = HostRead<Real>(mesh->coords());
  }
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);
  printf("max ctrl pts v %f e %f\n", get_max(vert_ctrlPts), get_max(mesh->get_ctrlPts(1)));

  Real xi_start = 0.0;
  Real xi_end = 1.0;
  Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
  auto u_h = HostRead<Real>(Read<Real>(n_sample_pts, xi_start, delta_xi,
                                      "samplePts"));

  HostWrite<Real> host_coords(n_sample_pts*nedge*dim);
  LO wireframe_nedge = (n_sample_pts - 1)*nedge;
  HostWrite<LO> host_ev2v(wireframe_nedge*2);
  std::vector<int> edge_vertices[1];
  edge_vertices[0].reserve(wireframe_nedge*2);
  LO count_wireframe_vtx = 0;

  auto gid = HostRead<LO>(mesh->get_array<LO>(1, "class_id"));
  auto gdim = HostRead<I8>(mesh->get_array<I8>(1, "class_dim"));

    for (LO i = 0; i < nedge; ++i) {
    auto v0 = ev2v_h[i*2];
    auto v1 = ev2v_h[i*2 + 1];
    Real cx0 = vert_ctrlPts_h[v0*dim + 0];
    Real cy0 = vert_ctrlPts_h[v0*dim + 1];
    Real cz0 = vert_ctrlPts_h[v0*dim + 2];
    Real cx2 = vert_ctrlPts_h[v1*dim + 0];
    Real cy2 = vert_ctrlPts_h[v1*dim + 1];
    Real cz2 = vert_ctrlPts_h[v1*dim + 2];

    //testing for fgrid, rmv later
    /*
    Real r_min=1e16;
    Real r = std::sqrt(cx0*cx0+cy0*cy0);
    if (r < r_min) r_min = r;
    r = std::sqrt(cx2*cx2+cy2*cy2);
    if (r < r_min) r_min = r;
    if ((r_min > 0.8) && (gdim[i] <= 2) && (gid[i] !=3) && (gid[i] !=13))  {
    */
    //

    Real cx1 = ctrlPts_h[i*dim + 0];
    Real cy1 = ctrlPts_h[i*dim + 1];
    Real cz1 = ctrlPts_h[i*dim + 2];

    for (LO i = 0; i < u_h.size(); ++i) {
      auto x_bezier = cx0*B0_quad(u_h[i]) + cx1*B1_quad(u_h[i]) +
                      cx2*B2_quad(u_h[i]);
      auto y_bezier = cy0*B0_quad(u_h[i]) + cy1*B1_quad(u_h[i]) +
                      cy2*B2_quad(u_h[i]);
      auto z_bezier = cz0*B0_quad(u_h[i]) + cz1*B1_quad(u_h[i]) +
                      cz2*B2_quad(u_h[i]);

      host_coords[count_wireframe_vtx*dim + 0] = x_bezier;
      host_coords[count_wireframe_vtx*dim + 1] = y_bezier;
      host_coords[count_wireframe_vtx*dim + 2] = z_bezier;

      edge_vertices[0].push_back(count_wireframe_vtx);
      if ((i > 0) && (i < (u_h.size() - 1))) {
        edge_vertices[0].push_back(count_wireframe_vtx);
      }

      ++count_wireframe_vtx;
    }
    //}
  }

  for (int i = 0; i < wireframe_nedge*2; ++i) {
    host_ev2v[i] = edge_vertices[0][static_cast<std::size_t>(i)];
  }

  wireframe_mesh->set_parting(OMEGA_H_ELEM_BASED);
  wireframe_mesh->set_dim(dim);
  wireframe_mesh->set_family(OMEGA_H_SIMPLEX);
  wireframe_mesh->set_verts(n_sample_pts*nedge);

  wireframe_mesh->add_coords(Reals(host_coords.write()));
  wireframe_mesh->set_ents(1, Adj(LOs(host_ev2v.write())));

  return;
}

void build_quadratic_wireframe_2d(Mesh* mesh, Mesh* wireframe_mesh,
                               LO n_sample_pts) {
  auto dim = mesh->dim();
  auto nedge = mesh->nedges();
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto vert_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(0));
  auto vert_ctrlPts = mesh->get_ctrlPts(0);
  if (!vert_ctrlPts.exists()) {
    vert_ctrlPts_h = HostRead<Real>(mesh->coords());
  }
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);

  Real xi_start = 0.0;
  Real xi_end = 1.0;
  Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
  auto u_h = HostRead<Real>(Read<Real>(n_sample_pts, xi_start, delta_xi,
                                      "samplePts"));

  HostWrite<Real> host_coords(n_sample_pts*nedge*dim);
  LO wireframe_nedge = (n_sample_pts - 1)*nedge;
  HostWrite<LO> host_ev2v(wireframe_nedge*2);
  std::vector<int> edge_vertices[1];
  edge_vertices[0].reserve(wireframe_nedge*2);
  LO count_wireframe_vtx = 0;

  for (LO i = 0; i < nedge; ++i) {
    auto v0 = ev2v_h[i*2];
    auto v1 = ev2v_h[i*2 + 1];

    /*
    if (dim == 2) {
      Vector<2> sample_pt;
      auto c0 = get_vector<2>(vert_ctrlPts_h, v0);
      auto c1 = get_vector<2>(ctrlPts_h, i);
      auto c2 = get_vector<2>(vert_ctrlPts_h, v1);
      for (LO i = 0; i < u_h.size(); ++i) {
        for (Int j = 0; j < dim; ++j) {
          c1[j] = (p1[j] - B0_quad(xi_1)*c0[j] - B2_quad(xi_1)*c2[j])/B1_quad(xi_1);
        }
      }
      set_vector(host_coords, count_wireframe_vtx, c1);
    }
    */

    Real cx0 = vert_ctrlPts_h[v0*dim + 0];
    Real cy0 = vert_ctrlPts_h[v0*dim + 1];
    //Real cz0 = vert_ctrlPts_h[v0*dim + 2];
    Real cx1 = ctrlPts_h[i*dim + 0];
    Real cy1 = ctrlPts_h[i*dim + 1];
    //Real cz1 = ctrlPts_h[i*dim + 2];
    Real cx2 = vert_ctrlPts_h[v1*dim + 0];
    Real cy2 = vert_ctrlPts_h[v1*dim + 1];
    //Real cz2 = vert_ctrlPts_h[v1*dim + 2];

    for (LO i = 0; i < u_h.size(); ++i) {
      auto x_bezier = cx0*B0_quad(u_h[i]) + cx1*B1_quad(u_h[i]) +
                      cx2*B2_quad(u_h[i]);
      auto y_bezier = cy0*B0_quad(u_h[i]) + cy1*B1_quad(u_h[i]) +
                      cy2*B2_quad(u_h[i]);
      //auto z_bezier = cz0*B0_quad(u_h[i]) + cz1*B1_quad(u_h[i]) +
       //               cz2*B2_quad(u_h[i]);

      host_coords[count_wireframe_vtx*dim + 0] = x_bezier;
      host_coords[count_wireframe_vtx*dim + 1] = y_bezier;
      //host_coords[count_wireframe_vtx*dim + 2] = z_bezier;

      edge_vertices[0].push_back(count_wireframe_vtx);
      if ((i > 0) && (i < (u_h.size() - 1))) {
        edge_vertices[0].push_back(count_wireframe_vtx);
      }

      ++count_wireframe_vtx;
    }
  }

  for (int i = 0; i < wireframe_nedge*2; ++i) {
    host_ev2v[i] = edge_vertices[0][static_cast<std::size_t>(i)];
  }

  wireframe_mesh->set_parting(OMEGA_H_ELEM_BASED);
  wireframe_mesh->set_dim(dim);
  wireframe_mesh->set_family(OMEGA_H_SIMPLEX);
  wireframe_mesh->set_verts(n_sample_pts*nedge);

  wireframe_mesh->add_coords(Reals(host_coords.write()));
  wireframe_mesh->set_ents(1, Adj(LOs(host_ev2v.write())));

  return;
}

void build_cubic_wireframe_2d(Mesh* mesh, Mesh* wireframe_mesh,
                              LO n_sample_pts) {
  auto dim = mesh->dim();
  auto nedge = mesh->nedges();
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto vert_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(0));
  auto vert_ctrlPts = mesh->get_ctrlPts(0);
  if (!vert_ctrlPts.exists()) {
    vert_ctrlPts_h = HostRead<Real>(mesh->coords());
  }
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);

  HostWrite<Real> host_coords(n_sample_pts*nedge*dim);
  LO wireframe_nedge = (n_sample_pts - 1)*nedge;
  HostWrite<LO> host_ev2v(wireframe_nedge*2);
  std::vector<int> edge_vertices[1];
  edge_vertices[0].reserve(wireframe_nedge*2);
  LO count_wireframe_vtx = 0;

  for (LO i = 0; i < nedge; ++i) {
    auto v0 = ev2v_h[i*2];
    auto v1 = ev2v_h[i*2 + 1];

    Real xi_start = 0.0;
    Real xi_end = 1.0;
    Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
    auto u_h = HostRead<Real>(Read<Real>(n_sample_pts, xi_start, delta_xi,
          "samplePts"));
 
    Real cx0 = vert_ctrlPts_h[v0*dim + 0];
    Real cy0 = vert_ctrlPts_h[v0*dim + 1];
    //Real cz0 = vert_ctrlPts_h[v0*dim + 2];
    Real cx1 = ctrlPts_h[i*pts_per_edge*dim + 0];
    Real cy1 = ctrlPts_h[i*pts_per_edge*dim + 1];
    //Real cz1 = ctrlPts_h[i*pts_per_edge*dim + 2];
    Real cx2 = ctrlPts_h[i*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
    Real cy2 = ctrlPts_h[i*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
    //Real cz2 = ctrlPts_h[i*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
    Real cx3 = vert_ctrlPts_h[v1*dim + 0];
    Real cy3 = vert_ctrlPts_h[v1*dim + 1];
    //Real cz3 = vert_ctrlPts_h[v0*dim + 2];

    for (LO i = 0; i < u_h.size(); ++i) {
      auto x_bezier = cx0*B0_cube(u_h[i]) + cx1*B1_cube(u_h[i]) +
                      cx2*B2_cube(u_h[i]) + cx3*B3_cube(u_h[i]);
      auto y_bezier = cy0*B0_cube(u_h[i]) + cy1*B1_cube(u_h[i]) +
                      cy2*B2_cube(u_h[i]) + cy3*B3_cube(u_h[i]);
      //auto z_bezier = cz0*B0_cube(u_h[i]) + cz1*B1_cube(u_h[i]) +
        //              cz2*B2_cube(u_h[i]) + cz3*B3_cube(u_h[i]);

      host_coords[count_wireframe_vtx*dim + 0] = x_bezier;
      host_coords[count_wireframe_vtx*dim + 1] = y_bezier;
      //host_coords[count_wireframe_vtx*dim + 2] = z_bezier;

      edge_vertices[0].push_back(count_wireframe_vtx);
      if ((i > 0) && (i < (u_h.size() - 1))) {
        edge_vertices[0].push_back(count_wireframe_vtx);
      }

      ++count_wireframe_vtx;
    }
  }

  for (int i = 0; i < wireframe_nedge*2; ++i) {
    host_ev2v[i] = edge_vertices[0][static_cast<std::size_t>(i)];
  }

  wireframe_mesh->set_parting(OMEGA_H_ELEM_BASED);
  wireframe_mesh->set_dim(dim);
  wireframe_mesh->set_family(OMEGA_H_SIMPLEX);
  wireframe_mesh->set_verts(n_sample_pts*nedge);

  wireframe_mesh->add_coords(Reals(host_coords.write()));
  wireframe_mesh->set_ents(1, Adj(LOs(host_ev2v.write())));

  return;
}

void build_cubic_wireframe_3d(Mesh* mesh, Mesh* wireframe_mesh,
                              LO n_sample_pts) {
  auto dim = mesh->dim();
  auto nedge = mesh->nedges();
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto vert_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(0));
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);

  HostWrite<Real> host_coords(n_sample_pts*nedge*dim);
  LO wireframe_nedge = (n_sample_pts - 1)*nedge;
  HostWrite<LO> host_ev2v(wireframe_nedge*2);
  std::vector<int> edge_vertices[1];
  edge_vertices[0].reserve(wireframe_nedge*2);
  LO count_wireframe_vtx = 0;

  for (LO i = 0; i < nedge; ++i) {
    auto v0 = ev2v_h[i*2];
    auto v1 = ev2v_h[i*2 + 1];

    Real xi_start = 0.0;
    Real xi_end = 1.0;
    Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
    auto u_h = HostRead<Real>(Read<Real>(n_sample_pts, xi_start, delta_xi,
          "samplePts"));
 
    Real cx0 = vert_ctrlPts_h[v0*dim + 0];
    Real cy0 = vert_ctrlPts_h[v0*dim + 1];
    Real cz0 = vert_ctrlPts_h[v0*dim + 2];
    Real cx1 = ctrlPts_h[i*pts_per_edge*dim + 0];
    Real cy1 = ctrlPts_h[i*pts_per_edge*dim + 1];
    Real cz1 = ctrlPts_h[i*pts_per_edge*dim + 2];
    Real cx2 = ctrlPts_h[i*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
    Real cy2 = ctrlPts_h[i*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
    Real cz2 = ctrlPts_h[i*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
    Real cx3 = vert_ctrlPts_h[v1*dim + 0];
    Real cy3 = vert_ctrlPts_h[v1*dim + 1];
    Real cz3 = vert_ctrlPts_h[v1*dim + 2];

    for (LO i = 0; i < u_h.size(); ++i) {
      auto x_bezier = cx0*B0_cube(u_h[i]) + cx1*B1_cube(u_h[i]) +
                      cx2*B2_cube(u_h[i]) + cx3*B3_cube(u_h[i]);
      auto y_bezier = cy0*B0_cube(u_h[i]) + cy1*B1_cube(u_h[i]) +
                      cy2*B2_cube(u_h[i]) + cy3*B3_cube(u_h[i]);
      auto z_bezier = cz0*B0_cube(u_h[i]) + cz1*B1_cube(u_h[i]) +
                      cz2*B2_cube(u_h[i]) + cz3*B3_cube(u_h[i]);

      host_coords[count_wireframe_vtx*dim + 0] = x_bezier;
      host_coords[count_wireframe_vtx*dim + 1] = y_bezier;
      host_coords[count_wireframe_vtx*dim + 2] = z_bezier;

      edge_vertices[0].push_back(count_wireframe_vtx);
      if ((i > 0) && (i < (u_h.size() - 1))) {
        edge_vertices[0].push_back(count_wireframe_vtx);
      }

      ++count_wireframe_vtx;
    }
  }

  for (int i = 0; i < wireframe_nedge*2; ++i) {
    host_ev2v[i] = edge_vertices[0][static_cast<std::size_t>(i)];
  }

  wireframe_mesh->set_parting(OMEGA_H_ELEM_BASED);
  wireframe_mesh->set_dim(dim);
  wireframe_mesh->set_family(OMEGA_H_SIMPLEX);
  wireframe_mesh->set_verts(n_sample_pts*nedge);

  wireframe_mesh->add_coords(Reals(host_coords.write()));
  wireframe_mesh->set_ents(1, Adj(LOs(host_ev2v.write())));

  return;
}

void build_quartic_wireframe(Mesh* mesh, Mesh* wireframe_mesh,
                             LO n_sample_pts) {
  auto nedge = mesh->nedges();
  auto coords_h = HostRead<Real>(mesh->coords());
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);

  I8 const dim = mesh->dim();
  Real xi_start = 0.0;
  Real xi_end = 1.0;
  Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
  auto u_h = HostRead<Real>(Read<Real>(n_sample_pts, xi_start, delta_xi,
                                      "samplePts"));

  HostWrite<Real> host_coords(n_sample_pts*nedge*dim);
  LO wireframe_nedge = (n_sample_pts - 1)*nedge;
  HostWrite<LO> host_ev2v(wireframe_nedge*2);
  std::vector<int> edge_vertices[1];
  edge_vertices[0].reserve(wireframe_nedge*2);
  LO count_wireframe_vtx = 0;

  for (LO i = 0; i < nedge; ++i) {
    auto v0 = ev2v_h[i*2];
    auto v1 = ev2v_h[i*2 + 1];

    if (dim == 3) {
      Real cx0 = coords_h[v0*dim + 0];
      Real cy0 = coords_h[v0*dim + 1];
      Real cz0 = coords_h[v0*dim + 2];
      Real cx1 = ctrlPts_h[i*pts_per_edge*dim + 0];
      Real cy1 = ctrlPts_h[i*pts_per_edge*dim + 1];
      Real cz1 = ctrlPts_h[i*pts_per_edge*dim + 2];
      Real cx2 = ctrlPts_h[i*pts_per_edge*dim + dim + 0];
      Real cy2 = ctrlPts_h[i*pts_per_edge*dim + dim + 1];
      Real cz2 = ctrlPts_h[i*pts_per_edge*dim + dim + 2];
      Real cx3 = ctrlPts_h[i*pts_per_edge*dim + dim + dim + 0];
      Real cy3 = ctrlPts_h[i*pts_per_edge*dim + dim + dim + 1];
      Real cz3 = ctrlPts_h[i*pts_per_edge*dim + dim + dim + 2];
      Real cx4 = coords_h[v1*dim + 0];
      Real cy4 = coords_h[v1*dim + 1];
      Real cz4 = coords_h[v1*dim + 2];

      for (LO i = 0; i < u_h.size(); ++i) {
        auto x_bezier = cx0*B0_quart(u_h[i]) + cx1*B1_quart(u_h[i]) +
          cx2*B2_quart(u_h[i]) + cx3*B3_quart(u_h[i]) +
          cx4*B4_quart(u_h[i]);
        auto y_bezier = cy0*B0_quart(u_h[i]) + cy1*B1_quart(u_h[i]) +
          cy2*B2_quart(u_h[i]) + cy3*B3_quart(u_h[i]) +
          cy4*B4_quart(u_h[i]);
        auto z_bezier = cz0*B0_quart(u_h[i]) + cz1*B1_quart(u_h[i]) +
          cz2*B2_quart(u_h[i]) + cz3*B3_quart(u_h[i]) +
          cz4*B4_quart(u_h[i]);

        host_coords[count_wireframe_vtx*dim + 0] = x_bezier;
        host_coords[count_wireframe_vtx*dim + 1] = y_bezier;
        host_coords[count_wireframe_vtx*dim + 2] = z_bezier;

        edge_vertices[0].push_back(count_wireframe_vtx);
        if ((i > 0) && (i < (u_h.size() - 1))) {
          edge_vertices[0].push_back(count_wireframe_vtx);
        }

        ++count_wireframe_vtx;
      }
    }
    else {
      OMEGA_H_CHECK (dim == 2);
      Real cx0 = coords_h[v0*dim + 0];
      Real cy0 = coords_h[v0*dim + 1];
      Real cx1 = ctrlPts_h[i*pts_per_edge*dim + 0];
      Real cy1 = ctrlPts_h[i*pts_per_edge*dim + 1];
      Real cx2 = ctrlPts_h[i*pts_per_edge*dim + dim + 0];
      Real cy2 = ctrlPts_h[i*pts_per_edge*dim + dim + 1];
      Real cx3 = ctrlPts_h[i*pts_per_edge*dim + dim + dim + 0];
      Real cy3 = ctrlPts_h[i*pts_per_edge*dim + dim + dim + 1];
      Real cx4 = coords_h[v1*dim + 0];
      Real cy4 = coords_h[v1*dim + 1];

      for (LO i = 0; i < u_h.size(); ++i) {
        auto x_bezier = cx0*B0_quart(u_h[i]) + cx1*B1_quart(u_h[i]) +
          cx2*B2_quart(u_h[i]) + cx3*B3_quart(u_h[i]) +
          cx4*B4_quart(u_h[i]);
        auto y_bezier = cy0*B0_quart(u_h[i]) + cy1*B1_quart(u_h[i]) +
          cy2*B2_quart(u_h[i]) + cy3*B3_quart(u_h[i]) +
          cy4*B4_quart(u_h[i]);

        host_coords[count_wireframe_vtx*dim + 0] = x_bezier;
        host_coords[count_wireframe_vtx*dim + 1] = y_bezier;

        edge_vertices[0].push_back(count_wireframe_vtx);
        if ((i > 0) && (i < (u_h.size() - 1))) {
          edge_vertices[0].push_back(count_wireframe_vtx);
        }

        ++count_wireframe_vtx;
      }
    }
  }

  for (int i = 0; i < wireframe_nedge*2; ++i) {
    host_ev2v[i] = edge_vertices[0][static_cast<std::size_t>(i)];
  }

  wireframe_mesh->set_parting(OMEGA_H_ELEM_BASED);
  wireframe_mesh->set_dim(dim);
  wireframe_mesh->set_family(OMEGA_H_SIMPLEX);
  wireframe_mesh->set_verts(n_sample_pts*nedge);

  wireframe_mesh->add_coords(Reals(host_coords.write()));
  wireframe_mesh->set_ents(1, Adj(LOs(host_ev2v.write())));

  return;
}

void build_quadratic_curveVtk_3d(Mesh* mesh, Mesh* curveVtk_mesh,
                              LO n_sample_pts) {
  auto nface = mesh->nfaces();
  //auto coords_h = HostRead<Real>(mesh->coords());
  auto coords_h = HostRead<Real>(mesh->get_ctrlPts(0));
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto fv2v_h = HostRead<LO>(mesh->ask_down(2, 0).ab2b);
  auto fe2e_h = HostRead<LO>(mesh->get_adj(2, 1).ab2b);

  I8 dim = 3;
  Real xi_start = 0.0;
  Real xi_end = 1.0;
  Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
  double xi[n_sample_pts][n_sample_pts][FACE] = {0.0};
  LO count_curveVtk_mesh_vtx_perTri = 0;

  for (LO i = 0; i < n_sample_pts; ++i) {
    for (LO j = 0; j < n_sample_pts - i; ++j) {
      xi[i][j][0] = i*delta_xi;
      xi[i][j][1] = j*delta_xi;
      ++count_curveVtk_mesh_vtx_perTri;
    }
  }

  HostWrite<Real> host_coords(count_curveVtk_mesh_vtx_perTri*nface*dim);
  LO faces_perTri = (n_sample_pts - 1)*(n_sample_pts - 1);
  LO curveVtk_mesh_nface = faces_perTri*nface;
  HostWrite<LO> host_fv2v(curveVtk_mesh_nface*3);
  std::vector<int> face_vertices[1];
  face_vertices[0].reserve(curveVtk_mesh_nface*3);
  
  auto gid = HostRead<LO>(mesh->get_array<LO>(2, "class_id"));
  auto gdim = HostRead<I8>(mesh->get_array<I8>(2, "class_dim"));

  LO count_curveVtk_mesh_vtx = 0;
  for (LO face = 0; face < nface; ++face) {
    auto v0 = fv2v_h[face*3];
    auto v1 = fv2v_h[face*3 + 1];
    auto v2 = fv2v_h[face*3 + 2];
    auto e0 = fe2e_h[face*3];
    auto e1 = fe2e_h[face*3 + 1];
    auto e2 = fe2e_h[face*3 + 2];

    Real cx00 = coords_h[v0*dim + 0];
    Real cy00 = coords_h[v0*dim + 1];
    Real cz00 = coords_h[v0*dim + 2];

    Real cx20 = coords_h[v1*dim + 0];
    Real cy20 = coords_h[v1*dim + 1];
    Real cz20 = coords_h[v1*dim + 2];

    Real cx02 = coords_h[v2*dim + 0];
    Real cy02 = coords_h[v2*dim + 1];
    Real cz02 = coords_h[v2*dim + 2];

    //testing for fgrid, rmv later
    /*
    Real r_min=1e16;
    Real r = std::sqrt(cx00*cx00+cy00*cy00);
    if (r < r_min) r_min = r;
    r = std::sqrt(cx20*cx20+cy20*cy20);
    if (r < r_min) r_min = r;
    r = std::sqrt(cx02*cx02+cy02*cy02);
    if (r < r_min) r_min = r;
    if ((r_min > 0.8) && (gdim[face] == 2) && (gid[face] !=3) && (gid[face] !=13))  {
    */
    //

    Real cx10 = ctrlPts_h[e0*dim + 0];
    Real cy10 = ctrlPts_h[e0*dim + 1];
    Real cz10 = ctrlPts_h[e0*dim + 2];

    Real cx11 = ctrlPts_h[e1*dim + 0];
    Real cy11 = ctrlPts_h[e1*dim + 1];
    Real cz11 = ctrlPts_h[e1*dim + 2];
 
    Real cx01 = ctrlPts_h[e2*dim + 0];
    Real cy01 = ctrlPts_h[e2*dim + 1];
    Real cz01 = ctrlPts_h[e2*dim + 2];

    for (LO i = 0; i < n_sample_pts; ++i) {
      for (LO j = 0; j < n_sample_pts - i; ++j) {
        auto x_bezier = cx00*B00_quad(xi[i][j][0], xi[i][j][1]) +
                        cx10*B10_quad(xi[i][j][0], xi[i][j][1]) +
                        cx20*B20_quad(xi[i][j][0], xi[i][j][1]) +
                        cx11*B11_quad(xi[i][j][0], xi[i][j][1]) +
                        cx02*B02_quad(xi[i][j][0], xi[i][j][1]) +
                        cx01*B01_quad(xi[i][j][0], xi[i][j][1]);
        auto y_bezier = cy00*B00_quad(xi[i][j][0], xi[i][j][1]) +
                        cy10*B10_quad(xi[i][j][0], xi[i][j][1]) +
                        cy20*B20_quad(xi[i][j][0], xi[i][j][1]) +
                        cy11*B11_quad(xi[i][j][0], xi[i][j][1]) +
                        cy02*B02_quad(xi[i][j][0], xi[i][j][1]) +
                        cy01*B01_quad(xi[i][j][0], xi[i][j][1]);
        auto z_bezier = cz00*B00_quad(xi[i][j][0], xi[i][j][1]) +
                        cz10*B10_quad(xi[i][j][0], xi[i][j][1]) +
                        cz20*B20_quad(xi[i][j][0], xi[i][j][1]) +
                        cz11*B11_quad(xi[i][j][0], xi[i][j][1]) +
                        cz02*B02_quad(xi[i][j][0], xi[i][j][1]) +
                        cz01*B01_quad(xi[i][j][0], xi[i][j][1]);

        host_coords[count_curveVtk_mesh_vtx*dim + 0] = x_bezier;
        host_coords[count_curveVtk_mesh_vtx*dim + 1] = y_bezier;
        host_coords[count_curveVtk_mesh_vtx*dim + 2] = z_bezier;

        if ((i < n_sample_pts - 1) && (j < n_sample_pts - i - 1)) {
          face_vertices[0].push_back(count_curveVtk_mesh_vtx);
          face_vertices[0].push_back(count_curveVtk_mesh_vtx + n_sample_pts-i);
          face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
          if (i > 0) {
            face_vertices[0].push_back(count_curveVtk_mesh_vtx);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1 - (n_sample_pts-(i-1)));
          }
        }

        ++count_curveVtk_mesh_vtx;
      }
    }
    //}
  }

  for (int i = 0; i < curveVtk_mesh_nface*3; ++i) {
    host_fv2v[i] = face_vertices[0][static_cast<std::size_t>(i)];
  }

  curveVtk_mesh->set_parting(OMEGA_H_ELEM_BASED);
  curveVtk_mesh->set_dim(dim);
  curveVtk_mesh->set_family(OMEGA_H_SIMPLEX);
  curveVtk_mesh->set_verts(count_curveVtk_mesh_vtx_perTri*nface);
  curveVtk_mesh->add_coords(Reals(host_coords.write()));
  curveVtk_mesh->set_down(2, 0, LOs(host_fv2v.write()));

  return;
}

void build_cubic_curveVtk_2d(Mesh* mesh, Mesh* curveVtk_mesh,
                          LO n_sample_pts) {
  auto nface = mesh->nfaces();
  auto coords_h = HostRead<Real>(mesh->coords());
  auto vert_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(0));
  auto vert_ctrlPts = mesh->get_ctrlPts(0);
  if (vert_ctrlPts.exists()) {
    coords_h = vert_ctrlPts_h;
  }
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto face_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(2));
  auto fv2v_h = HostRead<LO>(mesh->ask_down(2, 0).ab2b);
  auto fe2e_h = HostRead<LO>(mesh->get_adj(2, 1).ab2b);
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);

  I8 dim = mesh->dim();
  Real xi_start = 0.0;
  Real xi_end = 1.0;
  Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
  double xi[n_sample_pts][n_sample_pts][FACE] = {0.0};
  LO count_curveVtk_mesh_vtx_perTri = 0;

  for (LO i = 0; i < n_sample_pts; ++i) {
    for (LO j = 0; j < n_sample_pts - i; ++j) {
      xi[i][j][0] = i*delta_xi;
      xi[i][j][1] = j*delta_xi;
      ++count_curveVtk_mesh_vtx_perTri;
    }
  }

  HostWrite<Real> host_coords(count_curveVtk_mesh_vtx_perTri*nface*dim);
  LO faces_perTri = (n_sample_pts - 1)*(n_sample_pts - 1);
  LO curveVtk_mesh_nface = faces_perTri*nface;
  HostWrite<LO> host_fv2v(curveVtk_mesh_nface*3);
  std::vector<int> face_vertices[1];
  face_vertices[0].reserve(curveVtk_mesh_nface*3);
  
  LO count_curveVtk_mesh_vtx = 0;
  for (LO face = 0; face < nface; ++face) {
    auto v0 = fv2v_h[face*3];
    auto v1 = fv2v_h[face*3 + 1];
    auto v2 = fv2v_h[face*3 + 2];
    auto e0 = fe2e_h[face*3];
    auto e1 = fe2e_h[face*3 + 1];
    auto e2 = fe2e_h[face*3 + 2];

    I8 e0_flip = -1;
    I8 e1_flip = -1;
    I8 e2_flip = -1;
    auto e0v0 = ev2v_h[e0*2 + 0];
    auto e0v1 = ev2v_h[e0*2 + 1];
    auto e1v0 = ev2v_h[e1*2 + 0];
    auto e1v1 = ev2v_h[e1*2 + 1];
    auto e2v0 = ev2v_h[e2*2 + 0];
    auto e2v1 = ev2v_h[e2*2 + 1];
    if ((e0v0 == v1) && (e0v1 == v0)) {
      e0_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      e1_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      e2_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e2v0 == v2) && (e2v1 == v0));
    }

    Real cx00 = coords_h[v0*dim + 0];
    Real cy00 = coords_h[v0*dim + 1];
    //Real cz00 = coords_h[v0*dim + 2];
    Real cx30 = coords_h[v1*dim + 0];
    Real cy30 = coords_h[v1*dim + 1];
    //Real cz30 = coords_h[v1*dim + 2];
    Real cx03 = coords_h[v2*dim + 0];
    Real cy03 = coords_h[v2*dim + 1];
    //Real cz03 = coords_h[v2*dim + 2];

    Real cx10 = ctrlPts_h[e0*pts_per_edge*dim + 0];
    Real cy10 = ctrlPts_h[e0*pts_per_edge*dim + 1];
    //Real cz10 = ctrlPts_h[e0*pts_per_edge*dim + 2];
    Real cx20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
    Real cy20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
    //Real cz20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
    if (e0_flip > 0) {
      auto tempx = cx10;
      auto tempy = cy10;
      //auto tempz = cz10;
      cx10 = cx20;
      cy10 = cy20;
      //cz10 = cz20;
      cx20 = tempx;
      cy20 = tempy;
      //cz20 = tempz;
    }

    Real cx21 = ctrlPts_h[e1*pts_per_edge*dim + 0];
    Real cy21 = ctrlPts_h[e1*pts_per_edge*dim + 1];
    //Real cz21 = ctrlPts_h[e1*pts_per_edge*dim + 2];
    Real cx12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
    Real cy12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
    //Real cz12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
    if (e1_flip > 0) {
      auto tempx = cx21;
      auto tempy = cy21;
      //auto tempz = cz21;
      cx21 = cx12;
      cy21 = cy12;
      //cz21 = cz12;
      cx12 = tempx;
      cy12 = tempy;
      //cz12 = tempz;
    }

    Real cx02 = ctrlPts_h[e2*pts_per_edge*dim + 0];
    Real cy02 = ctrlPts_h[e2*pts_per_edge*dim + 1];
    //Real cz02 = ctrlPts_h[e2*pts_per_edge*dim + 2];
    Real cx01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
    Real cy01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
    //Real cz01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
    if (e2_flip > 0) {
      auto tempx = cx02;
      auto tempy = cy02;
      //auto tempz = cz02;
      cx02 = cx01;
      cy02 = cy01;
      //cz02 = cz01;
      cx01 = tempx;
      cy01 = tempy;
      //cz01 = tempz;
    }

    Real cx11 = face_ctrlPts_h[face*dim + 0];
    Real cy11 = face_ctrlPts_h[face*dim + 1];
    //Real cz11 = face_ctrlPts_h[face*dim + 2];

    for (LO i = 0; i < n_sample_pts; ++i) {
      for (LO j = 0; j < n_sample_pts - i; ++j) {
        auto x_bezier = cx00*B00_cube(xi[i][j][0], xi[i][j][1]) +
                        cx10*B10_cube(xi[i][j][0], xi[i][j][1]) +
                        cx20*B20_cube(xi[i][j][0], xi[i][j][1]) +
                        cx30*B30_cube(xi[i][j][0], xi[i][j][1]) +
                        cx21*B21_cube(xi[i][j][0], xi[i][j][1]) +
                        cx12*B12_cube(xi[i][j][0], xi[i][j][1]) +
                        cx03*B03_cube(xi[i][j][0], xi[i][j][1]) +
                        cx02*B02_cube(xi[i][j][0], xi[i][j][1]) +
                        cx01*B01_cube(xi[i][j][0], xi[i][j][1]) +
                        cx11*B11_cube(xi[i][j][0], xi[i][j][1]);
        auto y_bezier = cy00*B00_cube(xi[i][j][0], xi[i][j][1]) +
                        cy10*B10_cube(xi[i][j][0], xi[i][j][1]) +
                        cy20*B20_cube(xi[i][j][0], xi[i][j][1]) +
                        cy30*B30_cube(xi[i][j][0], xi[i][j][1]) +
                        cy21*B21_cube(xi[i][j][0], xi[i][j][1]) +
                        cy12*B12_cube(xi[i][j][0], xi[i][j][1]) +
                        cy03*B03_cube(xi[i][j][0], xi[i][j][1]) +
                        cy02*B02_cube(xi[i][j][0], xi[i][j][1]) +
                        cy01*B01_cube(xi[i][j][0], xi[i][j][1]) +
                        cy11*B11_cube(xi[i][j][0], xi[i][j][1]);

        host_coords[count_curveVtk_mesh_vtx*dim + 0] = x_bezier;
        host_coords[count_curveVtk_mesh_vtx*dim + 1] = y_bezier;

        if ((i < n_sample_pts - 1) && (j < n_sample_pts - i - 1)) {
          face_vertices[0].push_back(count_curveVtk_mesh_vtx);
          face_vertices[0].push_back(count_curveVtk_mesh_vtx + n_sample_pts-i);
          face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
          if (i > 0) {
            face_vertices[0].push_back(count_curveVtk_mesh_vtx);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1 - (n_sample_pts-(i-1)));
          }
        }

        ++count_curveVtk_mesh_vtx;
      }
    }
  }

  for (int i = 0; i < curveVtk_mesh_nface*3; ++i) {
    host_fv2v[i] = face_vertices[0][static_cast<std::size_t>(i)];
  }

  curveVtk_mesh->set_parting(OMEGA_H_ELEM_BASED);
  curveVtk_mesh->set_dim(dim);
  curveVtk_mesh->set_family(OMEGA_H_SIMPLEX);
  curveVtk_mesh->set_verts(count_curveVtk_mesh_vtx_perTri*nface);
  curveVtk_mesh->add_coords(Reals(host_coords.write()));
  curveVtk_mesh->set_down(2, 0, LOs(host_fv2v.write()));

  return;
}

void build_cubic_curveVtk_3d(Mesh* mesh, Mesh* curveVtk_mesh,
                          LO n_sample_pts) {
  auto nface = mesh->nfaces();
  auto coords_h = HostRead<Real>(mesh->coords());
  auto vert_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(0));
  auto vert_ctrlPts = mesh->get_ctrlPts(0);
  if (vert_ctrlPts.exists()) {
    coords_h = vert_ctrlPts_h;
  }
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto face_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(2));
  auto fv2v_h = HostRead<LO>(mesh->ask_down(2, 0).ab2b);
  auto fe2e_h = HostRead<LO>(mesh->get_adj(2, 1).ab2b);
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);

  I8 dim = mesh->dim();
  Real xi_start = 0.0;
  Real xi_end = 1.0;
  Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
  double xi[n_sample_pts][n_sample_pts][FACE] = {0.0};
  LO count_curveVtk_mesh_vtx_perTri = 0;

  for (LO i = 0; i < n_sample_pts; ++i) {
    for (LO j = 0; j < n_sample_pts - i; ++j) {
      xi[i][j][0] = i*delta_xi;
      xi[i][j][1] = j*delta_xi;
      ++count_curveVtk_mesh_vtx_perTri;
    }
  }

  HostWrite<Real> host_coords(count_curveVtk_mesh_vtx_perTri*nface*dim);
  LO faces_perTri = (n_sample_pts - 1)*(n_sample_pts - 1);
  LO curveVtk_mesh_nface = faces_perTri*nface;
  HostWrite<LO> host_fv2v(curveVtk_mesh_nface*3);
  std::vector<int> face_vertices[1];
  face_vertices[0].reserve(curveVtk_mesh_nface*3);
  
  LO count_curveVtk_mesh_vtx = 0;
  for (LO face = 0; face < nface; ++face) {
    auto v0 = fv2v_h[face*3];
    auto v1 = fv2v_h[face*3 + 1];
    auto v2 = fv2v_h[face*3 + 2];
    auto e0 = fe2e_h[face*3];
    auto e1 = fe2e_h[face*3 + 1];
    auto e2 = fe2e_h[face*3 + 2];

    I8 e0_flip = -1;
    I8 e1_flip = -1;
    I8 e2_flip = -1;
    auto e0v0 = ev2v_h[e0*2 + 0];
    auto e0v1 = ev2v_h[e0*2 + 1];
    auto e1v0 = ev2v_h[e1*2 + 0];
    auto e1v1 = ev2v_h[e1*2 + 1];
    auto e2v0 = ev2v_h[e2*2 + 0];
    auto e2v1 = ev2v_h[e2*2 + 1];
    if ((e0v0 == v1) && (e0v1 == v0)) {
      e0_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      e1_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      e2_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e2v0 == v2) && (e2v1 == v0));
    }

    Real cx00 = coords_h[v0*dim + 0];
    Real cy00 = coords_h[v0*dim + 1];
    Real cz00 = coords_h[v0*dim + 2];
    Real cx30 = coords_h[v1*dim + 0];
    Real cy30 = coords_h[v1*dim + 1];
    Real cz30 = coords_h[v1*dim + 2];
    Real cx03 = coords_h[v2*dim + 0];
    Real cy03 = coords_h[v2*dim + 1];
    Real cz03 = coords_h[v2*dim + 2];

    Real cx10 = ctrlPts_h[e0*pts_per_edge*dim + 0];
    Real cy10 = ctrlPts_h[e0*pts_per_edge*dim + 1];
    Real cz10 = ctrlPts_h[e0*pts_per_edge*dim + 2];
    Real cx20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
    Real cy20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
    Real cz20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
    if (e0_flip > 0) {
      auto tempx = cx10;
      auto tempy = cy10;
      auto tempz = cz10;
      cx10 = cx20;
      cy10 = cy20;
      cz10 = cz20;
      cx20 = tempx;
      cy20 = tempy;
      cz20 = tempz;
    }

    Real cx21 = ctrlPts_h[e1*pts_per_edge*dim + 0];
    Real cy21 = ctrlPts_h[e1*pts_per_edge*dim + 1];
    Real cz21 = ctrlPts_h[e1*pts_per_edge*dim + 2];
    Real cx12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
    Real cy12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
    Real cz12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
    if (e1_flip > 0) {
      auto tempx = cx21;
      auto tempy = cy21;
      auto tempz = cz21;
      cx21 = cx12;
      cy21 = cy12;
      cz21 = cz12;
      cx12 = tempx;
      cy12 = tempy;
      cz12 = tempz;
    }

    Real cx02 = ctrlPts_h[e2*pts_per_edge*dim + 0];
    Real cy02 = ctrlPts_h[e2*pts_per_edge*dim + 1];
    Real cz02 = ctrlPts_h[e2*pts_per_edge*dim + 2];
    Real cx01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
    Real cy01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
    Real cz01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
    if (e2_flip > 0) {
      auto tempx = cx02;
      auto tempy = cy02;
      auto tempz = cz02;
      cx02 = cx01;
      cy02 = cy01;
      cz02 = cz01;
      cx01 = tempx;
      cy01 = tempy;
      cz01 = tempz;
    }

    Real cx11 = face_ctrlPts_h[face*dim + 0];
    Real cy11 = face_ctrlPts_h[face*dim + 1];
    Real cz11 = face_ctrlPts_h[face*dim + 2];

    for (LO i = 0; i < n_sample_pts; ++i) {
      for (LO j = 0; j < n_sample_pts - i; ++j) {
        auto x_bezier = cx00*B00_cube(xi[i][j][0], xi[i][j][1]) +
                        cx10*B10_cube(xi[i][j][0], xi[i][j][1]) +
                        cx20*B20_cube(xi[i][j][0], xi[i][j][1]) +
                        cx30*B30_cube(xi[i][j][0], xi[i][j][1]) +
                        cx21*B21_cube(xi[i][j][0], xi[i][j][1]) +
                        cx12*B12_cube(xi[i][j][0], xi[i][j][1]) +
                        cx03*B03_cube(xi[i][j][0], xi[i][j][1]) +
                        cx02*B02_cube(xi[i][j][0], xi[i][j][1]) +
                        cx01*B01_cube(xi[i][j][0], xi[i][j][1]) +
                        cx11*B11_cube(xi[i][j][0], xi[i][j][1]);
        auto y_bezier = cy00*B00_cube(xi[i][j][0], xi[i][j][1]) +
                        cy10*B10_cube(xi[i][j][0], xi[i][j][1]) +
                        cy20*B20_cube(xi[i][j][0], xi[i][j][1]) +
                        cy30*B30_cube(xi[i][j][0], xi[i][j][1]) +
                        cy21*B21_cube(xi[i][j][0], xi[i][j][1]) +
                        cy12*B12_cube(xi[i][j][0], xi[i][j][1]) +
                        cy03*B03_cube(xi[i][j][0], xi[i][j][1]) +
                        cy02*B02_cube(xi[i][j][0], xi[i][j][1]) +
                        cy01*B01_cube(xi[i][j][0], xi[i][j][1]) +
                        cy11*B11_cube(xi[i][j][0], xi[i][j][1]);
        auto z_bezier = cz00*B00_cube(xi[i][j][0], xi[i][j][1]) +
                        cz10*B10_cube(xi[i][j][0], xi[i][j][1]) +
                        cz20*B20_cube(xi[i][j][0], xi[i][j][1]) +
                        cz30*B30_cube(xi[i][j][0], xi[i][j][1]) +
                        cz21*B21_cube(xi[i][j][0], xi[i][j][1]) +
                        cz12*B12_cube(xi[i][j][0], xi[i][j][1]) +
                        cz03*B03_cube(xi[i][j][0], xi[i][j][1]) +
                        cz02*B02_cube(xi[i][j][0], xi[i][j][1]) +
                        cz01*B01_cube(xi[i][j][0], xi[i][j][1]) +
                        cz11*B11_cube(xi[i][j][0], xi[i][j][1]);

        host_coords[count_curveVtk_mesh_vtx*dim + 0] = x_bezier;
        host_coords[count_curveVtk_mesh_vtx*dim + 1] = y_bezier;
        host_coords[count_curveVtk_mesh_vtx*dim + 2] = z_bezier;

        if ((i < n_sample_pts - 1) && (j < n_sample_pts - i - 1)) {
          face_vertices[0].push_back(count_curveVtk_mesh_vtx);
          face_vertices[0].push_back(count_curveVtk_mesh_vtx + n_sample_pts-i);
          face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
          if (i > 0) {
            face_vertices[0].push_back(count_curveVtk_mesh_vtx);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1 - (n_sample_pts-(i-1)));
          }
        }

        ++count_curveVtk_mesh_vtx;
      }
    }
  }

  for (int i = 0; i < curveVtk_mesh_nface*3; ++i) {
    host_fv2v[i] = face_vertices[0][static_cast<std::size_t>(i)];
  }

  curveVtk_mesh->set_parting(OMEGA_H_ELEM_BASED);
  curveVtk_mesh->set_dim(dim);
  curveVtk_mesh->set_family(OMEGA_H_SIMPLEX);
  curveVtk_mesh->set_verts(count_curveVtk_mesh_vtx_perTri*nface);
  curveVtk_mesh->add_coords(Reals(host_coords.write()));
  curveVtk_mesh->set_down(2, 0, LOs(host_fv2v.write()));

  return;
}

void build_cubic_cavities_3d(Mesh* mesh, Mesh* curveVtk_mesh,
                          LO n_sample_pts) {
  auto coords_h = HostRead<Real>(mesh->coords());
  auto vert_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(0));
  auto face_crvVis_h = HostRead<LO>(mesh->get_array<LO>(2, "face_crvVis"));
  LO nface = 0;
  for (LO i = 0; i < face_crvVis_h.size(); ++i) {
    if (face_crvVis_h[i] == 1) ++nface;
  }
  //printf("building %d nfaces of cavities\n", nface);

  auto vert_ctrlPts = mesh->get_ctrlPts(0);
  if (vert_ctrlPts.exists()) {
    coords_h = vert_ctrlPts_h;
  }
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto face_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(2));
  auto fv2v_h = HostRead<LO>(mesh->ask_down(2, 0).ab2b);
  auto fe2e_h = HostRead<LO>(mesh->get_adj(2, 1).ab2b);
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);

  I8 dim = mesh->dim();
  Real xi_start = 0.0;
  Real xi_end = 1.0;
  Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
  double xi[n_sample_pts][n_sample_pts][FACE] = {0.0};
  LO count_curveVtk_mesh_vtx_perTri = 0;

  for (LO i = 0; i < n_sample_pts; ++i) {
    for (LO j = 0; j < n_sample_pts - i; ++j) {
      xi[i][j][0] = i*delta_xi;
      xi[i][j][1] = j*delta_xi;
      ++count_curveVtk_mesh_vtx_perTri;
    }
  }

  HostWrite<Real> host_coords(count_curveVtk_mesh_vtx_perTri*nface*dim);
  LO faces_perTri = (n_sample_pts - 1)*(n_sample_pts - 1);
  LO curveVtk_mesh_nface = faces_perTri*nface;
  HostWrite<LO> host_fv2v(curveVtk_mesh_nface*3);
  std::vector<int> face_vertices[1];
  face_vertices[0].reserve(curveVtk_mesh_nface*3);
  
  LO count_curveVtk_mesh_vtx = 0;
  for (LO face = 0; face < mesh->nfaces(); ++face) {
    if (face_crvVis_h[face] == 1) {
      auto v0 = fv2v_h[face*3];
      auto v1 = fv2v_h[face*3 + 1];
      auto v2 = fv2v_h[face*3 + 2];
      auto e0 = fe2e_h[face*3];
      auto e1 = fe2e_h[face*3 + 1];
      auto e2 = fe2e_h[face*3 + 2];

      I8 e0_flip = -1;
      I8 e1_flip = -1;
      I8 e2_flip = -1;
      auto e0v0 = ev2v_h[e0*2 + 0];
      auto e0v1 = ev2v_h[e0*2 + 1];
      auto e1v0 = ev2v_h[e1*2 + 0];
      auto e1v1 = ev2v_h[e1*2 + 1];
      auto e2v0 = ev2v_h[e2*2 + 0];
      auto e2v1 = ev2v_h[e2*2 + 1];
      if ((e0v0 == v1) && (e0v1 == v0)) {
        e0_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
      }
      if ((e1v0 == v2) && (e1v1 == v1)) {
        e1_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
      }
      if ((e2v0 == v0) && (e2v1 == v2)) {
        e2_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e2v0 == v2) && (e2v1 == v0));
      }

      Real cx00 = coords_h[v0*dim + 0];
      Real cy00 = coords_h[v0*dim + 1];
      Real cz00 = coords_h[v0*dim + 2];
      Real cx30 = coords_h[v1*dim + 0];
      Real cy30 = coords_h[v1*dim + 1];
      Real cz30 = coords_h[v1*dim + 2];
      Real cx03 = coords_h[v2*dim + 0];
      Real cy03 = coords_h[v2*dim + 1];
      Real cz03 = coords_h[v2*dim + 2];

      Real cx10 = ctrlPts_h[e0*pts_per_edge*dim + 0];
      Real cy10 = ctrlPts_h[e0*pts_per_edge*dim + 1];
      Real cz10 = ctrlPts_h[e0*pts_per_edge*dim + 2];
      Real cx20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
      Real cy20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      Real cz20 = ctrlPts_h[e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        auto tempz = cz10;
        cx10 = cx20;
        cy10 = cy20;
        cz10 = cz20;
        cx20 = tempx;
        cy20 = tempy;
        cz20 = tempz;
      }

      Real cx21 = ctrlPts_h[e1*pts_per_edge*dim + 0];
      Real cy21 = ctrlPts_h[e1*pts_per_edge*dim + 1];
      Real cz21 = ctrlPts_h[e1*pts_per_edge*dim + 2];
      Real cx12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
      Real cy12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      Real cz12 = ctrlPts_h[e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e1_flip > 0) {
        auto tempx = cx21;
        auto tempy = cy21;
        auto tempz = cz21;
        cx21 = cx12;
        cy21 = cy12;
        cz21 = cz12;
        cx12 = tempx;
        cy12 = tempy;
        cz12 = tempz;
      }

      Real cx02 = ctrlPts_h[e2*pts_per_edge*dim + 0];
      Real cy02 = ctrlPts_h[e2*pts_per_edge*dim + 1];
      Real cz02 = ctrlPts_h[e2*pts_per_edge*dim + 2];
      Real cx01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
      Real cy01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      Real cz01 = ctrlPts_h[e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e2_flip > 0) {
        auto tempx = cx02;
        auto tempy = cy02;
        auto tempz = cz02;
        cx02 = cx01;
        cy02 = cy01;
        cz02 = cz01;
        cx01 = tempx;
        cy01 = tempy;
        cz01 = tempz;
      }

      Real cx11 = face_ctrlPts_h[face*dim + 0];
      Real cy11 = face_ctrlPts_h[face*dim + 1];
      Real cz11 = face_ctrlPts_h[face*dim + 2];

      for (LO i = 0; i < n_sample_pts; ++i) {
        for (LO j = 0; j < n_sample_pts - i; ++j) {
          auto x_bezier = cx00*B00_cube(xi[i][j][0], xi[i][j][1]) +
            cx10*B10_cube(xi[i][j][0], xi[i][j][1]) +
            cx20*B20_cube(xi[i][j][0], xi[i][j][1]) +
            cx30*B30_cube(xi[i][j][0], xi[i][j][1]) +
            cx21*B21_cube(xi[i][j][0], xi[i][j][1]) +
            cx12*B12_cube(xi[i][j][0], xi[i][j][1]) +
            cx03*B03_cube(xi[i][j][0], xi[i][j][1]) +
            cx02*B02_cube(xi[i][j][0], xi[i][j][1]) +
            cx01*B01_cube(xi[i][j][0], xi[i][j][1]) +
            cx11*B11_cube(xi[i][j][0], xi[i][j][1]);
          auto y_bezier = cy00*B00_cube(xi[i][j][0], xi[i][j][1]) +
            cy10*B10_cube(xi[i][j][0], xi[i][j][1]) +
            cy20*B20_cube(xi[i][j][0], xi[i][j][1]) +
            cy30*B30_cube(xi[i][j][0], xi[i][j][1]) +
            cy21*B21_cube(xi[i][j][0], xi[i][j][1]) +
            cy12*B12_cube(xi[i][j][0], xi[i][j][1]) +
            cy03*B03_cube(xi[i][j][0], xi[i][j][1]) +
            cy02*B02_cube(xi[i][j][0], xi[i][j][1]) +
            cy01*B01_cube(xi[i][j][0], xi[i][j][1]) +
            cy11*B11_cube(xi[i][j][0], xi[i][j][1]);
          auto z_bezier = cz00*B00_cube(xi[i][j][0], xi[i][j][1]) +
            cz10*B10_cube(xi[i][j][0], xi[i][j][1]) +
            cz20*B20_cube(xi[i][j][0], xi[i][j][1]) +
            cz30*B30_cube(xi[i][j][0], xi[i][j][1]) +
            cz21*B21_cube(xi[i][j][0], xi[i][j][1]) +
            cz12*B12_cube(xi[i][j][0], xi[i][j][1]) +
            cz03*B03_cube(xi[i][j][0], xi[i][j][1]) +
            cz02*B02_cube(xi[i][j][0], xi[i][j][1]) +
            cz01*B01_cube(xi[i][j][0], xi[i][j][1]) +
            cz11*B11_cube(xi[i][j][0], xi[i][j][1]);

          host_coords[count_curveVtk_mesh_vtx*dim + 0] = x_bezier;
          host_coords[count_curveVtk_mesh_vtx*dim + 1] = y_bezier;
          host_coords[count_curveVtk_mesh_vtx*dim + 2] = z_bezier;

          if ((i < n_sample_pts - 1) && (j < n_sample_pts - i - 1)) {
            face_vertices[0].push_back(count_curveVtk_mesh_vtx);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + n_sample_pts-i);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
            if (i > 0) {
              face_vertices[0].push_back(count_curveVtk_mesh_vtx);
              face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
              face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1 - (n_sample_pts-(i-1)));
            }
          }

          ++count_curveVtk_mesh_vtx;
        }
      }
    }
  }

  for (int i = 0; i < curveVtk_mesh_nface*3; ++i) {
    host_fv2v[i] = face_vertices[0][static_cast<std::size_t>(i)];
  }

  curveVtk_mesh->set_parting(OMEGA_H_ELEM_BASED);
  curveVtk_mesh->set_dim(dim);
  curveVtk_mesh->set_family(OMEGA_H_SIMPLEX);
  curveVtk_mesh->set_verts(count_curveVtk_mesh_vtx_perTri*nface);
  curveVtk_mesh->add_coords(Reals(host_coords.write()));
  curveVtk_mesh->set_down(2, 0, LOs(host_fv2v.write()));

  return;
}

/*
 * we know the Id of tets in a cav, iterate over them and count unique faces,
 * edges, verts doing this requires device to host data transfer
*/
void build_given_tets(Mesh* mesh, Mesh *full_mesh, Read<I8> build_tet, 
    Read<I8> build_face, Read<I8> build_edge, Read<I8> build_vert,
    LOs const keys2verts) {
  //TODO calculate classinfo correctly using make_cavity_class

  auto full_coords = full_mesh->coords();
  auto full_rv2v = full_mesh->ask_down(3, 0).ab2b;
  auto full_fv2v = full_mesh->ask_down(2, 0).ab2b;
  auto full_re2e = full_mesh->ask_down(3, 1).ab2b;
  auto full_ev2v = full_mesh->get_adj(1, 0).ab2b;
  auto full_rf2f = full_mesh->get_adj(3, 2).ab2b;
  auto full_classids_v = full_mesh->get_array<LO>(0, "class_id");
  auto full_classdim_v = full_mesh->get_array<I8>(0, "cav_classdim");
  auto full_classids_e = full_mesh->get_array<LO>(1, "class_id");
  auto full_classdim_e = full_mesh->get_array<I8>(1, "cav_classdim");
  auto full_classids_f = full_mesh->get_array<LO>(2, "class_id");
  auto full_classdim_f = full_mesh->get_array<I8>(2, "cav_classdim");
  auto full_classids_r = full_mesh->get_array<LO>(3, "class_id");
  auto full_vtxCtrlPts = full_mesh->get_ctrlPts(0);
  auto full_edgeCtrlPts = full_mesh->get_ctrlPts(EDGE);
  auto full_faceCtrlPts = full_mesh->get_ctrlPts(FACE);

  // we know which tets need to be counted
  // then count id of new verts, edges, faces, tets on the cpu
  // convert those ids to gpu and then transfer the coods, class, and topo
  // info on gpu
  auto build_tet_h = HostRead<I8>(build_tet);
  auto build_face_h = HostRead<I8>(build_face);
  auto build_edge_h = HostRead<I8>(build_edge);
  auto build_vert_h = HostRead<I8>(build_vert);
  HostWrite<I8> full2cav_tet_h(full_mesh->nregions());
  HostWrite<I8> full2cav_face_h(full_mesh->nfaces());
  HostWrite<I8> full2cav_edge_h(full_mesh->nedges());
  HostWrite<I8> full2cav_vert_h(full_mesh->nverts());

  make_cavity_class(full_mesh, keys2verts);//calcs class dim for all cavs

  LO numRegions = 0;
  for (LO t=0; t<build_tet.size(); ++t) {
    if (build_tet_h[t] > 0) {
      printf("cav tet id %d fullid %d\n", numRegions, t);
      full2cav_tet_h[t] = numRegions;
      ++numRegions;
    }
    else {
      full2cav_tet_h[t] = -1;
    }
  }
  LO numFaces = 0;
  for (LO f = 0; f< build_face.size(); ++f) {
    if (build_face_h[f] > 0) {
      printf("cav faceid %d fullid %d\n", numFaces, f);
      full2cav_face_h[f] = numFaces;
      ++numFaces;
    }
    else {
      full2cav_face_h[f] = -1;
    }
  }
  LO numEdges = 0;
  for (LO e = 0; e< build_edge.size(); ++e) {
    if (build_edge_h[e] > 0) {
      printf("cav edgeid %d fullid %d\n", numEdges, e);
      full2cav_edge_h[e] = numEdges;
      ++numEdges;
    }
    else {
      full2cav_edge_h[e] = -1;
    }
  }
  LO numVtx = 0;
  for (LO v = 0; v< build_vert.size(); ++v) {
    if (build_vert_h[v] > 0) {
      printf("cav vertid %d fullid %d\n", numVtx, v);
      full2cav_vert_h[v] = numVtx;
      ++numVtx;
    }
    else {
      full2cav_vert_h[v] = -1;
    }
  }
  auto full2cav_tet = full2cav_tet_h.write();
  auto full2cav_face = full2cav_face_h.write();
  auto full2cav_edge = full2cav_edge_h.write();
  auto full2cav_vert = full2cav_vert_h.write();

  mesh->set_parting(OMEGA_H_ELEM_BASED);
  Int max_dim=3;
  mesh->set_dim(max_dim);
  mesh->set_family(OMEGA_H_SIMPLEX);
  mesh->set_curved(1);
  mesh->set_max_order(3);
  fprintf(stderr, "tri=%d, tet=%d\n", numFaces, numRegions);

  Write<LO> rgn_vertices(4*numRegions);
  Write<LO> face_vertices(3*numFaces);
  Write<LO> edge_vertices(2*numEdges);
  Write<LO> class_ids_vtx(numVtx);
  Write<I8> class_dim_vtx(numVtx);
  Write<LO> class_ids_edge(numEdges);
  Write<I8> class_dim_edge(numEdges);
  Write<LO> class_ids_face(numFaces);
  Write<I8> class_dim_face(numFaces);
  Write<LO> class_ids_rgn(numRegions);
  Write<I8> class_dim_rgn(numRegions);

  Write<Real> coords(numVtx*max_dim);
  Write<Real> vtxPt_coords(numVtx*max_dim);
  auto transfer_v = OMEGA_H_LAMBDA (LO const v) {
    LO const cav_v = full2cav_vert[v];
    if (cav_v >= 0) {
      printf("v=%d, cavv=%d, numV %d\n", v, cav_v, numVtx);
      for (I8 j=0; j<max_dim; j++) {
        coords[cav_v*max_dim + j] = full_coords[v*max_dim + j];
        vtxPt_coords[cav_v*max_dim + j] = full_vtxCtrlPts[v*max_dim + j];
      }
      class_ids_vtx[cav_v] = full_classids_v[v];
      class_dim_vtx[cav_v] = full_classdim_v[v];
    }
  };
  parallel_for(full_mesh->nverts(), std::move(transfer_v));
  mesh->set_verts(numVtx);
  mesh->add_coords(Reals(coords));
  mesh->add_tag<ClassId>(0, "class_id", 1, Read<ClassId>(class_ids_vtx));
  mesh->add_tag<I8>(0, "class_dim", 1, Read<I8>(class_dim_vtx));

  fprintf(stderr, "ok0\n");
  I8 edge_numPts = 2;
  Write<Real> edgePt_coords(numEdges*edge_numPts*max_dim);
  auto transfer_e = OMEGA_H_LAMBDA (LO const e) {
    LO const cav_e = full2cav_edge[e];
    if (cav_e >= 0) {
      for(I8 j=0; j<max_dim; ++j) {
        edgePt_coords[cav_e*2*max_dim+j] = full_edgeCtrlPts[e*2*max_dim+j];
        edgePt_coords[cav_e*2*max_dim+max_dim+j] = 
          full_edgeCtrlPts[e*2*max_dim+max_dim+j];
      }
      for (I8 j=0; j<2; ++j) {
        printf("j %d\n", j);
        LO const vtx = full_ev2v[e*2+j];
        LO const cav_v = full2cav_vert[vtx];
        edge_vertices[cav_e*2+j] = cav_v;
      }
      class_ids_edge[cav_e] = full_classids_e[e];
      class_dim_edge[cav_e] = full_classdim_e[e];
    }
  };
  parallel_for(full_mesh->nedges(), std::move(transfer_e));
  mesh->set_ents(1, Adj(LOs(edge_vertices)));
  mesh->add_tag<ClassId>(1, "class_id", 1,
      Read<ClassId>(class_ids_edge));
  mesh->add_tag<I8>(1, "class_dim", 1,
      Read<I8>(class_dim_edge));
  fprintf(stderr, "ok1\n");

  Write<Real> facePt_coords(numFaces*max_dim);
  auto transfer_f = OMEGA_H_LAMBDA (LO const f) {
    LO const cav_f = full2cav_face[f];
    if (cav_f >= 0) {
      for(int j=0; j<max_dim; j++) {
        facePt_coords[cav_f*max_dim + j] = full_faceCtrlPts[f*max_dim + j];
      }
      for (int j=0; j<3; ++j) {
        LO const vtx = full_fv2v[f*3+j];
        LO const cav_v = full2cav_vert[vtx];
        face_vertices[cav_f*3+j] = cav_v;
      }
      class_ids_face[cav_f] = full_classids_f[f];
      class_dim_face[cav_f] = full_classdim_f[f];
    }
  };
  parallel_for(full_mesh->nfaces(), std::move(transfer_f));
  fprintf(stderr, "ok2\n");
  Adj edge2vert;
  Adj vert2edge;
  edge2vert = mesh->get_adj(1, 0);
  fprintf(stderr, "ok2.1\n");
  vert2edge = mesh->ask_up(0, 1);
  fprintf(stderr, "ok2.2\n");
  auto tri2verts = Read<LO>(face_vertices);
  Adj down;
  {
    down = reflect_down(tri2verts, edge2vert.ab2b, vert2edge,
	OMEGA_H_SIMPLEX, 2, 1);
    mesh->set_ents(2, down);
  }
  mesh->add_tag<ClassId>(2, "class_id", 1,
      Read<ClassId>(class_ids_face));
  mesh->add_tag<I8>(2, "class_dim", 1,
      Read<I8>(class_dim_face));
  fprintf(stderr, "ok3\n");

  auto transfer_r = OMEGA_H_LAMBDA (LO const r) {
    LO const cav_r = full2cav_tet[r];
    if (cav_r >= 0) {
      for (int j=0; j<4; ++j) {
        LO const vtx = full_rv2v[r*4+j];
        LO const cav_v = full2cav_vert[vtx];
        rgn_vertices[cav_r*4+j] = cav_v;
      }
      class_ids_rgn[cav_r] = full_classids_r[r];
      class_dim_rgn[cav_r] = 3;
    }
  };
  parallel_for(full_mesh->nregions(), std::move(transfer_r));

  Adj tri2vert;
  Adj vert2tri;
  tri2vert = mesh->ask_down(2, 0);
  vert2tri = mesh->ask_up(0, 2);
  auto tet2verts = Read<LO>(rgn_vertices);
  fprintf(stderr, "ok4\n");
  {
    down = reflect_down(tet2verts, tri2vert.ab2b, vert2tri,
	OMEGA_H_SIMPLEX, 3, 2);
    mesh->set_ents(3, down);
  }
  mesh->add_tag<ClassId>(3, "class_id", 1,
      Read<ClassId>(class_ids_rgn));
  mesh->add_tag<I8>(3, "class_dim", 1, Read<I8>(class_dim_rgn));

  fprintf(stderr, "ok5\n");
  for (LO i=0; i<=mesh->dim(); ++i) {
    if (!mesh->has_tag(i, "global")) {
      mesh->add_tag(i, "global", 1, Omega_h::GOs(mesh->nents(i), 0, 1));
    }
  }
  mesh->add_tags_for_ctrlPts();
  mesh->set_tag_for_ctrlPts(0, Reals(vtxPt_coords));
  mesh->set_tag_for_ctrlPts(1, Reals(edgePt_coords));
  mesh->set_tag_for_ctrlPts(2, Reals(facePt_coords));
  fprintf(stderr, "ok6\n");

  return;
}

void build_quartic_curveVtk(Mesh* mesh, Mesh* curveVtk_mesh,
                            LO n_sample_pts) {
  auto nface = mesh->nfaces();
  auto coords_h = HostRead<Real>(mesh->coords());
  auto ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(1));
  auto face_ctrlPts_h = HostRead<Real>(mesh->get_ctrlPts(2));
  auto fv2v_h = HostRead<LO>(mesh->ask_down(2, 0).ab2b);
  auto fe2e_h = HostRead<LO>(mesh->get_adj(2, 1).ab2b);
  auto ev2v_h = HostRead<LO>(mesh->get_adj(1, 0).ab2b);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);
  auto pts_per_face = mesh->n_internal_ctrlPts(2);

  I8 const dim = mesh->dim();
  Real xi_start = 0.0;
  Real xi_end = 1.0;
  Real delta_xi = (xi_end - xi_start)/(n_sample_pts - 1);
  double xi[n_sample_pts][n_sample_pts][FACE] = {0.0};
  LO count_curveVtk_mesh_vtx_perTri = 0;

  for (LO i = 0; i < n_sample_pts; ++i) {
    for (LO j = 0; j < n_sample_pts - i; ++j) {
      xi[i][j][0] = i*delta_xi;
      xi[i][j][1] = j*delta_xi;
      ++count_curveVtk_mesh_vtx_perTri;
    }
  }

  HostWrite<Real> host_coords(count_curveVtk_mesh_vtx_perTri*nface*dim);
  LO faces_perTri = (n_sample_pts - 1)*(n_sample_pts - 1);
  LO curveVtk_mesh_nface = faces_perTri*nface;
  HostWrite<LO> host_fv2v(curveVtk_mesh_nface*3);
  std::vector<int> face_vertices[1];
  
  LO count_curveVtk_mesh_vtx = 0;
  for (LO face = 0; face < nface; ++face) {
    auto v0 = fv2v_h[face*3];
    auto v1 = fv2v_h[face*3 + 1];
    auto v2 = fv2v_h[face*3 + 2];
    auto e0 = fe2e_h[face*3];
    auto e1 = fe2e_h[face*3 + 1];
    auto e2 = fe2e_h[face*3 + 2];

    I8 e0_flip = -1;
    I8 e1_flip = -1;
    I8 e2_flip = -1;
    auto e0v0 = ev2v_h[e0*2 + 0];
    auto e0v1 = ev2v_h[e0*2 + 1];
    auto e1v0 = ev2v_h[e1*2 + 0];
    auto e1v1 = ev2v_h[e1*2 + 1];
    auto e2v0 = ev2v_h[e2*2 + 0];
    auto e2v1 = ev2v_h[e2*2 + 1];
    if ((e0v0 == v1) && (e0v1 == v0)) {
      e0_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      e1_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      e2_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e2v0 == v2) && (e2v1 == v0));
    }

    if (dim == 3) {
      Real cx00 = coords_h[v0*dim + 0];
      Real cy00 = coords_h[v0*dim + 1];
      Real cz00 = coords_h[v0*dim + 2];
      Real cx40 = coords_h[v1*dim + 0];
      Real cy40 = coords_h[v1*dim + 1];
      Real cz40 = coords_h[v1*dim + 2];
      Real cx04 = coords_h[v2*dim + 0];
      Real cy04 = coords_h[v2*dim + 1];
      Real cz04 = coords_h[v2*dim + 2];

      Real cx10 = ctrlPts_h[e0*pts_per_edge*dim + 0];
      Real cy10 = ctrlPts_h[e0*pts_per_edge*dim + 1];
      Real cz10 = ctrlPts_h[e0*pts_per_edge*dim + 2];
      Real cx20 = ctrlPts_h[e0*pts_per_edge*dim + dim + 0];
      Real cy20 = ctrlPts_h[e0*pts_per_edge*dim + dim + 1];
      Real cz20 = ctrlPts_h[e0*pts_per_edge*dim + dim + 2];
      Real cx30 = ctrlPts_h[e0*pts_per_edge*dim + dim + dim + 0];
      Real cy30 = ctrlPts_h[e0*pts_per_edge*dim + dim + dim + 1];
      Real cz30 = ctrlPts_h[e0*pts_per_edge*dim + dim + dim + 2];
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        auto tempz = cz10;
        cx10 = cx30;
        cy10 = cy30;
        cz10 = cz30;
        cx30 = tempx;
        cy30 = tempy;
        cz30 = tempz;
      }

      Real cx31 = ctrlPts_h[e1*pts_per_edge*dim + 0];
      Real cy31 = ctrlPts_h[e1*pts_per_edge*dim + 1];
      Real cz31 = ctrlPts_h[e1*pts_per_edge*dim + 2];
      Real cx22 = ctrlPts_h[e1*pts_per_edge*dim + dim + 0];
      Real cy22 = ctrlPts_h[e1*pts_per_edge*dim + dim + 1];
      Real cz22 = ctrlPts_h[e1*pts_per_edge*dim + dim + 2];
      Real cx13 = ctrlPts_h[e1*pts_per_edge*dim + dim + dim + 0];
      Real cy13 = ctrlPts_h[e1*pts_per_edge*dim + dim + dim + 1];
      Real cz13 = ctrlPts_h[e1*pts_per_edge*dim + dim + dim + 2];
      if (e1_flip > 0) {
        auto tempx = cx31;
        auto tempy = cy31;
        auto tempz = cz31;
        cx31 = cx13;
        cy31 = cy13;
        cz31 = cz13;
        cx13 = tempx;
        cy13 = tempy;
        cz13 = tempz;
      }

      Real cx03 = ctrlPts_h[e2*pts_per_edge*dim + 0];
      Real cy03 = ctrlPts_h[e2*pts_per_edge*dim + 1];
      Real cz03 = ctrlPts_h[e2*pts_per_edge*dim + 2];
      Real cx02 = ctrlPts_h[e2*pts_per_edge*dim + dim + 0];
      Real cy02 = ctrlPts_h[e2*pts_per_edge*dim + dim + 1];
      Real cz02 = ctrlPts_h[e2*pts_per_edge*dim + dim + 2];
      Real cx01 = ctrlPts_h[e2*pts_per_edge*dim + dim + dim + 0];
      Real cy01 = ctrlPts_h[e2*pts_per_edge*dim + dim + dim + 1];
      Real cz01 = ctrlPts_h[e2*pts_per_edge*dim + dim + dim + 2];
      if (e2_flip > 0) {
        auto tempx = cx03;
        auto tempy = cy03;
        auto tempz = cz03;
        cx03 = cx01;
        cy03 = cy01;
        cz03 = cz01;
        cx01 = tempx;
        cy01 = tempy;
        cz01 = tempz;
      }

      Real cx11 = face_ctrlPts_h[face*pts_per_face*dim + 0];
      Real cy11 = face_ctrlPts_h[face*pts_per_face*dim + 1];
      Real cz11 = face_ctrlPts_h[face*pts_per_face*dim + 2];
      Real cx21 = face_ctrlPts_h[face*pts_per_face*dim + dim + 0];
      Real cy21 = face_ctrlPts_h[face*pts_per_face*dim + dim + 1];
      Real cz21 = face_ctrlPts_h[face*pts_per_face*dim + dim + 2];
      Real cx12 = face_ctrlPts_h[face*pts_per_face*dim + dim + dim + 0];
      Real cy12 = face_ctrlPts_h[face*pts_per_face*dim + dim + dim + 1];
      Real cz12 = face_ctrlPts_h[face*pts_per_face*dim + dim + dim + 2];

      for (LO i = 0; i < n_sample_pts; ++i) {
        for (LO j = 0; j < n_sample_pts - i; ++j) {
          auto x_bezier = cx00*B00_quart(xi[i][j][0], xi[i][j][1]) +
            cx10*B10_quart(xi[i][j][0], xi[i][j][1]) +
            cx20*B20_quart(xi[i][j][0], xi[i][j][1]) +
            cx30*B30_quart(xi[i][j][0], xi[i][j][1]) +
            cx40*B40_quart(xi[i][j][0], xi[i][j][1]) +
            cx31*B31_quart(xi[i][j][0], xi[i][j][1]) +
            cx22*B22_quart(xi[i][j][0], xi[i][j][1]) +
            cx13*B13_quart(xi[i][j][0], xi[i][j][1]) +
            cx04*B04_quart(xi[i][j][0], xi[i][j][1]) +
            cx03*B03_quart(xi[i][j][0], xi[i][j][1]) +
            cx02*B02_quart(xi[i][j][0], xi[i][j][1]) +
            cx01*B01_quart(xi[i][j][0], xi[i][j][1]) +
            cx11*B11_quart(xi[i][j][0], xi[i][j][1]) +
            cx21*B21_quart(xi[i][j][0], xi[i][j][1]) +
            cx12*B12_quart(xi[i][j][0], xi[i][j][1]);
          auto y_bezier = cy00*B00_quart(xi[i][j][0], xi[i][j][1]) +
            cy10*B10_quart(xi[i][j][0], xi[i][j][1]) +
            cy20*B20_quart(xi[i][j][0], xi[i][j][1]) +
            cy30*B30_quart(xi[i][j][0], xi[i][j][1]) +
            cy40*B40_quart(xi[i][j][0], xi[i][j][1]) +
            cy31*B31_quart(xi[i][j][0], xi[i][j][1]) +
            cy22*B22_quart(xi[i][j][0], xi[i][j][1]) +
            cy13*B13_quart(xi[i][j][0], xi[i][j][1]) +
            cy04*B04_quart(xi[i][j][0], xi[i][j][1]) +
            cy03*B03_quart(xi[i][j][0], xi[i][j][1]) +
            cy02*B02_quart(xi[i][j][0], xi[i][j][1]) +
            cy01*B01_quart(xi[i][j][0], xi[i][j][1]) +
            cy11*B11_quart(xi[i][j][0], xi[i][j][1]) +
            cy21*B21_quart(xi[i][j][0], xi[i][j][1]) +
            cy12*B12_quart(xi[i][j][0], xi[i][j][1]);
          auto z_bezier = cz00*B00_quart(xi[i][j][0], xi[i][j][1]) +
            cz10*B10_quart(xi[i][j][0], xi[i][j][1]) +
            cz20*B20_quart(xi[i][j][0], xi[i][j][1]) +
            cz30*B30_quart(xi[i][j][0], xi[i][j][1]) +
            cz40*B40_quart(xi[i][j][0], xi[i][j][1]) +
            cz31*B31_quart(xi[i][j][0], xi[i][j][1]) +
            cz22*B22_quart(xi[i][j][0], xi[i][j][1]) +
            cz13*B13_quart(xi[i][j][0], xi[i][j][1]) +
            cz04*B04_quart(xi[i][j][0], xi[i][j][1]) +
            cz03*B03_quart(xi[i][j][0], xi[i][j][1]) +
            cz02*B02_quart(xi[i][j][0], xi[i][j][1]) +
            cz01*B01_quart(xi[i][j][0], xi[i][j][1]) +
            cz11*B11_quart(xi[i][j][0], xi[i][j][1]) +
            cz21*B21_quart(xi[i][j][0], xi[i][j][1]) +
            cz12*B12_quart(xi[i][j][0], xi[i][j][1]);

          host_coords[count_curveVtk_mesh_vtx*dim + 0] = x_bezier;
          host_coords[count_curveVtk_mesh_vtx*dim + 1] = y_bezier;
          host_coords[count_curveVtk_mesh_vtx*dim + 2] = z_bezier;

          if ((i < n_sample_pts - 1) && (j < n_sample_pts - i - 1)) {
            face_vertices[0].push_back(count_curveVtk_mesh_vtx);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + n_sample_pts-i);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
            if (i > 0) {
              face_vertices[0].push_back(count_curveVtk_mesh_vtx);
              face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
              face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1 - (n_sample_pts-(i-1)));
            }
          }

          ++count_curveVtk_mesh_vtx;
        }
      }
    }
    else {
      OMEGA_H_CHECK(dim == 2);

      Real cx00 = coords_h[v0*dim + 0];
      Real cy00 = coords_h[v0*dim + 1];
      Real cx40 = coords_h[v1*dim + 0];
      Real cy40 = coords_h[v1*dim + 1];
      Real cx04 = coords_h[v2*dim + 0];
      Real cy04 = coords_h[v2*dim + 1];

      Real cx10 = ctrlPts_h[e0*pts_per_edge*dim + 0];
      Real cy10 = ctrlPts_h[e0*pts_per_edge*dim + 1];
      Real cx20 = ctrlPts_h[e0*pts_per_edge*dim + dim + 0];
      Real cy20 = ctrlPts_h[e0*pts_per_edge*dim + dim + 1];
      Real cx30 = ctrlPts_h[e0*pts_per_edge*dim + dim + dim + 0];
      Real cy30 = ctrlPts_h[e0*pts_per_edge*dim + dim + dim + 1];
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        cx10 = cx30;
        cy10 = cy30;
        cx30 = tempx;
        cy30 = tempy;
      }

      Real cx31 = ctrlPts_h[e1*pts_per_edge*dim + 0];
      Real cy31 = ctrlPts_h[e1*pts_per_edge*dim + 1];
      Real cx22 = ctrlPts_h[e1*pts_per_edge*dim + dim + 0];
      Real cy22 = ctrlPts_h[e1*pts_per_edge*dim + dim + 1];
      Real cx13 = ctrlPts_h[e1*pts_per_edge*dim + dim + dim + 0];
      Real cy13 = ctrlPts_h[e1*pts_per_edge*dim + dim + dim + 1];
      if (e1_flip > 0) {
        auto tempx = cx31;
        auto tempy = cy31;
        cx31 = cx13;
        cy31 = cy13;
        cx13 = tempx;
        cy13 = tempy;
      }

      Real cx03 = ctrlPts_h[e2*pts_per_edge*dim + 0];
      Real cy03 = ctrlPts_h[e2*pts_per_edge*dim + 1];
      Real cx02 = ctrlPts_h[e2*pts_per_edge*dim + dim + 0];
      Real cy02 = ctrlPts_h[e2*pts_per_edge*dim + dim + 1];
      Real cx01 = ctrlPts_h[e2*pts_per_edge*dim + dim + dim + 0];
      Real cy01 = ctrlPts_h[e2*pts_per_edge*dim + dim + dim + 1];
      if (e2_flip > 0) {
        auto tempx = cx03;
        auto tempy = cy03;
        cx03 = cx01;
        cy03 = cy01;
        cx01 = tempx;
        cy01 = tempy;
      }

      Real cx11 = face_ctrlPts_h[face*pts_per_face*dim + 0];
      Real cy11 = face_ctrlPts_h[face*pts_per_face*dim + 1];
      Real cx21 = face_ctrlPts_h[face*pts_per_face*dim + dim + 0];
      Real cy21 = face_ctrlPts_h[face*pts_per_face*dim + dim + 1];
      Real cx12 = face_ctrlPts_h[face*pts_per_face*dim + dim + dim + 0];
      Real cy12 = face_ctrlPts_h[face*pts_per_face*dim + dim + dim + 1];

      for (LO i = 0; i < n_sample_pts; ++i) {
        for (LO j = 0; j < n_sample_pts - i; ++j) {
          auto x_bezier = cx00*B00_quart(xi[i][j][0], xi[i][j][1]) +
            cx10*B10_quart(xi[i][j][0], xi[i][j][1]) +
            cx20*B20_quart(xi[i][j][0], xi[i][j][1]) +
            cx30*B30_quart(xi[i][j][0], xi[i][j][1]) +
            cx40*B40_quart(xi[i][j][0], xi[i][j][1]) +
            cx31*B31_quart(xi[i][j][0], xi[i][j][1]) +
            cx22*B22_quart(xi[i][j][0], xi[i][j][1]) +
            cx13*B13_quart(xi[i][j][0], xi[i][j][1]) +
            cx04*B04_quart(xi[i][j][0], xi[i][j][1]) +
            cx03*B03_quart(xi[i][j][0], xi[i][j][1]) +
            cx02*B02_quart(xi[i][j][0], xi[i][j][1]) +
            cx01*B01_quart(xi[i][j][0], xi[i][j][1]) +
            cx11*B11_quart(xi[i][j][0], xi[i][j][1]) +
            cx21*B21_quart(xi[i][j][0], xi[i][j][1]) +
            cx12*B12_quart(xi[i][j][0], xi[i][j][1]);
          auto y_bezier = cy00*B00_quart(xi[i][j][0], xi[i][j][1]) +
            cy10*B10_quart(xi[i][j][0], xi[i][j][1]) +
            cy20*B20_quart(xi[i][j][0], xi[i][j][1]) +
            cy30*B30_quart(xi[i][j][0], xi[i][j][1]) +
            cy40*B40_quart(xi[i][j][0], xi[i][j][1]) +
            cy31*B31_quart(xi[i][j][0], xi[i][j][1]) +
            cy22*B22_quart(xi[i][j][0], xi[i][j][1]) +
            cy13*B13_quart(xi[i][j][0], xi[i][j][1]) +
            cy04*B04_quart(xi[i][j][0], xi[i][j][1]) +
            cy03*B03_quart(xi[i][j][0], xi[i][j][1]) +
            cy02*B02_quart(xi[i][j][0], xi[i][j][1]) +
            cy01*B01_quart(xi[i][j][0], xi[i][j][1]) +
            cy11*B11_quart(xi[i][j][0], xi[i][j][1]) +
            cy21*B21_quart(xi[i][j][0], xi[i][j][1]) +
            cy12*B12_quart(xi[i][j][0], xi[i][j][1]);
          
          host_coords[count_curveVtk_mesh_vtx*dim + 0] = x_bezier;
          host_coords[count_curveVtk_mesh_vtx*dim + 1] = y_bezier;

          if ((i < n_sample_pts - 1) && (j < n_sample_pts - i - 1)) {
            face_vertices[0].push_back(count_curveVtk_mesh_vtx);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + n_sample_pts-i);
            face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
            if (i > 0) {
              face_vertices[0].push_back(count_curveVtk_mesh_vtx);
              face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1);
              face_vertices[0].push_back(count_curveVtk_mesh_vtx + 1 - (n_sample_pts-(i-1)));
            }
          }

          ++count_curveVtk_mesh_vtx;
        }
      }
    }
  }

  for (int i = 0; i < curveVtk_mesh_nface*3; ++i) {
    host_fv2v[i] = face_vertices[0][static_cast<std::size_t>(i)];
  }

  curveVtk_mesh->set_parting(OMEGA_H_ELEM_BASED);
  curveVtk_mesh->set_dim(dim);
  curveVtk_mesh->set_family(OMEGA_H_SIMPLEX);
  curveVtk_mesh->set_verts(count_curveVtk_mesh_vtx_perTri*nface);
  curveVtk_mesh->add_coords(Reals(host_coords.write()));
  curveVtk_mesh->set_down(2, 0, LOs(host_fv2v.write()));

  return;
}

}  // end namespace Omega_h
