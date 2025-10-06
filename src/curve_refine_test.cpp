#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_curve_validity_3d.hpp>
#include <Omega_h_refine.hpp>

using namespace Omega_h;

void test_adapt_inclusion(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read(
      "/lore/joshia5/Meshes/curved/inclusion_3p_sizes_2.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  //mesh.sync_tag(0, "metric");
  //mesh.sync_tag(0, "target_metric");
  vtk::FullWriter writer;

  /*
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "../omega_h/meshes/box_circleCut_4k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "../omega_h/meshes/box_circleCut_4k_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  */

  auto opts = AdaptOpts(&mesh);
  opts.should_swap = false;
  opts.should_coarsen = 0;
  opts.should_coarsen_slivers = 0;
  opts.should_filter_invalids = 0;
  opts.check_crv_qual = false;
  opts.verbosity = EXTRA_STATS;
  fprintf(stderr, "initial mesh size %d\n", mesh.nregions());
  I8 max_adapt_itr = 1;
  for (LO adapt_itr = 0; adapt_itr < max_adapt_itr; ++adapt_itr) {
    while (approach_metric(&mesh, opts) && mesh.nelems() < 10000) {
      adapt(&mesh, opts);
    }
  }
  auto qual = calc_crvQuality_3d(&mesh);
  /*
  writer = vtk::FullWriter(
      "../omega_h/meshes/boxCircle_aft.vtk", &mesh);
  writer.write();
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "../omega_h/meshes/box_circleCut_ref5k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "../omega_h/meshes/box_circleCut_ref5k_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  */
  return;
}

void test_adapt_rf(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read(
      "/lore/joshia5/Meshes/RF/assemble/v10_2rgn_11smallFeat_625k_p2.osh", comm);

  auto implied_metrics = get_implied_metrics(&mesh);
  LO const dim = 3;
  mesh.add_tag(VERT, "metric", symm_ncomps(dim), implied_metrics);
  mesh.add_tag<Real>(VERT, "target_metric", symm_ncomps(dim));
  
  auto target_metrics_w = Write<Real>
    (mesh.nverts() * symm_ncomps(dim));

  auto hd_hc = Write<Real>(mesh.nverts());
  Reals const coords = mesh.coords();
  auto ev2v = mesh.get_adj(1,0).ab2b;
  auto length_edg_w = Write<Real> (mesh.nedges());
  auto f1 = OMEGA_H_LAMBDA(LO e) {
    auto v0 = ev2v[e*2 + 0];
    auto v1 = ev2v[e*2 + 1];
    auto p0 = get_vector<dim>(coords, v0);
    auto p1 = get_vector<dim>(coords, v1);
    Real dist = 0.0;
    for (Int i = 0; i < dim; ++i) {
      dist += (p1[i] - p0[i])*(p1[i] - p0[i]);
    }
    dist = std::pow(dist, 0.5);
    length_edg_w[e] = dist;
  };
  parallel_for(mesh.nedges(), f1);
  mesh.add_tag(EDGE, "length_parent", 1, Reals(length_edg_w));
  ProjectFieldtoVertex (&mesh, "length_parent", 1);
  auto length_c = mesh.get_array<Real> (0, "length_parent");

  auto v_class_id = mesh.get_array<LO> (0, "class_id");

  auto f = OMEGA_H_LAMBDA(LO v) {
    auto h = Vector<dim>();
    for (Int i = 0; i < dim; ++i) {
      if (((v_class_id[v] == 5691) || (v_class_id[v] == 52842)) &&//52847
          ((coords[v*dim+1] > 0.5) && (coords[v*dim+2] > -0.2))) {
        printf("vtx %d on outer flux surf\n", v);
        h[i] = 0.5*length_c[v];
      }
      else {
        h[i] = length_c[v];
      }
      //h[i] = std::pow((error_des2/vtxError), 0.5)*length_c[v];
      //if ((v_class_id[v] == 190) ||  (v_class_id[v] == 186)) 
        //h[i] = length_c[v];
      //
      hd_hc[v] = h[i]/length_c[v];
    }
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh.nverts(), f);
  mesh.set_tag(VERT, "target_metric", Reals(target_metrics_w));
  mesh.add_tag(VERT, "hd_hc", 1, Reals(hd_hc));

/**//*set_target_metric<dim>(mesh, scale, pOmesh, error_des2);*/

  //mesh.set_parting(OMEGA_H_ELEM_BASED);
  //mesh.ask_lengths();
  //mesh.ask_qualities();
  //vtk::FullWriter writer;

  //printf("write mesh with size field\n");
  //binary::write("/lore/joshia5/Meshes/curved/inclusion_3p_sizes.osh", mesh);

  mesh.add_tag(0, "bezier_pts", 3, coords);
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  vtk::FullWriter writer;

  /*
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "../omega_h/meshes/box_circleCut_4k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "../omega_h/meshes/box_circleCut_4k_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  */

  auto opts = AdaptOpts(&mesh);
  opts.should_swap = false;
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = false;
  opts.should_filter_invalids = true;
  opts.check_crv_qual = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_allowed = 0.01;
  opts.max_length_allowed = 8.0;
  fprintf(stderr, "initial mesh size %d\n", mesh.nregions());
  I8 max_adapt_itr = 1;
  for (LO adapt_itr = 0; adapt_itr < max_adapt_itr; ++adapt_itr) {
    while (approach_metric(&mesh, opts)) {
    //while (approach_metric(&mesh, opts) && mesh.nelems() < 400000) 
      //approach_metric(&mesh, opts);
      adapt(&mesh, opts);
    }
  }
  //auto qual = calc_crvQuality_3d(&mesh);
  //writer = vtk::FullWriter(
    //  "../omega_h/meshes/boxCircle_aft.vtk", &mesh);
  //writer.write();
  
  printf("mesh size v e f r {%d,%d,%d,%d}\n", mesh.nverts(), mesh.nedges(),
      mesh.nfaces(), mesh.nelems());
  binary::write("/lore/joshia5/Meshes/RF/assemble/625k_ref1p5mil.osh",
      &mesh);
  /*
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh,4);
  std::string vtuPath = "/lore/joshia5/Meshes/RF/625k_ref_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh,4);
  vtuPath = "/lore/joshia5/Meshes/RF/625k_ref_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  */

  return; 
}

void uniref_p2_rf(Library *lib) {
  // uniform refinement
  auto comm = lib->world();

  auto mesh = binary::read(
      //"/lore/joshia5/Meshes/RF/testBall.osh", comm);
      "/lore/joshia5/Meshes/RF/assemble/v10_2rgn_12smallFeat_108k_p2.osh", comm);
      //"/lore/joshia5/Meshes/RF/assemble/v10_2rgn_12smallFeat_129k_p2.osh", comm);
  printf("mesh size v e f r {%d,%d,%d,%d}\n", mesh.nverts(), mesh.nedges(),
      mesh.nfaces(), mesh.nelems());
  Reals const coords = mesh.coords();
  if (!mesh.has_tag(0, "bezier_pts")) {
    mesh.add_tag(0, "bezier_pts", 3, coords);
  }
  LO const dim = 3;
  auto implied_metrics = get_implied_metrics(&mesh);
  mesh.add_tag(VERT, "metric", symm_ncomps(dim), implied_metrics);
  mesh.add_tag<Real>(VERT, "target_metric", symm_ncomps(dim));
  /*
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_quadratic_wireframe_3d(&mesh, &wireframe_mesh,3);
  std::string vtuPath = "/lore/joshia5/Meshes/RF/110k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_quadratic_curveVtk_3d(&mesh, &curveVtk_mesh,3);
  vtuPath = "/lore/joshia5/Meshes/RF/110k_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  */
  
  auto target_metrics_w = Write<Real>
    (mesh.nverts() * symm_ncomps(dim));

  auto hd_hc = Write<Real>(mesh.nverts());
  auto ev2v = mesh.get_adj(1,0).ab2b;
  auto length_edg_w = Write<Real> (mesh.nedges());
  auto f1 = OMEGA_H_LAMBDA(LO e) {
    auto v0 = ev2v[e*2 + 0];
    auto v1 = ev2v[e*2 + 1];
    auto p0 = get_vector<dim>(coords, v0);
    auto p1 = get_vector<dim>(coords, v1);
    Real dist = 0.0;
    for (Int i = 0; i < dim; ++i) {
      dist += (p1[i] - p0[i])*(p1[i] - p0[i]);
    }
    dist = std::pow(dist, 0.5);
    length_edg_w[e] = dist;
  };
  parallel_for(mesh.nedges(), f1);
  mesh.add_tag(EDGE, "length_parent", 1, Reals(length_edg_w));
  ProjectFieldtoVertex (&mesh, "length_parent", 1);
  auto length_c = mesh.get_array<Real> (0, "length_parent");

  auto f = OMEGA_H_LAMBDA(LO v) {
    auto h = Vector<dim>();
    for (Int i = 0; i < dim; ++i) {
      h[i] = 0.5*length_c[v];
      hd_hc[v] = h[i]/length_c[v];
    }
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh.nverts(), f);
  mesh.set_tag(VERT, "target_metric", Reals(target_metrics_w));
  mesh.add_tag(VERT, "hd_hc", 1, Reals(hd_hc));

  //mesh.set_parting(OMEGA_H_ELEM_BASED);
  //mesh.ask_lengths();
  //mesh.ask_qualities();
  //vtk::FullWriter writer;
  /*
  */

  calc_quad_ctrlPts_from_interpPts(&mesh);
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  /*
  vtk::FullWriter writer;
  writer = vtk::FullWriter("test_uniref_110k.vtk", &mesh);
  writer.write();
  */

  auto opts = AdaptOpts(&mesh);
  opts.should_swap = false;
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = false;
  opts.should_filter_invalids = true;
  opts.check_crv_qual = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_allowed = 0.001;
  opts.min_quality_desired = 0.001;
  opts.max_length_allowed = 6.0;
  opts.transfer_cands = true;
  fprintf(stderr, "initial mesh size %d\n", mesh.nregions());

  auto edge_is_cand = Bytes(mesh.nents(1), 1, "edge_is_cand");
  mesh.add_tag(EDGE, "candidate", 1, edge_is_cand);
  //opts.xfer_opts.type_map["candidate"] = OMEGA_H_INHERIT;
  while(refine(&mesh, opts));
  printf("refining...\n");
  /*
  I8 max_adapt_itr = 1;
  for (LO adapt_itr = 0; adapt_itr < max_adapt_itr; ++adapt_itr) {
    while (approach_metric(&mesh, opts)) {
      adapt(&mesh, opts);
    }
  }
  */
  
  printf("mesh size v e f r {%d,%d,%d,%d}\n", mesh.nverts(), mesh.nedges(),
      mesh.nfaces(), mesh.nelems());
  binary::write("/lore/joshia5/Meshes/RF/assemble/108kUniref862k_p2.osh",
      &mesh);
  /*
  wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_quadratic_wireframe_3d(&mesh, &wireframe_mesh,3);
  vtuPath = "/lore/joshia5/Meshes/RF/110kref_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_quadratic_curveVtk_3d(&mesh, &curveVtk_mesh,3);
  vtuPath = "/lore/joshia5/Meshes/RF/110kref_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  */

  return; 
}

void test_cubic_tet_quality(Library *lib) {
  auto mesh = Mesh(lib);
  auto comm = lib->world();
  binary::read(
      "/lore/joshia5/develop/mfem_omega/omega_h/meshes/Example_tet.osh",
      lib->world(), &mesh);
  mesh.set_curved(1);

  mesh.set_max_order(3);
  mesh.add_tags_for_ctrlPts();
  auto edge_nCtrlPts = mesh.n_internal_ctrlPts(1);
  auto coords = mesh.coords();
  auto ev2v = mesh.ask_down(1, 0).ab2b;
  auto fe2e = mesh.ask_down(2, 1).ab2b;
  auto rv2v = mesh.ask_down(3, 0).ab2b;
  auto re2e = mesh.ask_down(3, 1).ab2b;
  auto rf2f = mesh.ask_down(3, 2).ab2b;

  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  mesh.set_tag_for_ctrlPts(1, Reals({
                                     0.0, 1.0/3.0, 0.0,
                                     0.0, 2.0/3.0, 0.0,

                                     1.0/3.0, 0.0, 0.0,
                                     2.0/3.0, 0.0, 0.0,

                                     0.0, 0.0, 1.0/3.0,
                                     0.0, 0.0, 2.0/3.0,

                                     1.0/3.0, 2.0/3.0, 0.0,
                                     2.0/3.0, 1.0/3.0, 0.0,

                                     0.0, 1.0/3.0, 2.0/3.0,
                                     0.0, 2.0/3.0, 1.0/3.0,

                                     2.0/3.0, 0.0, 1.0/3.0,
                                     1.0/3.0, 0.0, 2.0/3.0
                                     }));
  mesh.set_tag_for_ctrlPts(2, Reals({
                                     1.0/3.0, 1.0/3.0, 0.0,
                                     1.0/3.0, 0.0, 1.0/3.0,
                                     1.0/3.0, 1.0/3.0, 1.0/3.0,
                                     0.0, 1.0/3.0, 1.0/3.0
                                     }));
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/tet_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/tet_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  add_implied_metric_tag(&mesh);
  auto qual = HostRead<Real>(calc_crvQuality_3d(&mesh));
  printf("regular tet quality %f\n", qual[0]);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  //test_adapt_inclusion(&lib);
  //test_cubic_tet_quality(&lib);
  //test_adapt_rf(&lib);
  uniref_p2_rf(&lib);

  return 0;
}
