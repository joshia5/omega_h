#include <iostream>
#include <fstream>
#include <math.h>

#include <Omega_h_timer.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_refine.hpp>
#include <Omega_h_array_ops.hpp>

using namespace Omega_h;

template <Int dim>
static void set_target_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(dim));

  auto ev2v = mesh->get_adj(1,0).ab2b;
  auto length_c = Write<Real> (mesh->nedges());
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
    length_c[e] = dist;
  };  
  parallel_for(mesh->nedges(), f1);

  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    Real const x_max = 0.9;
    Real const scale = 2.0;
    h[0] = length_c[v];
    if (coords[v*dim + 0] > 0.25) {
      h[0] = length_c[v]*0.5;
      //h[0] = length_c[v]*((x_max - coords[v*dim + 0])/x_max/scale);
      //h[0] = 0.25*length_c[v];
    }
    //if (coords[v*dim + 0] > 0.25) h[0] = length_c[v]*0.25;
    if (coords[v*dim + 0] > 0.65) h[0] = length_c[v]*0.25;
    //if (coords[v*dim + 0] > 0.75) h[0] = length_c[v]*0.0625;
    printf("h %f\n", h[0]);
    h[1] = h[0]; h[2] = h[0];
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), f);
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
  printf("set target metric\n");
}

template <Int dim>
void run_case(Mesh* mesh, char const* vtk_path, CommPtr comm) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = get_implied_metrics(mesh);
  mesh->add_tag(VERT, "metric", symm_ncomps(dim), implied_metrics);
  mesh->add_tag<Real>(VERT, "target_metric", symm_ncomps(dim));
  set_target_metric<dim>(mesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->ask_lengths();
  mesh->ask_qualities();
  vtk::FullWriter writer;
  if (vtk_path) {
    writer = vtk::FullWriter(vtk_path, mesh);
    writer.write();
  }
  auto opts = AdaptOpts(mesh);
  opts.verbosity = EXTRA_STATS;
  opts.max_length_allowed = opts.max_length_desired * 8;
  opts.should_smooth_snap = 0;
  opts.should_coarsen = 0;
  opts.should_swap = 0;
  opts.should_coarsen_slivers = 0;
  opts.check_crv_qual = 0;
  opts.min_quality_allowed = 0.01;
  //opts.min_quality_desired = 0.2;
  Now t0 = now();
  printf("entering adapt\n");
  auto nelems = mesh->nglobal_ents(mesh->dim());
  LO const max_itr = 3;
  for (LO itr = 0; itr<max_itr; ++itr) {
    while (refine_by_size(mesh, opts)) {
    printf("refining itr \n");

    auto wireframe_mesh = Mesh(comm->library());
    wireframe_mesh.set_comm(comm);
    build_cubic_wireframe_3d(mesh, &wireframe_mesh, 10);
    std::string vtuPath = 
      "/lore/joshia5/Meshes/curved/coreRgn-refWire.vtu";
    vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);

    auto curveVtk_mesh = Mesh(comm->library());
    curveVtk_mesh.set_comm(comm);
    build_cubic_curveVtk_3d(mesh, &curveVtk_mesh, 10);
    vtuPath = 
      "/lore/joshia5/Meshes/curved/coreRgn-refCurveVtk.vtu";
    vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

    }
    //printf("appr. met. %d refining itr %d\n", approach_metric(mesh, opts), itr);
  }
  /*
  while (approach_metric(mesh, opts)) {
    nelems = mesh->nglobal_ents(mesh->dim());
    printf("adapting\n");
    //if (nelems < 4000) refine_by_size(mesh, opts);
    //if (nelems < 4000) adapt(mesh, opts);
    //else break;
    if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);
  }
  */
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
}

void test_edgeRgn(Library *lib) {
  auto comm = lib->world();
  auto mesh = binary::read("/lore/joshia5/Models/RF/edge_region_v102rgn12/edgeRgn_75kp2.osh", comm);
  printf("read mesh\n");
  if (!mesh.has_tag(0, "bezier_pts")) 
    mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  calc_quad_ctrlPts_from_interpPts(&mesh);

  elevate_curve_order_2to3(&mesh);
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 10);
  std::string vtuPath = 
    "/lore/joshia5/Meshes/curved/edgeRgn75k-Wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);

  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 10);
  vtuPath = 
    "/lore/joshia5/Meshes/curved/edgeRgn75k-CurveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

/*
  run_case<3>(&mesh, NULL, comm);

  AdaptOpts opts(&mesh);
  auto nelems = mesh.nglobal_ents(mesh.dim());
  Now t0 = now();
    if (!mesh.has_tag(0, "metric")) {
      add_implied_metric_tag(&mesh);
      adapt_refine(&mesh, opts);
      nelems = mesh.nglobal_ents(mesh.dim());
      std::cout << "mesh now has " << nelems << " total elements\n";
    }
    auto metrics = mesh.get_array<double>(0, "metric");
    metrics = multiply_each_by(metrics, 1.2);
    auto const metric_ncomps =
      divide_no_remainder(metrics.size(), mesh.nverts());
    mesh.add_tag(0, "metric", metric_ncomps, metrics);
    refine_by_size(&mesh, opts);
    //adapt_refine(&mesh, opts);
    nelems = mesh.nglobal_ents(mesh.dim());
    std::cout << "mesh now has " << nelems << " total elements\n";
  }
  Now t1 = now();
  std::cout << "total refine time: " << (t1 - t0) << " seconds\n";
  //vtk::write_parallel("/lore/joshia5/Meshes/curved/kova_refined.vtk", &mesh, 2);

  wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 10);
  vtuPath = 
    "/lore/joshia5/Meshes/curved/coreRgn-refWire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);

  curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 10);
  vtuPath = 
    "/lore/joshia5/Meshes/curved/coreRgn-refCurveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  */
  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  /*
  if (argc != 6) {
    Omega_h_fail(
      "a.out <2d_in_osh> <2d_out_vtk> <3d_in_model-geomsim> <3d_in_mesh> <3d_out_vtk>\n");
  };
  char const* path_2d = nullptr;
  char const* path_2d_vtk = nullptr;
  path_2d = argv[1];
  path_2d_vtk = argv[2];

  char const* path_3d_g = nullptr;
  char const* path_3d_m = nullptr;
  char const* path_3d_vtk = nullptr;
  path_3d_g = argv[3];
  path_3d_m = argv[4];
  path_3d_vtk = argv[5];

  */
  test_edgeRgn(&lib);

  return 0;
}
