#include <Omega_h_adapt.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_curve_validity_3d.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>

using namespace Omega_h;

template <Int dim>
static void set_target_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto y = coords[v * dim + (dim - 2)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim; ++i)
    h[i] = 0.1+0.3*std::abs(y-0.05);
    //h[0] = 0.3; h[1] = 0.3; h[2] = 0.3;
    //for (Int i = 0; i < dim - 1; ++i) h[i] = 0.1;
    //h[dim - 2] = 0.001 + 0.198 * std::abs(y - 0.5);
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), f);
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
}

template <Int dim>
void run_case(Mesh* mesh, char const* vtk_path) {
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
  opts.length_histogram_max = 2.0;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  opts.should_smooth_snap = 0;
  opts.should_coarsen = 1;
  opts.should_swap = 0;
  opts.should_coarsen_slivers = 0;
  opts.check_crv_qual = 0;
  opts.min_quality_allowed = 0.1;
  opts.min_quality_desired = 0.25;
  Now t0 = now();
  auto desired_group_nelems = 2000;
  while (approach_metric(mesh, opts)) {
    auto nelems = mesh->nglobal_ents(mesh->dim());
    if (nelems < 8000) adapt(mesh, opts);
    else break;
    if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);
    /*
    auto wireframe_mesh = Mesh(comm->library());
    wireframe_mesh.set_comm(comm);
    build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 10);
    std::string vtuPath = 
      "/lore/joshia5/Meshes/curved/Kova-123tetwire_ansiov2.vtu";
    vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);

    auto curveVtk_mesh = Mesh(comm->library());
    curveVtk_mesh.set_comm(comm);
    build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 10);
    vtuPath = 
      "/lore/joshia5/Meshes/curved/Kova-123tetcurveVtk_anisov2.vtu";
    vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
    //if (vtk_path) writer.write();
    */
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
}

void test_cyl_grv(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read(
      "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/cyl_grv-5k.sms",
      "/users/joshia5/lore.scorec.rpi.edu/Models/curved/cyl_grv.smd",
      comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  /*
  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 8);
  std::string vtuPath =
    "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/cyl_grv-5k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &cubic_curveVtk_mesh, 8);
  vtuPath =
    "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/cyl_grv_5k.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris_aft = checkValidity_3d(&mesh, LOs(mesh.nregions(), 0, 1));

  while (-1) {
  }
  auto opts = AdaptOpts(&mesh);
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = true;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_desired = 0.99;
  opts.min_quality_allowed = 0.98;
  */
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(), 1));
  auto valid_tris_bef = checkValidity_3d(&mesh, LOs(mesh.nregions(), 0, 1));
  auto qual = calc_crvQuality_3d(&mesh);
  run_case<3>(&mesh, NULL);
  /*
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    adapt(&mesh, opts);
  }
  */

  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 8);
  std::string vtuPath =
    "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/cyl_grv-shock_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &cubic_curveVtk_mesh, 8);
  vtuPath =
    "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/cyl_grv_shock.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris_aft = checkValidity_3d(&mesh, LOs(mesh.nregions(), 0, 1));

  return;
}

void test_annulus3d(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read(
      "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/annulus3d-24.sms",
      "/users/joshia5/lore.scorec.rpi.edu/Models/curved/annulus3d.smd",
      comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  /*
  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/annulus3d_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus3d.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  */

  auto opts = AdaptOpts(&mesh);
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = true;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_desired = 0.99;
  opts.min_quality_allowed = 0.98;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(), 1));
  auto valid_tris_bef = checkValidity_3d(&mesh, LOs(mesh.nregions(), 0, 1));
  auto qual = calc_crvQuality_3d(&mesh);
  //auto quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    adapt(&mesh, opts);
  }

  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/annulus-3d-fs_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath =
    "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/annulus-3d-fs.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris_aft = checkValidity_3d(&mesh, LOs(mesh.nregions(), 0, 1));
  /*
  //quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-8-swap_full.vtk", &mesh);
  writer.write();
  */
  //mesh.ask_qualities();

  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  
  test_cyl_grv(&lib);
  //test_annulus3d(&lib);
  
  return 0;
}
