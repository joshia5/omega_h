#include <Omega_h_adapt.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_curve_validity_3d.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>

using namespace Omega_h;

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
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/annulus-3d-fs_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
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
  
  test_annulus3d(&lib);
  
  return 0;
}
