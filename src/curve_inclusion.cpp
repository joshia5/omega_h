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
  opts.should_swap = 1;
  opts.should_coarsen = 1;
  opts.should_coarsen_slivers = 1;
  opts.should_filter_invalids = 0;
  opts.check_crv_qual = false;
  opts.verbosity = EXTRA_STATS;
  fprintf(stderr, "initial mesh %d tet\n", mesh.nregions());
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

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  test_adapt_inclusion(&lib);

  return 0;
}
