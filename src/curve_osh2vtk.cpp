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

void rf(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read(
      "/users/joshia5/lore.scorec.rpi.edu/Meshes/RF/assemble/v102rgn12_643k-fscreen5mm_p2.osh", comm);
  printf("mesh size v e f r {%d,%d,%d,%d}\n", mesh.nverts(), mesh.nedges(),
      mesh.nfaces(), mesh.nelems());
  Reals const coords = mesh.coords();
  if (!mesh.has_tag(0, "bezier_pts")) {
    mesh.add_tag(0, "bezier_pts", 3, coords);
  }
  LO const dim = 3;
  
  calc_quad_ctrlPts_from_interpPts(&mesh);
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_quadratic_wireframe_3d(&mesh, &wireframe_mesh,10);
  std::string vtuPath = "/users/joshia5/lore.scorec.rpi.edu/Meshes/RF/interp_noise/643kfscreen5mm_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_quadratic_curveVtk_3d(&mesh, &curveVtk_mesh,10);
  vtuPath = "/users/joshia5/lore.scorec.rpi.edu/Meshes/RF/interp_noise/643kfscreen5mm_crvVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  return; 
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  rf(&lib);

  return 0;
}
