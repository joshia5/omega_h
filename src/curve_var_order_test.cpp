#include <Omega_h_build.hpp>
#include <Omega_h_swap.hpp>
#include <Omega_h_library.hpp>
//#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_bezier_interp.hpp>
//#include <Omega_h_for.hpp>
#include <Omega_h_curve_coarsen.hpp>
#include <Omega_h_curve_validity_3d.hpp>

using namespace Omega_h;

void test_annulus(Library *lib) {

  auto comm = lib->world();

  auto mesh = meshsim::read("/lore/joshia5/Meshes/curved/annulus-3k.sms",
                            "/lore/joshia5/Models/curved/annulus.smd", comm);
 
  var_order_2to1(&mesh);
  //calc_quad_ctrlPts_from_interpPts(&mesh);
  //elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto rc = mesh.ask_revClass(1);
  LO nedges_rc = rc.ab2b.size();
  printf("num. mesh elms %d\n", nedges_rc);

  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  test_annulus(&lib);

  return 0;
}
