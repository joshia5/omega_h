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

#include <Eigen/Core>
#include <iostream>
#include <LBFGS.h>
#include <ObjectiveFunction.hpp>

using Eigen::VectorXf;
using Eigen::MatrixXf;
namespace lbfgs = LBFGSpp;

class Rosenbrock
{
  private:
    int n;
  public:
    Rosenbrock(int n_) : n(n_) {}
    float operator()(const VectorXf& x, VectorXf& grad)
    {
      float fx = 0.0;
      for(int i = 0; i < n; i += 2)
      {
        float t1 = 1.0 - x[i];
        float t2 = 10 * (x[i + 1] - x[i] * x[i]);
        grad[i + 1] = 20 * t2;
        grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
        fx += t1 * t1 + t2 * t2;
      }
      return fx;
    }
};

void test_rosenbrock() {

  const int n = 10;
  lbfgs::LBFGSParam<float> param;
  lbfgs::LBFGSSolver<float> solver(param);
  Rosenbrock fun(n);

  VectorXf x = VectorXf::Zero(n);
  float fx;
  int niter = solver.minimize(fun, x, fx);

  std::cout << niter << " iterations" << std::endl;
  std::cout << "x = \n" << x.transpose() << std::endl;
  std::cout << "f(x) = " << fx << std::endl;
  std::cout << "grad = " << solver.final_grad().transpose() << std::endl;
  std::cout << "||grad|| = " << solver.final_grad_norm() << std::endl;

  return;
}

void test_annuluswithEigen(Library *lib) {

  auto comm = lib->world();

  auto mesh_4tri = meshsim::read("/lore/joshia5/Meshes/curved/annulus-120d-4.sms",
                            "/lore/joshia5/Models/curved/annulus-120cut.smd", comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh_4tri);
  elevate_curve_order_2to3(&mesh_4tri);
  mesh_4tri.add_tag<Real>(0, "bezier_pts", mesh_4tri.dim(), mesh_4tri.coords());
  auto const dim = mesh_4tri.dim();

  auto mesh = Mesh(lib);
  mesh.set_comm(comm);
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  mesh.set_dim(dim);
  mesh.set_family(OMEGA_H_SIMPLEX);
  LO const numVtx=4;
  LO const numEdge=5;
  LO const numFace=2;
  mesh.set_verts(numVtx);

  HostWrite<Real> host_coords(numVtx*dim);
  host_coords = {
    0,0.25,
    0,0.5,
    -0.433013,-0.25,
    -0.216506,-0.125
                };
  HostWrite<LO> class_id_0(numVtx);
  class_id_0 = {12,11,16,14};
  HostWrite<I8> class_dim_0(numVtx);// all 0s
  class_dim_0 = {0,0,0,0};
  mesh.add_coords(Reals(host_coords.write()));
  mesh.add_tag<ClassId>(0, "class_id", 1,
                           Read<LO>(class_id_0.write()));
  mesh.add_tag<I8>(0, "class_dim", 1,
                      Read<I8>(class_dim_0.write()));

  HostWrite<LO> ev2v(numEdge*2);
  ev2v = {0,1,  1,2,  2,3,  3,0, 1,3};
  mesh.set_ents(1, Adj(Read<LO>(ev2v.write())));
  HostWrite<LO> class_id_1(numEdge);
  class_id_1 = {13,7,17,15,2};
  HostWrite<I8> class_dim_1(numEdge);// all 0s
  class_dim_1 = {1,1,1,1,2};
  mesh.add_tag<ClassId>(1, "class_id", 1,
                           Read<LO>(class_id_1.write()));
  mesh.add_tag<I8>(1, "class_dim", 1,
                      Read<I8>(class_dim_1.write()));
 
  HostWrite<LO> fv2v(numFace*3);
  fv2v = {0,1,3,  2,3,1};
  //mesh.set_down(2, 0, LOs(fv2v.write()));//cant set if using set_ents
  Adj edge2vert = mesh.get_adj(1, 0);
  Adj vert2edge = mesh.ask_up(0, 1);
  Adj f2e = reflect_down(LOs(fv2v.write()), edge2vert.ab2b, vert2edge,
      OMEGA_H_SIMPLEX, 2, 1);
  mesh.set_ents(2, f2e);
  printf("ok5\n");
  HostWrite<LO> class_id_2(numFace);
  class_id_2 = {2,2};
  HostWrite<I8> class_dim_2(numFace);// all 0s
  class_dim_2 = {2,2};
  mesh.add_tag<ClassId>(2, "class_id", 1,
                           Read<LO>(class_id_2.write()));
  mesh.add_tag<I8>(2, "class_dim", 1,
                      Read<I8>(class_dim_2.write()));
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-120-2tri_full.vtk", &mesh);
  writer.write();
  mesh.set_curved(1);
  mesh.set_max_order(3);
  mesh.add_tags_for_ctrlPts();
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  printf("making 2tri curve mesh\n");

  HostWrite<Real> edgePt_coords(numEdge*dim*2);
  // init straight sided
  for (LO i=0; i<numEdge; ++i) {
    LO const v0 = ev2v[i*2 + 0];
    LO const v1 = ev2v[i*2 + 1];
    Real len = 0.;
    for (LO d=0; d<dim; ++d) {
      len += std::pow((host_coords[v1*dim + d] - host_coords[v0*dim + d]),2);
    }
    len = std::pow(len, 0.5);
    for (LO d=0; d<dim; ++d) {
      edgePt_coords[i*dim*2 + d] = host_coords[v0*dim + d] + 
        (host_coords[v1*dim + d] - host_coords[v0*dim + d])/3.;
      edgePt_coords[i*dim*2 + dim + d] = host_coords[v0*dim + d] + 
        (host_coords[v1*dim + d] - host_coords[v0*dim + d])*2./3.;
    }
    Vector<2> c0,c3,p1,p2;
    c0[0] = host_coords[v0*dim + 0]; 
    c0[1] = host_coords[v0*dim + 1]; 
    c3[0] = host_coords[v1*dim + 0]; 
    c3[1] = host_coords[v1*dim + 1]; 
    if (i == 1) {
      /*
      p1[0] = 0.5*std::cos(PI*130./180.);
      p1[1] = 0.5*std::sin(PI*130./180.);
      p2[0] = 0.5*std::cos(PI*170./180.);
      p2[1] = 0.5*std::sin(PI*170./180.);

      auto const c1_c2 = curve_interpToCtrl_pts_2d(3, c0, c3, p1, p2);
      edgePt_coords[i*dim*2 + 0] = c1_c2[0];
      edgePt_coords[i*dim*2 + 1] = c1_c2[1];
      edgePt_coords[i*dim*2 + dim + 0] = c1_c2[2];
      edgePt_coords[i*dim*2 + dim + 1] = c1_c2[3];
      */

      edgePt_coords[i*dim*2 + 0] = 0.5*std::cos(PI*130./180.);
      edgePt_coords[i*dim*2 + 1] = 0.5*std::sin(PI*130./180.);

      edgePt_coords[i*dim*2 + dim + 0] = 0.5*std::cos(PI*170./180.);
      edgePt_coords[i*dim*2 + dim + 1] = 0.5*std::sin(PI*170./180.);
    }
    if (i == 3) {
      edgePt_coords[i*dim*2 + 0] = 0.25*std::cos(PI*170./180.);
      edgePt_coords[i*dim*2 + 1] = 0.25*std::sin(PI*170./180.);

      edgePt_coords[i*dim*2 + dim + 0] = 0.25*std::cos(PI*130./180.);
      edgePt_coords[i*dim*2 + dim + 1] = 0.25*std::sin(PI*130./180.);
    }
  }
  HostWrite<Real> facePt_coords(numFace*dim);
  for (LO i=0; i<numFace; ++i) {
    LO const v0 = fv2v[i*2 + 0];
    LO const v1 = fv2v[i*2 + 1];
    LO const v2 = fv2v[i*2 + 2];
    for (LO d=0; d<dim; ++d) {
      facePt_coords[i*dim + d] = (host_coords[v0*dim + d] + 
        host_coords[v1*dim + d] + host_coords[v2*dim + d])/3.;
    }
  }
  mesh.set_tag_for_ctrlPts(1, Reals(edgePt_coords.write()));
  mesh.set_tag_for_ctrlPts(2, Reals(facePt_coords.write()));

  auto opts = AdaptOpts(&mesh);
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = false;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_desired = 0.99;
  opts.min_quality_allowed = 0.98;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(), 1));
  auto valid_tris_bef = checkValidity_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  askWorstQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&mesh, opts);
  }
  askWorstQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  const int n = 8;//number of internal ctrl pts*dim // for 2d 2tri, 1edge case =8
  lbfgs::LBFGSParam<double> param;
  lbfgs::LBFGSSolver<double> solver(param);
  ObjectiveFunction foo(n, &mesh); //another class that has a constructor
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  double fx;
  solver.minimize(foo, x, fx);
  //returns func value
  //

  /*
  std::cout << niter << " iterations" << std::endl;
  std::cout << "x = \n" << x.transpose() << std::endl;
  std::cout << "f(x) = " << fx << std::endl;
  std::cout << "grad = " << solver.final_grad().transpose() << std::endl;
  std::cout << "||grad|| = " << solver.final_grad_norm() << std::endl;
  */

  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  test_rosenbrock();

  //test_annuluswithEigen(&lib);

  return 0;
}
