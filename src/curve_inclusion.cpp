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
      "/users/joshia5/lore.scorec.rpi.edu/Meshes/curved/inclusion_3p_sizes_2.osh"
      , comm);
                            
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

  auto old_mesh = mesh;

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
    //while (approach_metric(&mesh, opts) && mesh.nelems() < 10000) {
      //approach_metric(&mesh, opts, 1);
      adapt(&mesh, opts);
    //}
  }
  auto qual = calc_crvQuality_3d(&mesh);

  I8 const n_edge_pts = 2;
  I8 const dim = 3;
  auto const nedge = mesh.nedges();
  auto const ev2v = mesh.get_adj(1, 0).ab2b;
  auto const vertCtrlPts = mesh.get_ctrlPts(0);
  auto const old_edgeCtrlPts = mesh.get_ctrlPts(1);
  Write<Real> edge_ctrlPts(nedge*n_edge_pts*dim, INT8_MAX);
  auto const edge_gdim = mesh.get_array<I8>(1, "class_dim");
  auto const edge_gid = mesh.get_array<LO>(1, "class_id");
  auto stbdr_edge_points = OMEGA_H_LAMBDA(LO e) {
    if ((edge_gdim[e] == 2) && (edge_gid[e] != 190)) {
      auto const v0 = ev2v[e*2 + 0];
      auto const v1 = ev2v[e*2 + 1];
      for (LO j=0; j<dim; ++j) {
        edge_ctrlPts[e*n_edge_pts*dim + j] = vertCtrlPts[v0*dim + j] +
          (vertCtrlPts[v1*dim + j] - vertCtrlPts[v0*dim + j])/3.0;
        edge_ctrlPts[e*n_edge_pts*dim + dim + j] = vertCtrlPts[v0*dim + j] +
          (vertCtrlPts[v1*dim + j] - vertCtrlPts[v0*dim + j])*(2.0/3.0);
      }
    }
    else {
      for (LO j=0; j<dim*n_edge_pts; ++j) {
        edge_ctrlPts[e*n_edge_pts*dim + j] = old_edgeCtrlPts[e*n_edge_pts*dim + j];
      }
    }
  };
  parallel_for(nedge, std::move(stbdr_edge_points), "stbdr_edge_points");
  mesh.set_tag_for_ctrlPts(1, Reals(edge_ctrlPts));
  
  qual = calc_crvQuality_3d(&mesh);
  
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&old_mesh, &wireframe_mesh);
  std::string vtuPath = "../../../Meshes/curved/inclusion_adpt3k_wire"
    + std::to_string(comm->rank())
    +".vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&old_mesh, &curveVtk_mesh);
  vtuPath = "../../../Meshes/curved/inclusion_adpt3k_curveVtk"
    + std::to_string(comm->rank())
    +".vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
 
  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  test_adapt_inclusion(&lib);

  return 0;
}
