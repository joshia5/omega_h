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

/*
void test_collapse_boxCircle(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read(
      "../omega_h/meshes/box_circleCut_4k.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

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

  auto opts = AdaptOpts(&mesh);
  opts.should_swap = false;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(100)));
  fprintf(stderr, "initial mesh size %d\n", mesh.nregions());
  I8 max_adapt_itr = 4;
  for (LO adapt_itr = 0; adapt_itr < max_adapt_itr; ++adapt_itr) {
    coarsen_by_size(&mesh, opts);
    //coarsen_slivers(&mesh, opts);
  }
  mesh.ask_qualities();
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "../omega_h/meshes/boxCircle_aft.vtk", &mesh);
  writer.write();
  return;
}

void test_antenna(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("/lore/joshia5/Meshes/curved/antenna_6k.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/antenna_6k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/antenna_6k_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  
  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(100)));
  coarsen_by_size(&mesh, opts);
  mesh.ask_qualities();
  return;
}

void test_collapse_cubicSlab(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("/lore/joshia5/Meshes/curved/cubic_slab-case1.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/cubic_slab_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/cubic_slab_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
 
  auto opts = AdaptOpts(&mesh);
  opts.min_quality_allowed = 0.00001;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(10000)));
  coarsen_by_size(&mesh, opts);
  mesh.ask_qualities();
  return;
}

void test_sphere(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("/lore/joshia5/Meshes/curved/sphere_8.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/sphere_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/sphere_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
 
  //binary::write("/lore/joshia5/Meshes/curved/sphere_8_invalid.osh", &mesh);
  
  //auto opts = AdaptOpts(&mesh);
  //opts.min_quality_allowed = 0.00001;
  //mesh.add_tag<Real>(VERT, "metric", 1);
  //mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
  //      metric_eigenvalue_from_length(10000)));
  //coarsen_by_size(&mesh, opts);
  //mesh.ask_qualities();
  
  return;
}
*/

void test_swap_kova(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("../omega_h/meshes/KovaGeomSim-quadratic_123tet.osh", comm);
  printf("initial ntets %d\n", mesh.nregions());
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  vtk::FullWriter writer;
  writer = vtk::FullWriter("/lore/joshia5/Meshes/curved/kovaCoarsen_bef.vtk", &mesh);
  writer.write();
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);

  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
 
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 5);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/Kova_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 5);
  vtuPath = "/lore/joshia5/Meshes/curved/Kova_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

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
        //metric_eigenvalue_from_length(1)));
  fprintf(stderr, "start swapping\n");
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&mesh, opts);
  }
  fprintf(stderr, "finish swapping kova\n");
  check_validity_all_tet(&mesh);

  wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 5);
  vtuPath = "/lore/joshia5/Meshes/curved/Kova-swap_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 5);
  vtuPath = "/lore/joshia5/Meshes/curved/Kova-swap.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  return;
}

void test_disc_swap(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_100_order2.sms",
                            "/users/joshia5/Models/curved/disk_semi_geomsim_100.smd", comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/disc100_cubic_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/disc100_cubic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);

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

  wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 10);
  vtuPath =
    "/lore/joshia5/Meshes/curved/disc100-swap_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 10);
  vtuPath = "/lore/joshia5/Meshes/curved/disc100-swap.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris = checkValidity_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  fprintf(stderr, "finish swapping disc-2d\n");
  //mesh.ask_qualities();

  return;
}

void test_annulus3d_swap(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read("/lore/joshia5/Meshes/curved/annulus3d-24.sms",
                            "/lore/joshia5/Models/curved/annulus3d.smd", comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

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
  auto valid_tris_bef = checkValidity_3d(&mesh, LOs(mesh.nregions(), 0, 1));
  auto qual = calc_crvQuality_3d(&mesh);
  /*
  //auto quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&mesh, opts);
  }

  wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-8-swap_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-8-swap.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris_aft = checkValidity_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  //quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  fprintf(stderr, "finish swapping annulus-3d\n");
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-8-swap_full.vtk", &mesh);
  writer.write();
  */
  //mesh.ask_qualities();

  return;
}

void test_annulus120_swap(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read("/lore/joshia5/Meshes/curved/annulus-120d-4.sms",
                            "/lore/joshia5/Models/curved/annulus-120cut.smd", comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  printf("read %d elem mesh p3\n", mesh.nelems());
  auto const dim = mesh.dim();

  auto cubic_2tri = Mesh(lib);
  cubic_2tri.set_comm(comm);
  cubic_2tri.set_parting(OMEGA_H_ELEM_BASED);
  cubic_2tri.set_dim(dim);
  cubic_2tri.set_family(OMEGA_H_SIMPLEX);
  LO const numVtx=4;
  LO const numEdge=5;
  LO const numFace=2;
  cubic_2tri.set_verts(numVtx);

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
  cubic_2tri.add_coords(Reals(host_coords.write()));
  cubic_2tri.add_tag<ClassId>(0, "class_id", 1,
                           Read<LO>(class_id_0.write()));
  cubic_2tri.add_tag<I8>(0, "class_dim", 1,
                      Read<I8>(class_dim_0.write()));

  HostWrite<LO> ev2v(numEdge*2);
  ev2v = {0,1,  1,2,  2,3,  3,0, 1,3};
  cubic_2tri.set_ents(1, Adj(Read<LO>(ev2v.write())));
  HostWrite<LO> class_id_1(numEdge);
  class_id_1 = {13,7,17,15,2};
  HostWrite<I8> class_dim_1(numEdge);// all 0s
  class_dim_1 = {1,1,1,1,2};
  cubic_2tri.add_tag<ClassId>(1, "class_id", 1,
                           Read<LO>(class_id_1.write()));
  cubic_2tri.add_tag<I8>(1, "class_dim", 1,
                      Read<I8>(class_dim_1.write()));
  
  HostWrite<LO> fv2v(numFace*3);
  fv2v = {0,1,3,  2,3,1};
  //cubic_2tri.set_down(2, 0, LOs(fv2v.write()));//cant set if using set_ents
  Adj edge2vert = cubic_2tri.get_adj(1, 0);
  Adj vert2edge = cubic_2tri.ask_up(0, 1);
  Adj f2e = reflect_down(LOs(fv2v.write()), edge2vert.ab2b, vert2edge,
      OMEGA_H_SIMPLEX, 2, 1);
  cubic_2tri.set_ents(2, f2e);
  printf("ok5\n");
  HostWrite<LO> class_id_2(numFace);
  class_id_2 = {2,2};
  HostWrite<I8> class_dim_2(numFace);// all 0s
  class_dim_2 = {2,2};
  cubic_2tri.add_tag<ClassId>(2, "class_id", 1,
                           Read<LO>(class_id_2.write()));
  cubic_2tri.add_tag<I8>(2, "class_dim", 1,
                      Read<I8>(class_dim_2.write()));
  for (LO i = 0; i <= cubic_2tri.dim(); ++i) {
    if (!cubic_2tri.has_tag(i, "global")) {
      cubic_2tri.add_tag(i, "global", 1, Omega_h::GOs(cubic_2tri.nents(i), 0, 1));
    }
  }
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-120-2tri_full.vtk", &cubic_2tri);
  writer.write();
  cubic_2tri.set_curved(1);
  cubic_2tri.set_max_order(3);
  cubic_2tri.add_tags_for_ctrlPts();
  cubic_2tri.add_tag<Real>(0, "bezier_pts", cubic_2tri.dim(), cubic_2tri.coords());
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
  cubic_2tri.set_tag_for_ctrlPts(1, Reals(edgePt_coords.write()));
  cubic_2tri.set_tag_for_ctrlPts(2, Reals(facePt_coords.write()));
  
  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&cubic_2tri, &wireframe_mesh, 9);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-120_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&cubic_2tri, &cubic_curveVtk_mesh, 9);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-120_crvVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  binary::write("/lore/joshia5/Meshes/curved/annulus-120.osh", &cubic_2tri);
  /*
  */

  /*
  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-8_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-8.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);

  */
  auto opts = AdaptOpts(&cubic_2tri);
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = false;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_desired = 0.99;
  opts.min_quality_allowed = 0.98;
  cubic_2tri.add_tag<Real>(VERT, "metric", 1);
  cubic_2tri.set_tag(VERT, "metric", Reals(cubic_2tri.nverts(), 1));
  auto valid_tris_bef = checkValidity_2d(&cubic_2tri, LOs(cubic_2tri.nfaces(), 0, 1), 2);
  //auto quals = askQuality_2d(&cubic_2tri, LOs(cubic_2tri.nfaces(), 0, 1), 2);
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&cubic_2tri, opts);
  }

  wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&cubic_2tri, &wireframe_mesh, 20);
  vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-120-swap_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&cubic_2tri, &cubic_curveVtk_mesh, 20);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-120-swap.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris_aft = checkValidity_2d(&cubic_2tri, LOs(cubic_2tri.nfaces(), 0, 1), 2);
  //quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  binary::write("/lore/joshia5/Meshes/curved/annulus-120_swap.osh", &cubic_2tri);

  fprintf(stderr, "finish swapping annulus-2d\n");
  //vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-120-swap_full.vtk", &cubic_2tri);
  writer.write();
  /*
  */

  return;
}

void test_annulus_swap(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read("/lore/joshia5/Meshes/curved/annulus-8.sms",
                            "/lore/joshia5/Models/curved/annulus-8.smd", comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 20);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-8-modif_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 20);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-8-modif.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-8_full.vtk", &mesh);
  writer.write();

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
  auto quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  for (LO adapt_itr = 0; adapt_itr < 2; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&mesh, opts);
  }

  wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 20);
  vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-8-swap_2itr_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 20);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-8-swap_2itr.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris_aft = checkValidity_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  //quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  fprintf(stderr, "finish swapping annulus-2d\n");
  //vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-8-swap_2itr_full.vtk", &mesh);
  writer.write();
  binary::write("/lore/joshia5/Meshes/curved/annulus-8_swap_2itr.osh", &mesh);
  //mesh.ask_qualities();

  return;
}

void test_annulus_swap_p4(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read("/lore/joshia5/Meshes/curved/annulus-8.sms",
                            "/lore/joshia5/Models/curved/annulus-8.smd", comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  elevate_curve_order_3to4(&mesh);
  fprintf(stderr, "elevated to order 4\n");
  if (!mesh.has_tag(0, "bezier_pts")) {
    mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  }

  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_quartic_wireframe(&mesh, &wireframe_mesh, 5);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-8p4_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  fprintf(stderr, "generated order 4 wireframe\n");
  auto quartic_curveVtk_mesh = Mesh(lib);
  quartic_curveVtk_mesh.set_comm(comm);
  build_quartic_curveVtk(&mesh, &quartic_curveVtk_mesh, 5);
  vtuPath = "/lore/joshia5/Meshes/curved/annulusp4-8.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &quartic_curveVtk_mesh, 2);

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
  /*
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&mesh, opts);
  }

  wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-8-swap_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-8-swap.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  fprintf(stderr, "finish swapping annulus-2d\n");
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-8-swap_full.vtk", &mesh);
  writer.write();
  //mesh.ask_qualities();
  */

  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  //test_annulus120_swap(&lib); //2d 120 degree cut
  test_annulus_swap(&lib); //2d
  //test_annulus_swap_p4(&lib); //2d

  //test_annulus3d_swap(&lib); //2d

  //test_disc_swap(&lib); //2d
  //test_swap_kova(&lib); //3d
  
  return 0;
}
