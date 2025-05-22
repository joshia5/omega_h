#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_atomics.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_scalar.hpp"

namespace Omega_h {

Real B0_quart(Real u) {
  return (1.0-u)*(1.0-u)*(1.0-u)*(1.0-u);
}

Real B1_quart(Real u) {
  return 4.0*u*(1.0-u)*(1.0-u)*(1.0-u);
}

Real B2_quart(Real u) {
  return 6.0*u*u*(1.0-u)*(1.0-u);
}

Real B3_quart(Real u) {
  return 4.0*u*u*u*(1.0-u);
}

Real B4_quart(Real u) {
  return u*u*u*u;
}

Real B00_quad(Real u, Real v) {
  return (1.0-u-v)*(1.0-u-v);
}

Real B10_quad(Real u, Real v) {
  return 2.0*u*(1.0-u-v);
}

Real B20_quad(Real u, Real v) {
  return u*u + 0.0*v;
}

Real B11_quad(Real u, Real v) {
  return 2.0*u*v;
}

Real B02_quad(Real u, Real v) {
  return v*v + 0.0*u;
}

Real B01_quad(Real u, Real v) {
  return 2.0*v*(1.0-u-v);
}

Real B00_cube(Real u, Real v) {
  return (1.0-u-v)*(1.0-u-v)*(1.0-u-v);
}

Real B10_cube(Real u, Real v) {
  return 3.0*u*(1.0-u-v)*(1.0-u-v);
}

Real B20_cube(Real u, Real v) {
  return 3.0*u*u*(1.0-u-v);
}

Real B30_cube(Real u, Real v) {
  return u*u*u + 0.0*v;
}

Real B21_cube(Real u, Real v) {
  return 3.0*u*u*v;
}

Real B12_cube(Real u, Real v) {
  return 3.0*u*v*v;
}

Real B03_cube(Real u, Real v) {
  return v*v*v + 0.0*u;
}

Real B02_cube(Real u, Real v) {
  return 3.0*v*v*(1.0-u-v);
}

Real B01_cube(Real u, Real v) {
  return 3.0*v*(1.0-u-v)*(1.0-u-v);
}

Real B11_cube(Real u, Real v) {
  return 6.0*u*v*(1.0-u-v);
}

Real B00_quart(Real u, Real v) {
  return (1.0-u-v)*(1.0-u-v)*(1.0-u-v)*(1.0-u-v);
}

Real B10_quart(Real u, Real v) {
  return 4.0*u*(1.0-u-v)*(1.0-u-v)*(1.0-u-v);
}

Real B20_quart(Real u, Real v) {
  return 6.0*u*u*(1.0-u-v)*(1.0-u-v);
}

Real B30_quart(Real u, Real v) {
  return 4.0*u*u*u*(1.0-u-v);
}

Real B40_quart(Real u, Real v) {
  return u*u*u*u + 0.0*v;
}

Real B31_quart(Real u, Real v) {
  return 4.0*u*u*u*v;
}

Real B22_quart(Real u, Real v) {
  return 6.0*u*u*v*v;
}

Real B13_quart(Real u, Real v) {
  return 4.0*u*v*v*v;
}

Real B04_quart(Real u, Real v) {
  return v*v*v*v + 0.0*u;
}

Real B03_quart(Real u, Real v) {
  return 4.0*v*v*v*(1.0-u-v);
}

Real B02_quart(Real u, Real v) {
  return 6.0*v*v*(1.0-u-v)*(1.0-u-v);
}

Real B01_quart(Real u, Real v) {
  return 4.0*v*(1.0-u-v)*(1.0-u-v)*(1.0-u-v);
}

Real B11_quart(Real u, Real v) {
  return 12.0*u*v*(1.0-u-v)*(1.0-u-v);
}

Real B21_quart(Real u, Real v) {
  return 12.0*u*u*v*(1.0-u-v);
}

Real B12_quart(Real u, Real v) {
  return 12.0*u*v*v*(1.0-u-v);
}

void calc_quad_ctrlPts_from_interpPts(Mesh *mesh) {
  auto coords = mesh->coords();
  auto interpPts = mesh->get_ctrlPts(1);
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto dim = mesh->dim();
  auto nedge = mesh->nedges();
  Real xi_1 = 0.5;
  Write<Real> new_pts(nedge*dim, 0.0);
  //auto const order = mesh->get_max_order();

  auto f = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];

    //TODO make this a templated fn like cal_pts_tmpl
    if (dim == 2) {
      Vector<2> c1;
      auto c0 = get_vector<2>(coords, v0);
      auto p1 = get_vector<2>(interpPts, i);
      auto c2 = get_vector<2>(coords, v1);
      for (Int j = 0; j < dim; ++j) {
        c1[j] = (p1[j] - B0_quad(xi_1)*c0[j] - B2_quad(xi_1)*c2[j])/B1_quad(xi_1);
      }
      set_vector(new_pts, i, c1);
    }
    else {
      OMEGA_H_CHECK(dim == 3);
      Vector<3> c1;
      auto c0 = get_vector<3>(coords, v0);
      auto p1 = get_vector<3>(interpPts, i);
      auto c2 = get_vector<3>(coords, v1);
      for (Int j = 0; j < dim; ++j) {
        c1[j] = (p1[j] - B0_quad(xi_1)*c0[j] - B2_quad(xi_1)*c2[j])/B1_quad(xi_1);
      }
      set_vector(new_pts, i, c1);
    }
  };
  parallel_for(nedge, std::move(f));

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
  return;
}

void elevate_curve_order_2to3(Mesh* mesh) {
  I8 new_order = 3;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;

  mesh->change_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  auto calc_edge_pts = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    if (dim == 3) {
      Vector<3> c1;
      Vector<3> c2;
      for (LO d = 0; d < dim; ++d) {
        c1[d] = (1.0/3.0)*coords[v0*dim + d] +
          (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
        c2[d] = (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
          (1.0/3.0)*coords[v1*dim + d];
        new_pts[i*n_new_pts*dim + d] = c1[d];
        new_pts[i*n_new_pts*dim + dim + d] = c2[d];
      }
    }
    else {
      OMEGA_H_CHECK (dim == 2);
      Vector<2> c1;
      Vector<2> c2;
      for (LO d = 0; d < dim; ++d) {
        c1[d] = (1.0/3.0)*coords[v0*dim + d] +
          (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
        c2[d] = (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
          (1.0/3.0)*coords[v1*dim + d];
        new_pts[i*n_new_pts*dim + d] = c1[d];
        new_pts[i*n_new_pts*dim + dim + d] = c2[d];
      }
    }

    // adding for a specific test case, rm after
    /*
    if (i == 12) {
      printf("edge %d ctrl pts ",i);
      for (LO d = 0; d < dim; ++d) {
        printf ("%f ", new_pts[i*n_new_pts*dim + d]);
      }
      printf("\n");
      for (LO d = 0; d < dim; ++d) {
        printf ("%f ", new_pts[i*n_new_pts*dim + dim + d]);
      }
      printf("\n");
      new_pts[i*n_new_pts*dim + dim-2] =-0.08;
      new_pts[i*n_new_pts*dim + dim-1] = 0.025;
    }
    if (i == 6) {
      printf("edge %d ctrl pts ",i);
      for (LO d = 0; d < dim; ++d) {
        printf ("%f ", new_pts[i*n_new_pts*dim + d]);
      }
      printf("\n");
      for (LO d = 0; d < dim; ++d) {
        printf ("%f ", new_pts[i*n_new_pts*dim + dim + d]);
      }
      printf("\n");
      new_pts[i*n_new_pts*dim + dim + dim-2] =-0.025;
      new_pts[i*n_new_pts*dim + dim + dim-1] = 0.065;
    }
    if (i == 7) {
      printf("edge %d ctrl pts ",i);
      for (LO d = 0; d < dim; ++d) {
        printf ("%f ", new_pts[i*n_new_pts*dim + d]);
      }
      printf("\n");
      for (LO d = 0; d < dim; ++d) {
        printf ("%f ", new_pts[i*n_new_pts*dim + dim + d]);
      }
      printf("\n");
      new_pts[i*n_new_pts*dim + dim + dim-2] = 0.06;
      new_pts[i*n_new_pts*dim + dim + dim-1] = 0.025;
    }
    if (i == 11) {
      printf("edge %d ctrl pts ",i);
      for (LO d = 0; d < dim; ++d) {
        printf ("%f ", new_pts[i*n_new_pts*dim + d]);
      }
      printf("\n");
      for (LO d = 0; d < dim; ++d) {
        printf ("%f ", new_pts[i*n_new_pts*dim + dim + d]);
      }
      printf("\n");
      new_pts[i*n_new_pts*dim + dim + dim-2] = 0.075;
      new_pts[i*n_new_pts*dim + dim + dim-1] = 0.025;
    }
    */
  };
  parallel_for(nedge, calc_edge_pts);
  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));

  auto nface = mesh->nfaces();
  n_new_pts = mesh->n_internal_ctrlPts(2);
  Write<Real> face_pts(nface*n_new_pts*dim, 0.0);
  if (!mesh->has_tag(2, "interp_pts")) {
    auto calc_face_pts = OMEGA_H_LAMBDA (LO i) {
      auto e0 = fe2e[i*3];
      auto e1 = fe2e[i*3 + 1];
      auto e2 = fe2e[i*3 + 2];
      for (LO d = 0; d < dim; ++d) {
        face_pts[i*n_new_pts*dim + d] =
          (1.0/3.0)*old_ctrl_pts[e0*old_n_ctrl_pts*dim + d] +
          (1.0/3.0)*old_ctrl_pts[e1*old_n_ctrl_pts*dim + d] +
          (1.0/3.0)*old_ctrl_pts[e2*old_n_ctrl_pts*dim + d];
      }
    };
    parallel_for(nface, calc_face_pts);
  }
  mesh->set_tag_for_ctrlPts(2, Reals(face_pts));

  return;
}

void var_order_2to1(Mesh* mesh) {
  //I8 new_order = 2;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;

  //mesh->change_max_order(new_order);

  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  //Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);

  Write<LO> n_quadratic_edges(1, 0);
  Write<LO> edge_order(nedge, 1);
  auto calc_edge_order = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    if (dim == 3) {
      Omega_h_fail("working on dim 2\n");
      /*
      Vector<3> c1;
      Vector<3> c2;
      for (LO d = 0; d < dim; ++d) {
        c1[d] = (1.0/3.0)*coords[v0*dim + d] +
          (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
        c2[d] = (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
          (1.0/3.0)*coords[v1*dim + d];
        new_pts[i*n_new_pts*dim + d] = c1[d];
        new_pts[i*n_new_pts*dim + dim + d] = c2[d];
      }
      */
    }
    else {
      OMEGA_H_CHECK (dim == 2);
      Vector<2> p0;
      Vector<2> p2;
      Vector<2> p1;
      Real m, c;
      for (LO d = 0; d < dim; ++d) {
        p0[d] = coords[v0*dim + d];
        p2[d] = coords[v1*dim + d];
        p1[d] = old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
      }
      m = (p2[1]-p0[1])/(p2[0]-p0[0]);
      c = p2[1] - (m*p2[0]);
      //printf("verif pt %1.15f eps %1.15f \n", (p0[1] - m*p0[0] - c) , EPSILON);
      OMEGA_H_CHECK((p0[1] - m*p0[0] - c) < EPSILON);

      if ((p1[1] - m*p1[0] - c) < EPSILON) {
        //printf("straight edge\n");
        edge_order[i] = 1;
      }
      else {
        //printf("quadratic edge\n");
        edge_order[i] = 2;
        atomic_increment(&n_quadratic_edges[0]);
      }
    }
  };
  parallel_for(nedge, calc_edge_order);
  printf("quadratic edges %d out of total %d\n",
      n_quadratic_edges[0], nedge);
  mesh->add_tag(1, "edge_order", 1, LOs(edge_order));
  //mesh->set_tag_for_ctrlPts(1, Reals(new_pts));

  auto nface = mesh->nfaces();
  Write<LO> face_order(nface, 1);
  auto calc_face_order = OMEGA_H_LAMBDA (LO i) {
    for (LO e = 0; e < 3; ++e) {
      auto adj_edge = fe2e[i*3 + e];
      if (edge_order[adj_edge] == 2) {
        face_order[i] = 2;
        break;
      }
    }
  };
  parallel_for(nface, calc_face_order);
  mesh->add_tag(2, "face_order", 1, LOs(face_order));

  if (dim == 3) {
    auto nrgn = mesh->nelems();
    auto rf2f = mesh->get_adj(3, 2).ab2b;
    Write<LO> rgn_order(nrgn, 1);
    auto calc_rgn_order = OMEGA_H_LAMBDA (LO i) {
      for (LO f = 0; f < 4; ++f) {
        auto adj_face = rf2f[i*4 + f];
        if (face_order[adj_face] == 2) {
          rgn_order[i] = 2;
          break;
        }
      }
    };
    parallel_for(nrgn, calc_rgn_order);
    mesh->add_tag(3, "rgn_order", 1, LOs(rgn_order));
  }

  return;
}

void elevate_curve_order_3to4(Mesh* mesh) {
  I8 new_order = 4;
  auto old_edge_ctrl_pts = mesh->get_ctrlPts(1);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);
  auto old_face_ctrlPts = mesh->get_ctrlPts(2);
  auto old_face_n_ctrl_pts = mesh->n_internal_ctrlPts(2);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;
  auto fv2v = mesh->ask_down(2, 0).ab2b;

  mesh->change_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  Write<Real> c1(dim, 0.0);
  Write<Real> c2(dim, 0.0);
  Write<Real> c3(dim, 0.0);
  auto calc_pts = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    for (LO d = 0; d < dim; ++d) {
      c1[d] = (1.0/4.0)*coords[v0*dim + d] +
              (3.0/4.0)*old_edge_ctrl_pts[i*pts_per_edge*dim + d];
      c2[d] = (1.0/2.0)*old_edge_ctrl_pts[i*pts_per_edge*dim + d] +
              (1.0/2.0)*old_edge_ctrl_pts[i*pts_per_edge*dim + 1*dim + d];
      c3[d] = (3.0/4.0)*old_edge_ctrl_pts[i*pts_per_edge*dim + 1*dim + d] +
              (1.0/4.0)*coords[v1*dim + d];
      new_pts[i*n_new_pts*dim + d] = c1[d];
      new_pts[i*n_new_pts*dim + 1*dim + d] = c2[d];
      new_pts[i*n_new_pts*dim + 2*dim + d] = c3[d];
    }
  };
  parallel_for(nedge, std::move(calc_pts));

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));

  auto nface = mesh->nfaces();
  auto n_new_face_pts = mesh->n_internal_ctrlPts(2);
  Write<Real> face_pts(nface*n_new_face_pts*dim, 0.0);
  fprintf(stderr, "elevated to quartic edges\n");

  auto calc_face_pts = OMEGA_H_LAMBDA (LO i) {
    auto e0 = fe2e[i*3];
    auto e1 = fe2e[i*3 + 1];
    auto e2 = fe2e[i*3 + 2];
    auto v0 = fv2v[i*3];
    auto v1 = fv2v[i*3 + 1];
    auto v2 = fv2v[i*3 + 2];

    I8 e0_flip = -1;
    I8 e1_flip = -1;
    I8 e2_flip = -1;

    auto e0v0 = ev2v[e0*2 + 0];
    auto e0v1 = ev2v[e0*2 + 1];
    auto e1v0 = ev2v[e1*2 + 0];
    auto e1v1 = ev2v[e1*2 + 1];
    auto e2v0 = ev2v[e2*2 + 0];
    auto e2v1 = ev2v[e2*2 + 1];
    if ((e0v0 == v1) && (e0v1 == v0)) {
      e0_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      e1_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      e2_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e2v0 == v2) && (e2v1 == v0));
    }

    if (dim == 2) {

      Real cx10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 0];
      Real cy10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 1];
      Real cx20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 0];
      Real cy20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 1];
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        cx10 = cx20;
        cy10 = cy20;
        cx20 = tempx;
        cy20 = tempy;
      }

      Real cx21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 0];
      Real cy21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 1];
      Real cx12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 0];
      Real cy12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 1];
      if (e1_flip > 0) {
        auto tempx = cx21;
        auto tempy = cy21;
        cx21 = cx12;
        cy21 = cy12;
        cx12 = tempx;
        cy12 = tempy;
      }

      Real cx02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 0];
      Real cy02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 1];
      Real cx01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 0];
      Real cy01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 1];
      if (e2_flip > 0) {
        auto tempx = cx02;
        auto tempy = cy02;
        cx02 = cx01;
        cy02 = cy01;
        cx01 = tempx;
        cy01 = tempy;
      }

      auto oldc11 = vector_2(old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 0],
          old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 1]);

      Vector<2> c11;
      Vector<2> c21;
      Vector<2> c12;

      c11[0] = (1.0/4.0)*(cx10 + cx01 + 2.0*oldc11[0]);
      c11[1] = (1.0/4.0)*(cy10 + cy01 + 2.0*oldc11[1]);
      c21[0] = (1.0/4.0)*(cx12 + cx02 + 2.0*oldc11[0]);
      c21[1] = (1.0/4.0)*(cy12 + cy02 + 2.0*oldc11[1]);
      c12[0] = (1.0/4.0)*(cx20 + cx21 + 2.0*oldc11[0]);
      c12[1] = (1.0/4.0)*(cy20 + cy21 + 2.0*oldc11[1]);

      for (LO d = 0; d < dim; ++d) {
        face_pts[i*n_new_face_pts*dim + d] = c11[d];
        face_pts[i*n_new_face_pts*dim + dim + d] = c21[d];
        face_pts[i*n_new_face_pts*dim + dim + dim + d] = c12[d];
      }
    }
    else {
      OMEGA_H_CHECK (dim == 3);

      Real cx10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 0];
      Real cy10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 1];
      Real cz10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 2];
      Real cx20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 0];
      Real cy20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 1];
      Real cz20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 2];
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        auto tempz = cz10;
        cx10 = cx20;
        cy10 = cy20;
        cz10 = cz20;
        cx20 = tempx;
        cy20 = tempy;
        cz20 = tempz;
      }

      Real cx21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 0];
      Real cy21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 1];
      Real cz21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 2];
      Real cx12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 0];
      Real cy12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 1];
      Real cz12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 2];
      if (e1_flip > 0) {
        auto tempx = cx21;
        auto tempy = cy21;
        auto tempz = cz21;
        cx21 = cx12;
        cy21 = cy12;
        cz21 = cz12;
        cx12 = tempx;
        cy12 = tempy;
        cz12 = tempz;
      }

      Real cx02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 0];
      Real cy02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 1];
      Real cz02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 2];
      Real cx01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 0];
      Real cy01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 1];
      Real cz01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 2];
      if (e2_flip > 0) {
        auto tempx = cx02;
        auto tempy = cy02;
        auto tempz = cz02;
        cx02 = cx01;
        cy02 = cy01;
        cz02 = cz01;
        cx01 = tempx;
        cy01 = tempy;
        cz01 = tempz;
      }

      auto oldc11 = vector_3(old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 0],
          old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 1],
          old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 2]);
      Vector<3> c11;
      Vector<3> c21;
      Vector<3> c12;

      c11[0] = (1.0/4.0)*(cx10 + cx01 + 2.0*oldc11[0]);
      c11[1] = (1.0/4.0)*(cy10 + cy01 + 2.0*oldc11[1]);
      c11[2] = (1.0/4.0)*(cz10 + cz01 + 2.0*oldc11[2]);
      c21[0] = (1.0/4.0)*(cx12 + cx02 + 2.0*oldc11[0]);
      c21[1] = (1.0/4.0)*(cy12 + cy02 + 2.0*oldc11[1]);
      c21[2] = (1.0/4.0)*(cz12 + cz02 + 2.0*oldc11[2]);
      c12[0] = (1.0/4.0)*(cx20 + cx21 + 2.0*oldc11[0]);
      c12[1] = (1.0/4.0)*(cy20 + cy21 + 2.0*oldc11[1]);
      c12[2] = (1.0/4.0)*(cz20 + cz21 + 2.0*oldc11[2]);

      for (LO d = 0; d < dim; ++d) {
        face_pts[i*n_new_face_pts*dim + d] = c11[d];
        face_pts[i*n_new_face_pts*dim + dim + d] = c21[d];
        face_pts[i*n_new_face_pts*dim + dim + dim + d] = c12[d];
      }
    }
  };
  parallel_for(nface, std::move(calc_face_pts));
  mesh->set_tag_for_ctrlPts(2, Reals(face_pts));

  //TODO calc pts inside region for (dim==3) mesh

  return;
}

void elevate_curve_order_4to5(Mesh* mesh) {
  I8 new_order = 5;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;

  mesh->change_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  Write<Real> c1(dim, 0.0);
  Write<Real> c2(dim, 0.0);
  Write<Real> c3(dim, 0.0);
  Write<Real> c4(dim, 0.0);
  auto calc_pts = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    for (LO d = 0; d < dim; ++d) {
      c1[d] = (1.0/5.0)*coords[v0*dim + d] +
              (4.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
      c2[d] = (2.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
              (3.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d];
      c3[d] = (3.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d] +
              (2.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 2*dim + d];
      c4[d] = (4.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 2*dim + d] +
              (1.0/5.0)*coords[v1*dim + d];
      new_pts[i*n_new_pts*dim + d] = c1[d];
      new_pts[i*n_new_pts*dim + 1*dim + d] = c2[d];
      new_pts[i*n_new_pts*dim + 2*dim + d] = c3[d];
      new_pts[i*n_new_pts*dim + 3*dim + d] = c4[d];
    }
  };
  parallel_for(nedge, std::move(calc_pts));

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
  
  auto nface = mesh->nfaces();
  auto n_new_face_pts = mesh->n_internal_ctrlPts(2);
  Write<Real> face_pts(nface*n_new_face_pts*dim, 0.0);
  fprintf(stderr, "elevated to quintic edges\n");

  /* TODO face pts for 3d
  auto calc_face_pts = OMEGA_H_LAMBDA (LO i) {
    auto e0 = fe2e[i*3];
    auto e1 = fe2e[i*3 + 1];
    auto e2 = fe2e[i*3 + 2];
    auto v0 = fv2v[i*3];
    auto v1 = fv2v[i*3 + 1];
    auto v2 = fv2v[i*3 + 2];

    I8 e0_flip = -1;
    I8 e1_flip = -1;
    I8 e2_flip = -1;

    auto e0v0 = ev2v[e0*2 + 0];
    auto e0v1 = ev2v[e0*2 + 1];
    auto e1v0 = ev2v[e1*2 + 0];
    auto e1v1 = ev2v[e1*2 + 1];
    auto e2v0 = ev2v[e2*2 + 0];
    auto e2v1 = ev2v[e2*2 + 1];
    if ((e0v0 == v1) && (e0v1 == v0)) {
      e0_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      e1_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      e2_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e2v0 == v2) && (e2v1 == v0));
    }

    if (dim == 2) {

      Real cx10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 0];
      Real cy10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 1];
      Real cx20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 0];
      Real cy20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 1];
      Real cx20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + dim + 0];
      Real cy20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + dim + 1];
      //TODO am stopping here  for today
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        cx10 = cx20;
        cy10 = cy20;
        cx20 = tempx;
        cy20 = tempy;
      }

      Real cx21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 0];
      Real cy21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 1];
      Real cx12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 0];
      Real cy12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 1];
      if (e1_flip > 0) {
        auto tempx = cx21;
        auto tempy = cy21;
        cx21 = cx12;
        cy21 = cy12;
        cx12 = tempx;
        cy12 = tempy;
      }

      Real cx02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 0];
      Real cy02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 1];
      Real cx01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 0];
      Real cy01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 1];
      if (e2_flip > 0) {
        auto tempx = cx02;
        auto tempy = cy02;
        cx02 = cx01;
        cy02 = cy01;
        cx01 = tempx;
        cy01 = tempy;
      }

      auto oldc11 = vector_2(old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 0],
          old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 1]);

      Vector<2> c11;
      Vector<2> c21;
      Vector<2> c12;

      c11[0] = (1.0/4.0)*(cx10 + cx01 + 2.0*oldc11[0]);
      c11[1] = (1.0/4.0)*(cy10 + cy01 + 2.0*oldc11[1]);
      c21[0] = (1.0/4.0)*(cx12 + cx02 + 2.0*oldc11[0]);
      c21[1] = (1.0/4.0)*(cy12 + cy02 + 2.0*oldc11[1]);
      c12[0] = (1.0/4.0)*(cx20 + cx21 + 2.0*oldc11[0]);
      c12[1] = (1.0/4.0)*(cy20 + cy21 + 2.0*oldc11[1]);

      for (LO d = 0; d < dim; ++d) {
        face_pts[i*n_new_face_pts*dim + d] = c11[d];
        face_pts[i*n_new_face_pts*dim + dim + d] = c21[d];
        face_pts[i*n_new_face_pts*dim + dim + dim + d] = c12[d];
      }
    }
    else {
      OMEGA_H_CHECK (dim == 3);

      Real cx10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 0];
      Real cy10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 1];
      Real cz10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 2];
      Real cx20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 0];
      Real cy20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 1];
      Real cz20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 2];
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        auto tempz = cz10;
        cx10 = cx20;
        cy10 = cy20;
        cz10 = cz20;
        cx20 = tempx;
        cy20 = tempy;
        cz20 = tempz;
      }

      Real cx21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 0];
      Real cy21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 1];
      Real cz21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 2];
      Real cx12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 0];
      Real cy12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 1];
      Real cz12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 2];
      if (e1_flip > 0) {
        auto tempx = cx21;
        auto tempy = cy21;
        auto tempz = cz21;
        cx21 = cx12;
        cy21 = cy12;
        cz21 = cz12;
        cx12 = tempx;
        cy12 = tempy;
        cz12 = tempz;
      }

      Real cx02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 0];
      Real cy02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 1];
      Real cz02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 2];
      Real cx01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 0];
      Real cy01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 1];
      Real cz01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 2];
      if (e2_flip > 0) {
        auto tempx = cx02;
        auto tempy = cy02;
        auto tempz = cz02;
        cx02 = cx01;
        cy02 = cy01;
        cz02 = cz01;
        cx01 = tempx;
        cy01 = tempy;
        cz01 = tempz;
      }

      auto oldc11 = vector_3(old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 0],
          old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 1],
          old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 2]);
      Vector<3> c11;
      Vector<3> c21;
      Vector<3> c12;

      c11[0] = (1.0/4.0)*(cx10 + cx01 + 2.0*oldc11[0]);
      c11[1] = (1.0/4.0)*(cy10 + cy01 + 2.0*oldc11[1]);
      c11[2] = (1.0/4.0)*(cz10 + cz01 + 2.0*oldc11[2]);
      c21[0] = (1.0/4.0)*(cx12 + cx02 + 2.0*oldc11[0]);
      c21[1] = (1.0/4.0)*(cy12 + cy02 + 2.0*oldc11[1]);
      c21[2] = (1.0/4.0)*(cz12 + cz02 + 2.0*oldc11[2]);
      c12[0] = (1.0/4.0)*(cx20 + cx21 + 2.0*oldc11[0]);
      c12[1] = (1.0/4.0)*(cy20 + cy21 + 2.0*oldc11[1]);
      c12[2] = (1.0/4.0)*(cz20 + cz21 + 2.0*oldc11[2]);

      for (LO d = 0; d < dim; ++d) {
        face_pts[i*n_new_face_pts*dim + d] = c11[d];
        face_pts[i*n_new_face_pts*dim + dim + d] = c21[d];
        face_pts[i*n_new_face_pts*dim + dim + dim + d] = c12[d];
      }
    }
  };
  parallel_for(nface, std::move(calc_face_pts));
  mesh->set_tag_for_ctrlPts(2, Reals(face_pts));
  */

  //TODO calc pts inside region for (dim==3) mesh


  return;
}

void elevate_curve_order_5to6(Mesh* mesh) {
  I8 new_order = 6;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;

  mesh->change_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  Write<Real> c1(dim, 0.0);
  Write<Real> c2(dim, 0.0);
  Write<Real> c3(dim, 0.0);
  Write<Real> c4(dim, 0.0);
  Write<Real> c5(dim, 0.0);
  auto calc_pts = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    for (LO d = 0; d < dim; ++d) {
      c1[d] = (1.0/6.0)*coords[v0*dim + d] +
              (5.0/6.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
      c2[d] = (1.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
              (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d];
      c3[d] = (1.0/2.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d] +
              (1.0/2.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 2*dim + d];
      c4[d] = (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 2*dim + d] +
              (1.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 3*dim + d];
      c5[d] = (5.0/6.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 3*dim + d] +
              (1.0/6.0)*coords[v1*dim + d];
      new_pts[i*n_new_pts*dim + d] = c1[d];
      new_pts[i*n_new_pts*dim + 1*dim + d] = c2[d];
      new_pts[i*n_new_pts*dim + 2*dim + d] = c3[d];
      new_pts[i*n_new_pts*dim + 3*dim + d] = c4[d];
      new_pts[i*n_new_pts*dim + 4*dim + d] = c5[d];
    }
  };
  parallel_for(nedge, std::move(calc_pts));

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
  return;
}

#define OMEGA_H_INST(T)
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

} // namespace Omega_h
