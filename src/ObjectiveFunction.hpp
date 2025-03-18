#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <iostream>
#include <vector>

//#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"

#include <Eigen/Core>
#include <LBFGS.h>
#include <iostream>
//namespace Omega_h {}

class ObjectiveFunction
{
  private:
    Mesh* mesh;
    int n;
  public:
    ObjectiveFunction(int n_, Mesh* mesh_) : mesh(mesh_) {}
    double operator() (const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
      //forget about Mesh for now, just do it for a given set of input
      //doubles, look at the body of the askWorstQual and come in at the point
      //where an array or whatever containing coordinate locations of control
      //points is provided to the detJ calc function

      grad[1] = 20;
      grad[0]     = -2.0;

      return askWorstQuality_2d(mesh, LOs(mesh->nelems(),0,1), 2);

    }
/*

##############body  of ask worst qual
  auto fv2v = mesh->ask_down(2, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto faceCtrlPts = mesh->get_ctrlPts(2);
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto order = mesh->get_max_order();
  OMEGA_H_CHECK(order == 3);
  LO const nnew_tris = new_tris.size();
 
  auto Qs = mesh->ask_qualities();

  Write<Real> Q(nnew_tris, -1e-10);
  //LO const ntri_pts = 10;

  auto check_worstQual = OMEGA_H_LAMBDA (LO n) {
    Few<Real, 400> tri_pts;//ntri_pts*dim=20
    //Few<Real, 20> tri_pts;//ntri_pts*dim=20
    auto tri = new_tris[n];

    //query the tri's down verts's ctrl pts and store
    for (LO j = 0; j < 3; ++j) {//3 is tri2vert degree
      if (mesh_dim == 2) {
        auto p = get_vector<2>(vertCtrlPts, fv2v[tri*3 + j]);
        for (LO k = 0; k < mesh_dim; ++k) {
          tri_pts[j*mesh_dim + k] = p[k];
        }
      }
      else {
        OMEGA_H_CHECK (mesh_dim == 3);
        auto p = get_vector<3>(vertCtrlPts, fv2v[tri*3 + j]);
        for (LO k = 0; k < mesh_dim; ++k) {
          tri_pts[j*mesh_dim + k] = p[k];
        }
      }
    }

    //query the tri's down edge's ctrl pts and store

    auto v0 = fv2v[tri*3 + 0];
    auto v1 = fv2v[tri*3 + 1];
    auto v2 = fv2v[tri*3 + 2];
    auto e0 = fe2e[tri*3 + 0];
    auto e1 = fe2e[tri*3 + 1];
    auto e2 = fe2e[tri*3 + 2];
    auto e0v0 = ev2v[e0*2 + 0];
    auto e0v1 = ev2v[e0*2 + 1];
    auto e1v0 = ev2v[e1*2 + 0];
    auto e1v1 = ev2v[e1*2 + 1];
    auto e2v0 = ev2v[e2*2 + 0];
    auto e2v1 = ev2v[e2*2 + 1];
    auto flip = vector_3(-1, -1, -1);
    if ((e0v0 == v1) && (e0v1 == v0)) {
      flip[0] = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      flip[1] = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      flip[2] = 1;
    }

    for (LO j = 0; j < 3; ++j) {
      LO index = 3;
      if (flip[j] == -1) {
        for (I8 d = 0; d < mesh_dim; ++d) {
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + d];
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + mesh_dim + d];
        }
      }
      else {
        //for flipped edges
        OMEGA_H_CHECK (flip[j] == 1);
        for (I8 d = 0; d < mesh_dim; ++d) {
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + mesh_dim + d];
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + d];
        }
      }
    }

    //query the face's ctrl pt and store
    for (I8 d = 0; d < mesh_dim; ++d) {
      LO index = 9;
      tri_pts[index*mesh_dim + d] = faceCtrlPts[tri*mesh_dim + d];
    }

    //TODO change to template for mesh_dim
    auto nodes_det = getTriJacDetNodes<200, 2>(order, tri_pts);

    auto const minJ = calcMinJacDet(nodes_det, 3);
    auto const maxJ = calcMaxJacDet(nodes_det, 3);
    printf("tri %d qs %f, minJ %f, maxJ %f, Q %f\n",n,Qs[n], minJ, maxJ, Q[n]);
    Q[n] = std::pow((minJ/maxJ), 1./2.);
    //Q[n] = std::pow((minJ/maxJ), 1./2.)*Qs[n];
    //#######################################

  };
  parallel_for(nnew_tris, std::move(check_worstQual));

  printf("\nmin qual %f\n\n",get_min(Reals(Q)));
  return get_min(Reals(Q));

*/

    /*
    std::vector<double> getGrad(Mesh *mesh, const std::vector<double> &_x) {
      double h;
      std::vector<double> x = _x;
      double eps = this->getTol();
      std::vector<double> g;
      for (size_t i = 0; i < x.size(); i++) {
	h = abs(x[i]) > eps ? eps * abs(x[i]) : eps;

	// forward diff
	x[i] += h;
	double ff = this->getValue(mesh, x);
	x[i] -= h;

	// backward diff
	x[i] -= h;
	double fb = this->getValue(mesh, x);
	x[i] += h;

	g.push_back( (ff - fb) / 2./ h );
      }
      return(g);
    }
    */
    ~ObjectiveFunction(){};
};

#endif
