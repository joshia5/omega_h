#include "Omega_h_curve_validity_3d.hpp"
#include "Omega_h_for.hpp"

namespace Omega_h {

OMEGA_H_INLINE Real calcMinJacDet(Few<Real, 84> nodes) {
  Real minJ = 1e10;
  //for (LO i = 0; i < 84; ++i) // commenting to check only verts
  for (LO i = 0; i < 4; ++i)
    minJ = min2(minJ,nodes[i]);
  return minJ;
}

OMEGA_H_INLINE Real calcMaxJacDet(Few<Real, 84> nodes) {
  Real maxJ = -1e10;
  //for (LO i = 0; i < 84; ++i) // commenting to check only verts
  for (LO i = 0; i < 4; ++i)
    maxJ = max2(maxJ,nodes[i]);
  return maxJ;
}

Reals calc_crvQuality_3d(Mesh *mesh) {
  auto rv2v = mesh->ask_down(3, 0).ab2b;
  auto re2e = mesh->ask_down(3, 1).ab2b;
  auto rf2f = mesh->get_adj(3, 2).ab2b;
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto faceCtrlPts = mesh->get_ctrlPts(2);
  auto nnew_tets = mesh->nregions();
  auto order = mesh->get_max_order();
  OMEGA_H_CHECK(order == 3);
  auto qs = mesh->ask_qualities();

  Write<Real> Q(nnew_tets, -1e-10);

  auto calc_quality = OMEGA_H_LAMBDA (LO n) {
    Few<Real, 60> tet_pts = collect_tet_pts(order, n, ev2v, rv2v, vertCtrlPts
        , edgeCtrlPts, faceCtrlPts, re2e, rf2f);
    Few<Real, 84> nodes_det = getTetJacDetNodes<84>(3, tet_pts);

    auto const minJ = calcMinJacDet(nodes_det);
    auto const maxJ = calcMaxJacDet(nodes_det);
    if (minJ > 0.) {
      Q[n] = std::pow(std::pow((minJ/maxJ), 1./3.), 1.0);
      //Q[n] = std::pow(std::pow((minJ/maxJ), 1./3.)*qs[n], 1.0);
      if (Q[n] < 0.01) printf("low quality %f for element %d\n", Q[n], n);
    }
    else {
      Q[n] = 0.;
      //printf("tet %d qs %f, minJ %f, maxJ %f, Q %f\n",n,qs[n], minJ, maxJ, Q[n]);
    }
  };
  parallel_for(nnew_tets, std::move(calc_quality));

  return Reals(Q);
}

}
