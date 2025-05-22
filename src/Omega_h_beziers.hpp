#ifndef OMEGA_H_BEZIERS_HPP
#define OMEGA_H_BEZIERS_HPP

#include "Omega_h_array.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_few.hpp"

namespace Omega_h {

OMEGA_H_INLINE LO binomial(int n, int i) {

  i = min2(n-i,i);

  if (i == 0)
    return 1;
  if (i == 1)
    return n;

  static LO const bn4[1] = {6};
  static LO const bn5[1] = {10};
  static LO const bn6[2] = {15,20};
  static LO const bn7[2] = {21,35};
  static LO const bn8[3] = {28,56,70};
  static LO const bn9[3] = {36,84,126};
  static LO const bn10[4] = {45,120,210,252};
  static LO const bn11[4] = {55,165,330,462};
  static LO const bn12[5] = {66,220,495,792,924};
  static LO const bn13[5] = {78,286,715,1287,1716};
  static LO const bn14[6] = {91,364,1001,2002,3003,3432};
  static LO const bn15[6] = {105,455,1365,3003,5005,6435};
  static LO const bn16[7] = {120,560,1820,4368,8008,11440,12870};
  static LO const bn17[7] = {136,680,2380,6188,12376,19448,24310};
  static LO const bn18[8] = {153,816,3060,8568,18564,31824,43758,48620};
  static LO const bn19[8] = {171,969,3876,11628,27132,50388,75582,92378};
  static LO const bn20[9] = {190,1140,4845,15504,38760,77520,125970,167960,
    184756};
  static LO const bn21[9] = {210,1330,5985,20349,54264,116280,203490,293930,
    352716};
  static LO const bn22[10] = {231,1540,7315,26334,74613,170544,319770,497420,
    646646,705432};
  static LO const bn23[10] = {253,1771,8855,33649,100947,245157,490314,817190,
    1144066,1352078};
  static LO const bn24[11] = {276,2024,10626,42504,134596,346104,735471,1307504,
    1961256,2496144,2704156};
  static LO const bn25[11] = {300,2300,12650,53130,177100,480700,1081575,2042975,
    3268760,4457400,5200300};
  static LO const bn26[12] = {325,2600,14950,65780,230230,657800,1562275,3124550,
    5311735,7726160,9657700,10400600};
  static LO const bn27[12] = {351,2925,17550,80730,296010,888030,2220075,4686825,
    8436285,13037895,17383860,20058300};
  static LO const bn28[13] = {378,3276,20475,98280,376740,1184040,3108105,6906900,
    13123110,21474180,30421755,37442160,40116600};
  static LO const bn29[13] = {406,3654,23751,118755,475020,1560780,4292145,10015005,
    20030010,34597290,51895935,67863915,77558760};
  static LO const bn30[14] = {435,4060,27405,142506,593775,2035800,5852925,14307150,
    30045015,54627300,86493225,119759850,145422675,155117520};
  static LO const bn31[14] = {465,4495,31465,169911,736281,2629575,7888725,20160075,
    44352165,84672315,141120525,206253075,265182525,300540195};
  static LO const bn32[15] = {496,4960,35960,201376,906192,3365856,10518300,28048800,
    64512240,129024480,225792840,347373600,471435600,565722720,601080390};
  static LO const bn33[15] = {528,5456,40920,237336,1107568,4272048,13884156,38567100,
    92561040,193536720,354817320,573166440,818809200,1037158320,1166803110};

  static LO const* const bnTable[34] = {0,0,0,0,bn4,bn5,bn6,bn7,bn8,
    bn9,bn10,bn11,bn12,bn13,bn14,bn15,bn16,bn17,bn18,bn19,bn20,bn21,bn22,bn23,
    bn24,bn25,bn26,bn27,bn28,bn29,bn30,bn31,bn32,bn33};

  return bnTable[n][i-2];
}

OMEGA_H_INLINE LO trinomial(LO n, LO i, LO j) {
  return binomial(n,i)*binomial(n-i,j);
}

OMEGA_H_INLINE LO quadnomial(LO n, LO i, LO j, LO k) {
  return binomial(n,i)*binomial(n-i,j)*binomial(n-i-j,k);
}

OMEGA_H_INLINE Real const_factorial(Int N) {
  switch (N) {
    case 0:
      return 1.0;
    case 1:
      return 1.0;
    case 2:
      return 2.0;
    case 3:
      return 6.0;
    case 4:
      return 24.0;
    case 5:
      return 120.0;
    case 6:
      return 720.0;
    case 7:
      return 5040.0;
    case 8:
      return 40320.0;
    case 9:
      return 362880.0;
    case 10:
      return 3628800.0;
  }
  return -1.0;
}

OMEGA_H_INLINE Real B0_quad(Real u) {
  return (1.0-u)*(1.0-u);
}

OMEGA_H_INLINE Real B1_quad(Real u) {
  return 2.0*u*(1.0-u);
}

OMEGA_H_INLINE Real B2_quad(Real u) {
  return u*u;
}

OMEGA_H_INLINE Real B0_cube(Real u) {
  return (1.0-u)*(1.0-u)*(1.0-u);
}

OMEGA_H_INLINE Real B1_cube(Real u) {
  return 3.0*u*(1.0-u)*(1.0-u);
}

OMEGA_H_INLINE Real B2_cube(Real u) {
  return 3.0*(1.0-u)*u*u;
}

OMEGA_H_INLINE Real B3_cube(Real u) {
  return u*u*u;
}

Real B0_quart(Real u);
Real B1_quart(Real u);
Real B2_quart(Real u);
Real B3_quart(Real u);
Real B4_quart(Real u);

Real B00_quad(Real u, Real v);
Real B10_quad(Real u, Real v);
Real B20_quad(Real u, Real v);
Real B11_quad(Real u, Real v);
Real B02_quad(Real u, Real v);
Real B01_quad(Real u, Real v);

Real B00_cube(Real u, Real v);
Real B10_cube(Real u, Real v);
Real B20_cube(Real u, Real v);
Real B30_cube(Real u, Real v);
Real B21_cube(Real u, Real v);
Real B12_cube(Real u, Real v);
Real B03_cube(Real u, Real v);
Real B02_cube(Real u, Real v);
Real B01_cube(Real u, Real v);
Real B11_cube(Real u, Real v);

Real B00_quart(Real u, Real v);
Real B10_quart(Real u, Real v);
Real B20_quart(Real u, Real v);
Real B30_quart(Real u, Real v);
Real B40_quart(Real u, Real v);
Real B31_quart(Real u, Real v);
Real B22_quart(Real u, Real v);
Real B13_quart(Real u, Real v);
Real B04_quart(Real u, Real v);
Real B03_quart(Real u, Real v);
Real B02_quart(Real u, Real v);
Real B01_quart(Real u, Real v);
Real B11_quart(Real u, Real v);
Real B21_quart(Real u, Real v);
Real B12_quart(Real u, Real v);

// Babushka chen points
// edge points
constexpr OMEGA_H_INLINE Real xi_1_quad() {
  return 0.5;
}

constexpr OMEGA_H_INLINE Real xi_1_cube() {
  return 0.2748043;
}
constexpr OMEGA_H_INLINE Real xi_2_cube() {
  return 0.7251957;
}

constexpr OMEGA_H_INLINE Real xi_1_quart() {
  return 0.1693976;
}
constexpr OMEGA_H_INLINE Real xi_2_quart() {
  return 0.5;
}
constexpr OMEGA_H_INLINE Real xi_3_quart() {
  return 0.8306024;
}

constexpr OMEGA_H_INLINE Real xi_1_quint() {
  return 0.1133573;
}
constexpr OMEGA_H_INLINE Real xi_2_quint() {
  return 0.3568239;
}
constexpr OMEGA_H_INLINE Real xi_3_quint() {
  return 0.6431761;
}
constexpr OMEGA_H_INLINE Real xi_4_quint() {
  return 0.8866427;
}

constexpr OMEGA_H_INLINE Real xi_1_hex() {
  return 0.0805979;
}
constexpr OMEGA_H_INLINE Real xi_2_hex() {
  return 0.2650895;
}
constexpr OMEGA_H_INLINE Real xi_3_hex() {
  return 0.5;
}
constexpr OMEGA_H_INLINE Real xi_4_hex() {
  return 0.7349105;
}
constexpr OMEGA_H_INLINE Real xi_5_hex() {
  return 0.9194021;
}

// face points
// order 3
OMEGA_H_INLINE Vector<2> xi_11_cube() {
  return vector_2(1.0/3.0, 1.0/3.0);
}
// order 4
OMEGA_H_INLINE Vector<2> xi_11_quart() {
  return vector_2(0.22088805, 0.22088805);
}
OMEGA_H_INLINE Vector<2> xi_21_quart() {
  return vector_2(0.5582239, 0.22088805);
}
OMEGA_H_INLINE Vector<2> xi_12_quart() {
  return vector_2(0.22088805, 0.5582239);
}
// order 6
OMEGA_H_INLINE Vector<2> xi_11_hex() {
  return vector_2(0.10971385,0.10971385);
}
OMEGA_H_INLINE Vector<2> xi_21_hex() {
  return vector_2(0.3157892,0.1256031);
}
OMEGA_H_INLINE Vector<2> xi_31_hex() {
  return vector_2(0.5586077,0.1256031);
}
OMEGA_H_INLINE Vector<2> xi_41_hex() {
  return vector_2(0.7805723,0.10971385);
}
OMEGA_H_INLINE Vector<2> xi_12_hex() {
  return vector_2(0.1256031,0.3157892);
}
OMEGA_H_INLINE Vector<2> xi_22_hex() {
  return vector_2(1./3.,1./3.);
}
OMEGA_H_INLINE Vector<2> xi_32_hex() {
  return vector_2(0.5586077,0.3157892);
}
OMEGA_H_INLINE Vector<2> xi_13_hex() {
  return vector_2(0.1256031,0.5586077);
}
OMEGA_H_INLINE Vector<2> xi_23_hex() {
  return vector_2(0.3157892,0.5586077);
}
OMEGA_H_INLINE Vector<2> xi_14_hex() {
  return vector_2(0.10971385,0.7805723);
}

// region points
// order 4
OMEGA_H_INLINE Vector<3> xi_111_quart() {
  return vector_3(0.25,0.25,0.25);
}
// order 6

/*
//TODO limititations for generalized order : 1. these fns will need to be
//dynamically sized, 2. lin eqn form and solve on gpu 
Reals curve_bezier_pts(LO const P) {
//OMEGA_H_INLINE Reals curve_bezier_pts(LO const P) {
  switch (P) {
    case 2: {
      return Reals({xi_1_quad()});
    }
    case 3: {
      return Reals({xi_1_cube(),xi_2_cube()});
    }
    case 4: {
      return Reals({xi_1_quart(),xi_2_quart(),xi_3_quart()});
    }
    case 5: {
      return Reals({xi_1_quint(),xi_2_quint(),xi_3_quint(),xi_4_quint()});
    }
    case 6: {
      return Reals({xi_1_hex(),xi_2_hex(),xi_3_hex(),xi_4_hex(),xi_5_hex()});
    }
  }
  OMEGA_H_NORETURN(Reals());
}

Reals triangle_bezier_pts(LO const P) {
//OMEGA_H_INLINE Reals triangle_bezier_pts(LO const P) {
  switch (P) {
    case 3: {
      return xi_11_cube();
    }
    case 4: {
      return Read<Real>(concat<Real>(concat<Real>(xi_11_quart(), xi_21_quart()), xi_12_quart()));
    }
    case 5: {
      return Reals({});
    }
    case 6: {
      return Reals({});
    }
  }
  OMEGA_H_NORETURN(Reals());
}

OMEGA_H_INLINE Reals tet_bezier_pts(LO const P) {
  switch (P) {
    case 4: {
      return xi_111_quart();
    }
    case 5: {
      return Reals({});
    }
    case 6: {
      return Reals({});
    }
  }
  OMEGA_H_NORETURN(Reals());
}
*/

OMEGA_H_INLINE Real
Bijk(LO const P, LO const i, LO const j, LO const k, Real const u,
     Real const v, Real const w) noexcept {
  LO const l = P - i - j - k;
  OMEGA_H_CHECK(l >= 0);
  Real const t = 1.0 - u - v - w;
  OMEGA_H_CHECK((t >= 0.0) && (t <= 1.0));
  Real resultant = 1.0;
  resultant = resultant*const_factorial(P);
  resultant = resultant*std::pow(u,i);
  resultant = resultant*std::pow(v,j);
  resultant = resultant*std::pow(w,k);
  resultant = resultant*std::pow(t,l);
  resultant = resultant/const_factorial(i);
  resultant = resultant/const_factorial(j);
  resultant = resultant/const_factorial(k);
  resultant = resultant/const_factorial(l);

  return resultant;
}

OMEGA_H_INLINE Real 
Bij(LO const P, LO const i, LO const j, Real const u, Real const v) noexcept {
  LO const k = P - i - j;
  OMEGA_H_CHECK((k >= 0) && (k <= P));
  Real const w = 1.0 - u - v;
  OMEGA_H_CHECK((w >= 0.0) && (w <= 1.0));
  Real resultant = 1.0;
  resultant = resultant*const_factorial(P);
  resultant = resultant*std::pow(u,i);
  resultant = resultant*std::pow(v,j);
  resultant = resultant*std::pow(w,k);
  resultant = resultant/const_factorial(i);
  resultant = resultant/const_factorial(j);
  resultant = resultant/const_factorial(k);

  return resultant;
}

OMEGA_H_INLINE Real
Bi(LO const P, LO const i, Real const u) noexcept {
  LO const j = P - i;
  OMEGA_H_CHECK((j >= 0) && (j <= P));
  Real const v = 1.0 - u;
  OMEGA_H_CHECK((v >= 0.0) && (v <= 1.0));
  Real resultant = 1.0;
  resultant = resultant*std::pow(u,i);
  resultant = resultant*std::pow(v,j);
  resultant = resultant*const_factorial(P);
  resultant = resultant/const_factorial(i);
  resultant = resultant/const_factorial(j);

  return resultant;
}

void elevate_curve_order_2to3(Mesh* mesh);
void elevate_curve_order_3to4(Mesh* mesh);
void elevate_curve_order_4to5(Mesh* mesh);
void elevate_curve_order_5to6(Mesh* mesh);

void var_order_2to1(Mesh* mesh);

void calc_quad_ctrlPts_from_interpPts(Mesh *mesh);

OMEGA_H_INLINE Few<Real,2> quadr_noKeyEdge_xi_values(
    LO old_vert, LO v0, LO v1, LO v2, LO old_edge, LO e0, LO e1, LO e2) {

  Few<Real, 2> xi1;
  if (old_vert == v0) {
    OMEGA_H_CHECK(old_edge == e1);
    xi1[0] = 0.25;
    xi1[1] = 0.25;
  }
  else if (old_vert == v1) {
    OMEGA_H_CHECK(old_edge == e2);
    xi1[0] = 0.5;
    xi1[1] = 0.25;
  }
  else if (old_vert == v2) {
    OMEGA_H_CHECK(old_edge == e0);
    xi1[0] = 0.25;
    xi1[1] = 0.5;
  }
  else {
  }

  return xi1;
}

OMEGA_H_INLINE Few<Real,4> cubic_noKeyEdge_xi_values(
    LO old_vert, LO v0, LO v1, LO v2, LO old_edge, LO e0, LO e1, LO e2) {

  Few<Real, 4> xi1_xi2;
  if (old_vert == v0) {
    OMEGA_H_CHECK(old_edge == e1);
    xi1_xi2[0] = 0.13740215;
    xi1_xi2[1] = 0.13740215;
    xi1_xi2[2] = 0.362597849;
    xi1_xi2[3] = 0.362597849;
  }
  else if (old_vert == v1) {
    OMEGA_H_CHECK(old_edge == e2);
    xi1_xi2[0] = 0.7251957;
    xi1_xi2[1] = 0.13740215;
    xi1_xi2[2] = 0.2748043;
    xi1_xi2[3] = 0.362597849;
  }
  else if (old_vert == v2) {
    OMEGA_H_CHECK(old_edge == e0);
    xi1_xi2[0] = 0.13740215;
    xi1_xi2[1] = 0.7251957;
    xi1_xi2[2] = 0.362597849;
    xi1_xi2[3] = 0.2748043;
  }
  else {
  }

  return xi1_xi2;
}

OMEGA_H_INLINE Few<Real,4> cubic_faceSplittingEdge_xi_values
  (LO const old_noKey_vert, LO const old_v0, LO const old_v1, LO const old_v2,
   LO const old_key_edge, LO const old_e0, LO const old_e1, LO const old_e2,
   LO const new_v0_f0, LO const new_v1_f0, LO const new_v2_f0,
   LO const new_v0_f1, LO const new_v1_f1, LO const new_v2_f1,
   LO const oldf_v0_new, LO const oldf_v1_new, LO const oldf_v2_new) {

  Few<Real, 4> xi1_xi2;
  I8 should_swap = -1;
  if (old_noKey_vert == old_v0) {
    OMEGA_H_CHECK(old_key_edge == old_e1);
    xi1_xi2[0] = 1.0/2.0;
    xi1_xi2[1] = 1.0/6.0;
    xi1_xi2[2] = 1.0/6.0;
    xi1_xi2[3] = 1.0/2.0;

    if ((oldf_v1_new == new_v0_f0) || (oldf_v1_new == new_v1_f0) || 
        (oldf_v1_new == new_v2_f0)) {
      should_swap = -1;
    }
    else {
      OMEGA_H_CHECK ((oldf_v1_new == new_v0_f1) || (oldf_v1_new == new_v1_f1) ||
                     (oldf_v1_new == new_v2_f1));
      should_swap = 1;
    }
  }
  else if (old_noKey_vert == old_v1) {
    OMEGA_H_CHECK(old_key_edge == old_e2);
    xi1_xi2[0] = 1.0/3.0;
    xi1_xi2[1] = 1.0/2.0;
    xi1_xi2[2] = 1.0/3.0;
    xi1_xi2[3] = 1.0/6.0;
    if ((oldf_v2_new == new_v0_f0) || (oldf_v2_new == new_v1_f0) ||
        (oldf_v2_new == new_v2_f0)) {
      should_swap = -1;
    }
    else {
      OMEGA_H_CHECK ((oldf_v2_new == new_v0_f1) || (oldf_v2_new == new_v1_f1) ||
                     (oldf_v2_new == new_v2_f1));
      should_swap = 1;
    }
  }
  else if (old_noKey_vert == old_v2) {
    OMEGA_H_CHECK(old_key_edge == old_e0);
    xi1_xi2[0] = 1.0/6.0;
    xi1_xi2[1] = 1.0/3.0;
    xi1_xi2[2] = 1.0/2.0;
    xi1_xi2[3] = 1.0/3.0;
    if ((oldf_v0_new == new_v0_f0) || (oldf_v0_new == new_v1_f0) ||
        (oldf_v0_new == new_v2_f0)) {
      should_swap = -1;
    }
    else {
      OMEGA_H_CHECK ((oldf_v0_new == new_v0_f1) || (oldf_v0_new == new_v1_f1) ||
                     (oldf_v0_new == new_v2_f1));
      should_swap = 1;
    }
  }
  else {
  }
  if (should_swap == 1) {
    swap2(xi1_xi2[0], xi1_xi2[2]);
    swap2(xi1_xi2[1], xi1_xi2[3]);
  }

  return xi1_xi2;
}

OMEGA_H_INLINE Vector<3> cubic_region_xi_values
  (LO const old_key_edge, LO const old_e0, LO const old_e1, LO const old_e2,
   LO const old_e3, LO const old_e4, LO const old_e5) {
  Vector<3> p11;
  if (old_key_edge == old_e0) {
    p11[0] = 1.0/6.0;
    p11[1] = 1.0/3.0;
    p11[2] = 1.0/3.0;
  }
  else if (old_key_edge == old_e1) {
    p11[0] = 1.0/6.0;
    p11[1] = 1.0/6.0;
    p11[2] = 1.0/3.0;
  }
  else if (old_key_edge == old_e2) {
    p11[0] = 1.0/3.0;
    p11[1] = 1.0/6.0;
    p11[2] = 1.0/3.0;
  }
  else if (old_key_edge == old_e3) {
    p11[0] = 1.0/3.0;
    p11[1] = 1.0/3.0;
    p11[2] = 1.0/6.0;
  }
  else if (old_key_edge == old_e4) {
    p11[0] = 1.0/6.0;
    p11[1] = 1.0/3.0;
    p11[2] = 1.0/6.0;
  }
  else if (old_key_edge == old_e5) {
    p11[0] = 1.0/3.0;
    p11[1] = 1.0/6.0;
    p11[2] = 1.0/6.0;
  }
  else {
  }

  return p11;
}

OMEGA_H_INLINE LO edge_is_flip(LO const e0v0, LO const e0v1, LO const v0,
                               LO const v1) {
  LO is_flip = -1;
  if ((e0v0 == v1) && (e0v1 == v0)) {
    is_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
  }
  return is_flip;
}

OMEGA_H_INLINE Int n_internal_ctrlPts(Int edim, Int max_order) {
  OMEGA_H_CHECK(max_order > 0);
  OMEGA_H_CHECK(edim >= 0);
  if (edim == 0) {
    return 1;
  }
  else if (edim == 1) {
    return max_order-1;
  }
  else if (edim == 2) {
    return ((max_order-1)*(max_order-2))/2;
  }
  else if (edim == 3) {
    return ((max_order-1)*(max_order-2)*(max_order-3))/6;
  }
  else {
    return -1;
  }
  return -1;
}

OMEGA_H_INLINE Int n_total_ctrlPts(Int edim, Int max_order) {
  OMEGA_H_CHECK(max_order > 0);
  OMEGA_H_CHECK(edim >= 0);
  if (edim == 0) {
    return 1;
  }
  else if (edim == 1) {
    return max_order+1;
  }
  else if (edim == 2) {
    return ((max_order+1)*(max_order+2))/2;
  }
  else if (edim == 3) {
    return ((max_order+1)*(max_order+2)*(max_order+3))/6;
  }
  else {
    return -1;
  }
  return -1;
}

inline Vector<3> rgn_parametricToParent_3dp2_h(
    LO const order, LO const old_rgn, HostRead<LO> old_rv2v,
    HostRead<Real> old_vertCtrlPts, HostRead<Real> old_edgeCtrlPts, 
    Vector<3> nodePt, HostRead<LO> old_re2e) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  auto old_rgn_v0 = old_rv2v[old_rgn*4 + 0];
  auto old_rgn_v1 = old_rv2v[old_rgn*4 + 1];
  auto old_rgn_v2 = old_rv2v[old_rgn*4 + 2];
  auto old_rgn_v3 = old_rv2v[old_rgn*4 + 3];
  Vector<3> c000, c200, c020, c002;
  for (LO d=0; d<3; ++d) {
    c000[d] = old_vertCtrlPts[old_rgn_v0*dim + d];
    c200[d] = old_vertCtrlPts[old_rgn_v1*dim + d];
    c020[d] = old_vertCtrlPts[old_rgn_v2*dim + d];
    c002[d] = old_vertCtrlPts[old_rgn_v3*dim + d];
  }

  auto old_rgn_e0 = old_re2e[old_rgn*6 + 0];
  auto old_rgn_e1 = old_re2e[old_rgn*6 + 1];
  auto old_rgn_e2 = old_re2e[old_rgn*6 + 2];
  auto old_rgn_e3 = old_re2e[old_rgn*6 + 3];
  auto old_rgn_e4 = old_re2e[old_rgn*6 + 4];
  auto old_rgn_e5 = old_re2e[old_rgn*6 + 5];

  auto pts_per_edge = n_edge_pts;
  Real cx100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 0];
  Real cy100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 1];
  Real cz100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 2];
  Real cx110 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 0];
  Real cy110 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 1];
  Real cz110 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 2];
  Real cx010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 0];
  Real cy010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 1];
  Real cz010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 2];
  Real cx001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 0];
  Real cy001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 1];
  Real cz001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 2];
  Real cx101 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 0];
  Real cy101 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 1];
  Real cz101 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 2];
  Real cx011 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 0];
  Real cy011 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 1];
  Real cz011 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 2];

  auto c100 = vector_3(cx100, cy100, cz100);
  auto c110 = vector_3(cx110, cy110, cz110);
  auto c010 = vector_3(cx010, cy010, cz010);
  auto c001 = vector_3(cx001, cy001, cz001);
  auto c101 = vector_3(cx101, cy101, cz101);
  auto c011 = vector_3(cx011, cy011, cz011);

  Vector<3> p11;
  for (LO j = 0; j < dim; ++j) {
    p11[j] = c000[j]*Bijk(order,0,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c200[j]*Bijk(order,2,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c020[j]*Bijk(order,0,2,0,nodePt[0],nodePt[1],nodePt[2]) +
      c002[j]*Bijk(order,0,0,2,nodePt[0],nodePt[1],nodePt[2]) +
      c100[j]*Bijk(order,1,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c110[j]*Bijk(order,1,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c010[j]*Bijk(order,0,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c001[j]*Bijk(order,0,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c101[j]*Bijk(order,1,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c011[j]*Bijk(order,0,1,1,nodePt[0],nodePt[1],nodePt[2]);
  }
  return p11;
}

inline Vector<3> rgn_parametricToParent_3dp2(
    LO const order, LO const old_rgn, Read<LO> old_rv2v,
    Read<Real> old_vertCtrlPts, Read<Real> old_edgeCtrlPts, 
    Vector<3> nodePt, HostRead<LO> old_re2e) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  auto old_rgn_v0 = old_rv2v[old_rgn*4 + 0];
  auto old_rgn_v1 = old_rv2v[old_rgn*4 + 1];
  auto old_rgn_v2 = old_rv2v[old_rgn*4 + 2];
  auto old_rgn_v3 = old_rv2v[old_rgn*4 + 3];
  Vector<3> c000, c200, c020, c002;
  for (LO d=0; d<3; ++d) {
    c000[d] = old_vertCtrlPts[old_rgn_v0*dim + d];
    c200[d] = old_vertCtrlPts[old_rgn_v1*dim + d];
    c020[d] = old_vertCtrlPts[old_rgn_v2*dim + d];
    c002[d] = old_vertCtrlPts[old_rgn_v3*dim + d];
  }

  auto old_rgn_e0 = old_re2e[old_rgn*6 + 0];
  auto old_rgn_e1 = old_re2e[old_rgn*6 + 1];
  auto old_rgn_e2 = old_re2e[old_rgn*6 + 2];
  auto old_rgn_e3 = old_re2e[old_rgn*6 + 3];
  auto old_rgn_e4 = old_re2e[old_rgn*6 + 4];
  auto old_rgn_e5 = old_re2e[old_rgn*6 + 5];

  auto pts_per_edge = n_edge_pts;
  Real cx100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 0];
  Real cy100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 1];
  Real cz100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 2];
  Real cx110 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 0];
  Real cy110 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 1];
  Real cz110 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 2];
  Real cx010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 0];
  Real cy010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 1];
  Real cz010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 2];
  Real cx001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 0];
  Real cy001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 1];
  Real cz001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 2];
  Real cx101 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 0];
  Real cy101 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 1];
  Real cz101 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 2];
  Real cx011 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 0];
  Real cy011 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 1];
  Real cz011 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 2];

  auto c100 = vector_3(cx100, cy100, cz100);
  auto c110 = vector_3(cx110, cy110, cz110);
  auto c010 = vector_3(cx010, cy010, cz010);
  auto c001 = vector_3(cx001, cy001, cz001);
  auto c101 = vector_3(cx101, cy101, cz101);
  auto c011 = vector_3(cx011, cy011, cz011);

  Vector<3> p11;
  for (LO j = 0; j < dim; ++j) {
    p11[j] = c000[j]*Bijk(order,0,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c200[j]*Bijk(order,2,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c020[j]*Bijk(order,0,2,0,nodePt[0],nodePt[1],nodePt[2]) +
      c002[j]*Bijk(order,0,0,2,nodePt[0],nodePt[1],nodePt[2]) +
      c100[j]*Bijk(order,1,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c110[j]*Bijk(order,1,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c010[j]*Bijk(order,0,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c001[j]*Bijk(order,0,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c101[j]*Bijk(order,1,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c011[j]*Bijk(order,0,1,1,nodePt[0],nodePt[1],nodePt[2]);
  }
  return p11;
}

inline Vector<3> rgn_parametricToParent_3d_h(
    LO const order, LO const old_rgn, HostRead<LO> old_ev2v, HostRead<LO> old_rv2v,
    HostRead<Real> old_vertCtrlPts, HostRead<Real> old_edgeCtrlPts, 
    HostRead<Real> old_faceCtrlPts, Vector<3> nodePt, HostRead<LO> old_re2e, 
    HostRead<LO> old_rf2f) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  auto old_rgn_v0 = old_rv2v[old_rgn*4 + 0];
  auto old_rgn_v1 = old_rv2v[old_rgn*4 + 1];
  auto old_rgn_v2 = old_rv2v[old_rgn*4 + 2];
  auto old_rgn_v3 = old_rv2v[old_rgn*4 + 3];
  Vector<3> c000, c300, c030, c003;
  for (LO d=0; d<3; ++d) {
    c000[d] = old_vertCtrlPts[old_rgn_v0*dim + d];
    c300[d] = old_vertCtrlPts[old_rgn_v1*dim + d];
    c030[d] = old_vertCtrlPts[old_rgn_v2*dim + d];
    c003[d] = old_vertCtrlPts[old_rgn_v3*dim + d];
  }

  auto old_rgn_e0 = old_re2e[old_rgn*6 + 0];
  auto old_rgn_e1 = old_re2e[old_rgn*6 + 1];
  auto old_rgn_e2 = old_re2e[old_rgn*6 + 2];
  auto old_rgn_e3 = old_re2e[old_rgn*6 + 3];
  auto old_rgn_e4 = old_re2e[old_rgn*6 + 4];
  auto old_rgn_e5 = old_re2e[old_rgn*6 + 5];

  auto old_rgn_f0 = old_rf2f[old_rgn*4 + 0];
  auto old_rgn_f1 = old_rf2f[old_rgn*4 + 1];
  auto old_rgn_f2 = old_rf2f[old_rgn*4 + 2];
  auto old_rgn_f3 = old_rf2f[old_rgn*4 + 3];
  Vector<3> c110, c101, c111, c011;
  for (LO d=0; d<3; ++d) {
    c110[d] = old_faceCtrlPts[old_rgn_f0*dim + d];
    c101[d] = old_faceCtrlPts[old_rgn_f1*dim + d];
    c111[d] = old_faceCtrlPts[old_rgn_f2*dim + d];
    c011[d] = old_faceCtrlPts[old_rgn_f3*dim + d];
  }

  LO e0_flip = -1;
  LO e1_flip = -1;
  LO e2_flip = -1;
  LO e3_flip = -1;
  LO e4_flip = -1;
  LO e5_flip = -1;
  {
    auto e0v0 = old_ev2v[old_rgn_e0*2 + 0];
    auto e0v1 = old_ev2v[old_rgn_e0*2 + 1];
    e0_flip = edge_is_flip(e0v0, e0v1, old_rgn_v0, old_rgn_v1);
    auto e1v0 = old_ev2v[old_rgn_e1*2 + 0];
    auto e1v1 = old_ev2v[old_rgn_e1*2 + 1];
    e1_flip = edge_is_flip(e1v0, e1v1, old_rgn_v1, old_rgn_v2);
    auto e2v0 = old_ev2v[old_rgn_e2*2 + 0];
    auto e2v1 = old_ev2v[old_rgn_e2*2 + 1];
    e2_flip = edge_is_flip(e2v0, e2v1, old_rgn_v2, old_rgn_v0);
    auto e3v0 = old_ev2v[old_rgn_e3*2 + 0];
    auto e3v1 = old_ev2v[old_rgn_e3*2 + 1];
    e3_flip = edge_is_flip(e3v0, e3v1, old_rgn_v0, old_rgn_v3);
    auto e4v0 = old_ev2v[old_rgn_e4*2 + 0];
    auto e4v1 = old_ev2v[old_rgn_e4*2 + 1];
    e4_flip = edge_is_flip(e4v0, e4v1, old_rgn_v1, old_rgn_v3);
    auto e5v0 = old_ev2v[old_rgn_e5*2 + 0];
    auto e5v1 = old_ev2v[old_rgn_e5*2 + 1];
    e5_flip = edge_is_flip(e5v0, e5v1, old_rgn_v2, old_rgn_v3);
  }

  auto pts_per_edge = n_edge_pts;
  Real cx100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 0];
  Real cy100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 1];
  Real cz100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 2];
  Real cx200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e0_flip > 0) {
    swap2(cx100, cx200);
    swap2(cy100, cy200);
    swap2(cz100, cz200);
  }
  Real cx210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 0];
  Real cy210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 1];
  Real cz210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 2];
  Real cx120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e1_flip > 0) {
    swap2(cx210, cx120);
    swap2(cy210, cy120);
    swap2(cz210, cz120);
  }
  Real cx020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 0];
  Real cy020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 1];
  Real cz020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 2];
  Real cx010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e2_flip > 0) {
    swap2(cx020, cx010);
    swap2(cy020, cy010);
    swap2(cz020, cz010);
  }
  Real cx001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 0];
  Real cy001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 1];
  Real cz001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 2];
  Real cx002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e3_flip > 0) {
    swap2(cx002, cx001);
    swap2(cy002, cy001);
    swap2(cz002, cz001);
  }
  Real cx201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 0];
  Real cy201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 1];
  Real cz201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 2];
  Real cx102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e4_flip > 0) {
    swap2(cx102, cx201);
    swap2(cy102, cy201);
    swap2(cz102, cz201);
  }
  Real cx021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 0];
  Real cy021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 1];
  Real cz021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 2];
  Real cx012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e5_flip > 0) {
    swap2(cx012, cx021);
    swap2(cy012, cy021);
    swap2(cz012, cz021);
  }
  auto c100 = vector_3(cx100, cy100, cz100);
  auto c200 = vector_3(cx200, cy200, cz200);
  auto c210 = vector_3(cx210, cy210, cz210);
  auto c120 = vector_3(cx120, cy120, cz120);
  auto c020 = vector_3(cx020, cy020, cz020);
  auto c010 = vector_3(cx010, cy010, cz010);
  auto c001 = vector_3(cx001, cy001, cz001);
  auto c002 = vector_3(cx002, cy002, cz002);
  auto c102 = vector_3(cx102, cy102, cz102);
  auto c201 = vector_3(cx201, cy201, cz201);
  auto c012 = vector_3(cx012, cy012, cz012);
  auto c021 = vector_3(cx021, cy021, cz021);

  Vector<3> p11;
  for (LO j = 0; j < dim; ++j) {
    p11[j] = c000[j]*Bijk(order,0,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c300[j]*Bijk(order,3,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c030[j]*Bijk(order,0,3,0,nodePt[0],nodePt[1],nodePt[2]) +
      c003[j]*Bijk(order,0,0,3,nodePt[0],nodePt[1],nodePt[2]) +
      c100[j]*Bijk(order,1,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c200[j]*Bijk(order,2,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c210[j]*Bijk(order,2,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c120[j]*Bijk(order,1,2,0,nodePt[0],nodePt[1],nodePt[2]) +
      c020[j]*Bijk(order,0,2,0,nodePt[0],nodePt[1],nodePt[2]) +
      c010[j]*Bijk(order,0,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c001[j]*Bijk(order,0,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c002[j]*Bijk(order,0,0,2,nodePt[0],nodePt[1],nodePt[2]) +
      c201[j]*Bijk(order,2,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c102[j]*Bijk(order,1,0,2,nodePt[0],nodePt[1],nodePt[2]) +
      c012[j]*Bijk(order,0,1,2,nodePt[0],nodePt[1],nodePt[2]) +
      c021[j]*Bijk(order,0,2,1,nodePt[0],nodePt[1],nodePt[2]) +
      c110[j]*Bijk(order,1,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c101[j]*Bijk(order,1,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c111[j]*Bijk(order,1,1,1,nodePt[0],nodePt[1],nodePt[2]) +
      c011[j]*Bijk(order,0,1,1,nodePt[0],nodePt[1],nodePt[2]);
  }
  return p11;
}

OMEGA_H_DEVICE Vector<3> rgn_parametricToParent_3d(
    LO const order, LO const old_rgn, LOs old_ev2v, LOs old_rv2v,
    Reals old_vertCtrlPts, Reals old_edgeCtrlPts, Reals old_faceCtrlPts,
    Vector<3> nodePt, LOs old_re2e, LOs old_rf2f) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  auto old_rgn_v0 = old_rv2v[old_rgn*4 + 0];
  auto old_rgn_v1 = old_rv2v[old_rgn*4 + 1];
  auto old_rgn_v2 = old_rv2v[old_rgn*4 + 2];
  auto old_rgn_v3 = old_rv2v[old_rgn*4 + 3];
  auto c000 = get_vector<3>(old_vertCtrlPts, old_rgn_v0);
  auto c300 = get_vector<3>(old_vertCtrlPts, old_rgn_v1);
  auto c030 = get_vector<3>(old_vertCtrlPts, old_rgn_v2);
  auto c003 = get_vector<3>(old_vertCtrlPts, old_rgn_v3);

  auto old_rgn_e0 = old_re2e[old_rgn*6 + 0];
  auto old_rgn_e1 = old_re2e[old_rgn*6 + 1];
  auto old_rgn_e2 = old_re2e[old_rgn*6 + 2];
  auto old_rgn_e3 = old_re2e[old_rgn*6 + 3];
  auto old_rgn_e4 = old_re2e[old_rgn*6 + 4];
  auto old_rgn_e5 = old_re2e[old_rgn*6 + 5];

  auto old_rgn_f0 = old_rf2f[old_rgn*4 + 0];
  auto old_rgn_f1 = old_rf2f[old_rgn*4 + 1];
  auto old_rgn_f2 = old_rf2f[old_rgn*4 + 2];
  auto old_rgn_f3 = old_rf2f[old_rgn*4 + 3];
  auto c110 = get_vector<3>(old_faceCtrlPts, old_rgn_f0);
  auto c101 = get_vector<3>(old_faceCtrlPts, old_rgn_f1);
  auto c111 = get_vector<3>(old_faceCtrlPts, old_rgn_f2);
  auto c011 = get_vector<3>(old_faceCtrlPts, old_rgn_f3);

  LO e0_flip = -1;
  LO e1_flip = -1;
  LO e2_flip = -1;
  LO e3_flip = -1;
  LO e4_flip = -1;
  LO e5_flip = -1;
  {
    auto e0v0 = old_ev2v[old_rgn_e0*2 + 0];
    auto e0v1 = old_ev2v[old_rgn_e0*2 + 1];
    e0_flip = edge_is_flip(e0v0, e0v1, old_rgn_v0, old_rgn_v1);
    auto e1v0 = old_ev2v[old_rgn_e1*2 + 0];
    auto e1v1 = old_ev2v[old_rgn_e1*2 + 1];
    e1_flip = edge_is_flip(e1v0, e1v1, old_rgn_v1, old_rgn_v2);
    auto e2v0 = old_ev2v[old_rgn_e2*2 + 0];
    auto e2v1 = old_ev2v[old_rgn_e2*2 + 1];
    e2_flip = edge_is_flip(e2v0, e2v1, old_rgn_v2, old_rgn_v0);
    auto e3v0 = old_ev2v[old_rgn_e3*2 + 0];
    auto e3v1 = old_ev2v[old_rgn_e3*2 + 1];
    e3_flip = edge_is_flip(e3v0, e3v1, old_rgn_v0, old_rgn_v3);
    auto e4v0 = old_ev2v[old_rgn_e4*2 + 0];
    auto e4v1 = old_ev2v[old_rgn_e4*2 + 1];
    e4_flip = edge_is_flip(e4v0, e4v1, old_rgn_v1, old_rgn_v3);
    auto e5v0 = old_ev2v[old_rgn_e5*2 + 0];
    auto e5v1 = old_ev2v[old_rgn_e5*2 + 1];
    e5_flip = edge_is_flip(e5v0, e5v1, old_rgn_v2, old_rgn_v3);
  }

  auto pts_per_edge = n_edge_pts;
  Real cx100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 0];
  Real cy100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 1];
  Real cz100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 2];
  Real cx200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e0_flip > 0) {
    swap2(cx100, cx200);
    swap2(cy100, cy200);
    swap2(cz100, cz200);
  }
  Real cx210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 0];
  Real cy210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 1];
  Real cz210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 2];
  Real cx120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e1_flip > 0) {
    swap2(cx210, cx120);
    swap2(cy210, cy120);
    swap2(cz210, cz120);
  }
  Real cx020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 0];
  Real cy020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 1];
  Real cz020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 2];
  Real cx010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e2_flip > 0) {
    swap2(cx020, cx010);
    swap2(cy020, cy010);
    swap2(cz020, cz010);
  }
  Real cx001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 0];
  Real cy001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 1];
  Real cz001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 2];
  Real cx002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e3_flip > 0) {
    swap2(cx002, cx001);
    swap2(cy002, cy001);
    swap2(cz002, cz001);
  }
  Real cx201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 0];
  Real cy201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 1];
  Real cz201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 2];
  Real cx102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e4_flip > 0) {
    swap2(cx102, cx201);
    swap2(cy102, cy201);
    swap2(cz102, cz201);
  }
  Real cx021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 0];
  Real cy021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 1];
  Real cz021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 2];
  Real cx012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e5_flip > 0) {
    swap2(cx012, cx021);
    swap2(cy012, cy021);
    swap2(cz012, cz021);
  }
  auto c100 = vector_3(cx100, cy100, cz100);
  auto c200 = vector_3(cx200, cy200, cz200);
  auto c210 = vector_3(cx210, cy210, cz210);
  auto c120 = vector_3(cx120, cy120, cz120);
  auto c020 = vector_3(cx020, cy020, cz020);
  auto c010 = vector_3(cx010, cy010, cz010);
  auto c001 = vector_3(cx001, cy001, cz001);
  auto c002 = vector_3(cx002, cy002, cz002);
  auto c102 = vector_3(cx102, cy102, cz102);
  auto c201 = vector_3(cx201, cy201, cz201);
  auto c012 = vector_3(cx012, cy012, cz012);
  auto c021 = vector_3(cx021, cy021, cz021);

  Vector<3> p11;
  for (LO j = 0; j < dim; ++j) {
    p11[j] = c000[j]*Bijk(order,0,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c300[j]*Bijk(order,3,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c030[j]*Bijk(order,0,3,0,nodePt[0],nodePt[1],nodePt[2]) +
      c003[j]*Bijk(order,0,0,3,nodePt[0],nodePt[1],nodePt[2]) +
      c100[j]*Bijk(order,1,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c200[j]*Bijk(order,2,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c210[j]*Bijk(order,2,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c120[j]*Bijk(order,1,2,0,nodePt[0],nodePt[1],nodePt[2]) +
      c020[j]*Bijk(order,0,2,0,nodePt[0],nodePt[1],nodePt[2]) +
      c010[j]*Bijk(order,0,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c001[j]*Bijk(order,0,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c002[j]*Bijk(order,0,0,2,nodePt[0],nodePt[1],nodePt[2]) +
      c201[j]*Bijk(order,2,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c102[j]*Bijk(order,1,0,2,nodePt[0],nodePt[1],nodePt[2]) +
      c012[j]*Bijk(order,0,1,2,nodePt[0],nodePt[1],nodePt[2]) +
      c021[j]*Bijk(order,0,2,1,nodePt[0],nodePt[1],nodePt[2]) +
      c110[j]*Bijk(order,1,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c101[j]*Bijk(order,1,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c111[j]*Bijk(order,1,1,1,nodePt[0],nodePt[1],nodePt[2]) +
      c011[j]*Bijk(order,0,1,1,nodePt[0],nodePt[1],nodePt[2]);
  }
  return p11;
}

OMEGA_H_DEVICE Vector<2> face_parametricToParent_2d(
    LO const order, LO const old_face, LOs old_ev2v, LOs old_fe2e,
    Reals old_vertCtrlPts, Reals old_edgeCtrlPts, Reals old_faceCtrlPts,
    Real nodePts_0, Real nodePts_1, LOs old_fv2v) {
  LO const dim = 2;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  LO const v0_old_face = old_fv2v[old_face*3 + 0];
  LO const v1_old_face = old_fv2v[old_face*3 + 1];
  LO const v2_old_face = old_fv2v[old_face*3 + 2];
  LO const old_face_e0 = old_fe2e[old_face*3 + 0];
  LO const old_face_e1 = old_fe2e[old_face*3 + 1];
  LO const old_face_e2 = old_fe2e[old_face*3 + 2];
  I8 e0_flip = -1;
  I8 e1_flip = -1;
  I8 e2_flip = -1;
  LO v1 = v1_old_face;
  LO v2 = v2_old_face;
  auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
  auto e0v1 = old_ev2v[old_face_e0*2 + 1];
  auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
  auto e1v1 = old_ev2v[old_face_e1*2 + 1];
  auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
  auto e2v1 = old_ev2v[old_face_e2*2 + 1];
  if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
    e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
  }
  if ((e1v0_old_face == v2) && (e1v1 == v1)) {
    e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
  }
  if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
    e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
  }

  Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
  Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
  Real cx30 = old_vertCtrlPts[v1*dim + 0];
  Real cy30 = old_vertCtrlPts[v1*dim + 1];
  Real cx03 = old_vertCtrlPts[v2*dim + 0];
  Real cy03 = old_vertCtrlPts[v2*dim + 1];

  auto pts_per_edge = n_edge_pts;
  Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
  Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
  Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e0_flip > 0) {
    swap2(cx10, cx20);
    swap2(cy10, cy20);
  }

  Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
  Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
  Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e1_flip > 0) {
    swap2(cx12, cx21);
    swap2(cy12, cy21);
  }

  Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
  Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
  Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e2_flip > 0) {
    swap2(cx02, cx01);
    swap2(cy02, cy01);
  }

  Real cx11 = old_faceCtrlPts[old_face*dim + 0];
  Real cy11 = old_faceCtrlPts[old_face*dim + 1];
  auto c00 = vector_2(cx00, cy00);
  auto c10 = vector_2(cx10, cy10);
  auto c20 = vector_2(cx20, cy20);
  auto c30 = vector_2(cx30, cy30);
  auto c21 = vector_2(cx21, cy21);
  auto c12 = vector_2(cx12, cy12);
  auto c03 = vector_2(cx03, cy03);
  auto c02 = vector_2(cx02, cy02);
  auto c01 = vector_2(cx01, cy01);
  auto c11 = vector_2(cx11, cy11);

  Vector<2> p11_w;
  for (LO k = 0; k < dim; ++k) {
    p11_w[k] = c00[k]*Bij(order, 0, 0, nodePts_0, nodePts_1) +
      c10[k]*Bij(order, 1, 0, nodePts_0, nodePts_1) +
      c20[k]*Bij(order, 2, 0, nodePts_0, nodePts_1) +
      c30[k]*Bij(order, 3, 0, nodePts_0, nodePts_1) +
      c21[k]*Bij(order, 2, 1, nodePts_0, nodePts_1) +
      c12[k]*Bij(order, 1, 2, nodePts_0, nodePts_1) +
      c03[k]*Bij(order, 0, 3, nodePts_0, nodePts_1) +
      c02[k]*Bij(order, 0, 2, nodePts_0, nodePts_1) +
      c01[k]*Bij(order, 0, 1, nodePts_0, nodePts_1) +
      c11[k]*Bij(order, 1, 1, nodePts_0, nodePts_1);
  }
  return p11_w;
}

inline Vector<2> face_parametricToParent_2d_h(
    LO const order, LO const old_face, HostRead<LO> old_ev2v,
    HostRead<LO> old_fe2e, HostRead<Real> old_vertCtrlPts, 
    HostRead<Real> old_edgeCtrlPts, HostRead<Real> old_faceCtrlPts,
    Real nodePts_0, Real nodePts_1, HostRead<LO> old_fv2v) {
  LO const dim = 2;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  LO const v0_old_face = old_fv2v[old_face*3 + 0];
  LO const v1_old_face = old_fv2v[old_face*3 + 1];
  LO const v2_old_face = old_fv2v[old_face*3 + 2];
  LO const old_face_e0 = old_fe2e[old_face*3 + 0];
  LO const old_face_e1 = old_fe2e[old_face*3 + 1];
  LO const old_face_e2 = old_fe2e[old_face*3 + 2];
  I8 e0_flip = -1;
  I8 e1_flip = -1;
  I8 e2_flip = -1;
  LO v1 = v1_old_face;
  LO v2 = v2_old_face;
  auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
  auto e0v1 = old_ev2v[old_face_e0*2 + 1];
  auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
  auto e1v1 = old_ev2v[old_face_e1*2 + 1];
  auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
  auto e2v1 = old_ev2v[old_face_e2*2 + 1];
  if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
    e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
  }
  if ((e1v0_old_face == v2) && (e1v1 == v1)) {
    e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
  }
  if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
    e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
  }

  Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
  Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
  Real cx30 = old_vertCtrlPts[v1*dim + 0];
  Real cy30 = old_vertCtrlPts[v1*dim + 1];
  Real cx03 = old_vertCtrlPts[v2*dim + 0];
  Real cy03 = old_vertCtrlPts[v2*dim + 1];

  auto pts_per_edge = n_edge_pts;
  Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
  Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
  Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e0_flip > 0) {
    swap2(cx10, cx20);
    swap2(cy10, cy20);
  }

  Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
  Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
  Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e1_flip > 0) {
    swap2(cx12, cx21);
    swap2(cy12, cy21);
  }

  Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
  Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
  Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e2_flip > 0) {
    swap2(cx02, cx01);
    swap2(cy02, cy01);
  }

  Real cx11 = old_faceCtrlPts[old_face*dim + 0];
  Real cy11 = old_faceCtrlPts[old_face*dim + 1];
  auto c00 = vector_2(cx00, cy00);
  auto c10 = vector_2(cx10, cy10);
  auto c20 = vector_2(cx20, cy20);
  auto c30 = vector_2(cx30, cy30);
  auto c21 = vector_2(cx21, cy21);
  auto c12 = vector_2(cx12, cy12);
  auto c03 = vector_2(cx03, cy03);
  auto c02 = vector_2(cx02, cy02);
  auto c01 = vector_2(cx01, cy01);
  auto c11 = vector_2(cx11, cy11);

  Vector<2> p11_w;
  for (LO k = 0; k < dim; ++k) {
    p11_w[k] = c00[k]*Bij(order, 0, 0, nodePts_0, nodePts_1) +
      c10[k]*Bij(order, 1, 0, nodePts_0, nodePts_1) +
      c20[k]*Bij(order, 2, 0, nodePts_0, nodePts_1) +
      c30[k]*Bij(order, 3, 0, nodePts_0, nodePts_1) +
      c21[k]*Bij(order, 2, 1, nodePts_0, nodePts_1) +
      c12[k]*Bij(order, 1, 2, nodePts_0, nodePts_1) +
      c03[k]*Bij(order, 0, 3, nodePts_0, nodePts_1) +
      c02[k]*Bij(order, 0, 2, nodePts_0, nodePts_1) +
      c01[k]*Bij(order, 0, 1, nodePts_0, nodePts_1) +
      c11[k]*Bij(order, 1, 1, nodePts_0, nodePts_1);
  }
  return p11_w;
}

OMEGA_H_DEVICE Few<LO, 3> calc_edge_flips(LO const v0, LO const v1,
    LO const v2, LO const e0_v0, LO const e0_v1, LO const e1_v0, 
    LO const e1_v1, LO const e2_v0, LO const e2_v1) {
  Few<LO, 3> edge_flips;
  LO e0_flip = -1;
  LO e1_flip = -1;
  LO e2_flip = -1;
  if ((e0_v0 == v1) && (e0_v1 == v0)) {
    e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0_v0 == v0) && (e0_v1 == v1));
  }
  if ((e1_v0 == v2) && (e1_v1 == v1)) {
    e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e1_v0 == v1) && (e1_v1 == v2));
  }
  if ((e2_v0 == v0) && (e2_v1 == v2)) {
    e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e2_v0 == v2) && (e2_v1 == v0));
  }
  edge_flips[0] = e0_flip;
  edge_flips[1] = e1_flip;
  edge_flips[2] = e2_flip;

  return edge_flips;
}

OMEGA_H_DEVICE Vector<3> face_parametricToParent_3d_p2(
    LO const order, LO const old_face, LOs old_fe2e,
    Reals old_vertCtrlPts, Reals old_edgeCtrlPts,
    Real nodePts_0, Real nodePts_1, LOs old_fv2v) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  LO const v0_old_face = old_fv2v[old_face*3 + 0];
  LO const v1_old_face = old_fv2v[old_face*3 + 1];
  LO const v2_old_face = old_fv2v[old_face*3 + 2];
  LO const old_face_e0 = old_fe2e[old_face*3 + 0];
  LO const old_face_e1 = old_fe2e[old_face*3 + 1];
  LO const old_face_e2 = old_fe2e[old_face*3 + 2];
  LO v1 = v1_old_face;
  LO v2 = v2_old_face;

  Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
  Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
  Real cz00 = old_vertCtrlPts[v0_old_face*dim + 2];
  Real cx20 = old_vertCtrlPts[v1*dim + 0];
  Real cy20 = old_vertCtrlPts[v1*dim + 1];
  Real cz20 = old_vertCtrlPts[v1*dim + 2];
  Real cx02 = old_vertCtrlPts[v2*dim + 0];
  Real cy02 = old_vertCtrlPts[v2*dim + 1];
  Real cz02 = old_vertCtrlPts[v2*dim + 2];

  auto pts_per_edge = n_edge_pts;
  Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
  Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
  Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
  Real cx11 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
  Real cy11 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
  Real cz11 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
  Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
  Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
  Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];

  auto c00 = vector_3(cx00, cy00, cz00);
  auto c10 = vector_3(cx10, cy10, cz10);
  auto c20 = vector_3(cx20, cy20, cz20);
  auto c11 = vector_3(cx11, cy11, cz11);
  auto c02 = vector_3(cx02, cy02, cz02);
  auto c01 = vector_3(cx01, cy01, cz01);

  Vector<3> p11_w;
  for (LO k = 0; k < 3; ++k) {
    p11_w[k] = c00[k]*Bij(order, 0, 0, nodePts_0, nodePts_1) +
      c10[k]*Bij(order, 1, 0, nodePts_0, nodePts_1) +
      c20[k]*Bij(order, 2, 0, nodePts_0, nodePts_1) +
      c11[k]*Bij(order, 1, 1, nodePts_0, nodePts_1) +
      c02[k]*Bij(order, 0, 2, nodePts_0, nodePts_1) +
      c01[k]*Bij(order, 0, 1, nodePts_0, nodePts_1);
  }
  return p11_w;
}

OMEGA_H_DEVICE Vector<3> face_parametricToParent_3d(
    LO const order, LO const old_face, LOs old_ev2v, LOs old_fe2e,
    Reals old_vertCtrlPts, Reals old_edgeCtrlPts, Reals old_faceCtrlPts,
    Real nodePts_0, Real nodePts_1, LOs old_fv2v) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  LO const v0_old_face = old_fv2v[old_face*3 + 0];
  LO const v1_old_face = old_fv2v[old_face*3 + 1];
  LO const v2_old_face = old_fv2v[old_face*3 + 2];
  LO const old_face_e0 = old_fe2e[old_face*3 + 0];
  LO const old_face_e1 = old_fe2e[old_face*3 + 1];
  LO const old_face_e2 = old_fe2e[old_face*3 + 2];
  I8 e0_flip = -1;
  I8 e1_flip = -1;
  I8 e2_flip = -1;
  LO v1 = v1_old_face;
  LO v2 = v2_old_face;
  auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
  auto e0v1 = old_ev2v[old_face_e0*2 + 1];
  auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
  auto e1v1 = old_ev2v[old_face_e1*2 + 1];
  auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
  auto e2v1 = old_ev2v[old_face_e2*2 + 1];
  if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
    e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
  }
  if ((e1v0_old_face == v2) && (e1v1 == v1)) {
    e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
  }
  if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
    e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
  }

  Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
  Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
  Real cz00 = old_vertCtrlPts[v0_old_face*dim + 2];
  Real cx30 = old_vertCtrlPts[v1*dim + 0];
  Real cy30 = old_vertCtrlPts[v1*dim + 1];
  Real cz30 = old_vertCtrlPts[v1*dim + 2];
  Real cx03 = old_vertCtrlPts[v2*dim + 0];
  Real cy03 = old_vertCtrlPts[v2*dim + 1];
  Real cz03 = old_vertCtrlPts[v2*dim + 2];

  auto pts_per_edge = n_edge_pts;
  Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
  Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
  Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
  Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
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

  Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
  Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
  Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
  Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
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

  Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
  Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
  Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
  Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
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

  Real cx11 = old_faceCtrlPts[old_face*dim + 0];
  Real cy11 = old_faceCtrlPts[old_face*dim + 1];
  Real cz11 = old_faceCtrlPts[old_face*dim + 2];
  auto c00 = vector_3(cx00, cy00, cz00);
  auto c10 = vector_3(cx10, cy10, cz10);
  auto c20 = vector_3(cx20, cy20, cz20);
  auto c30 = vector_3(cx30, cy30, cz30);
  auto c21 = vector_3(cx21, cy21, cz21);
  auto c12 = vector_3(cx12, cy12, cz12);
  auto c03 = vector_3(cx03, cy03, cz03);
  auto c02 = vector_3(cx02, cy02, cz02);
  auto c01 = vector_3(cx01, cy01, cz01);
  auto c11 = vector_3(cx11, cy11, cz11);

  Vector<3> p11_w;
  for (LO k = 0; k < 3; ++k) {
    p11_w[k] = c00[k]*Bij(order, 0, 0, nodePts_0, nodePts_1) +
      c10[k]*Bij(order, 1, 0, nodePts_0, nodePts_1) +
      c20[k]*Bij(order, 2, 0, nodePts_0, nodePts_1) +
      c30[k]*Bij(order, 3, 0, nodePts_0, nodePts_1) +
      c21[k]*Bij(order, 2, 1, nodePts_0, nodePts_1) +
      c12[k]*Bij(order, 1, 2, nodePts_0, nodePts_1) +
      c03[k]*Bij(order, 0, 3, nodePts_0, nodePts_1) +
      c02[k]*Bij(order, 0, 2, nodePts_0, nodePts_1) +
      c01[k]*Bij(order, 0, 1, nodePts_0, nodePts_1) +
      c11[k]*Bij(order, 1, 1, nodePts_0, nodePts_1);
  }
  return p11_w;
}

OMEGA_H_DEVICE Vector<3> face_blend_interp_3d(
    LO const order, LO const old_face, LOs old_ev2v, LOs old_fe2e,
    Reals old_vertCtrlPts, Reals old_edgeCtrlPts,
    LOs old_fv2v, Few<Real, 10> weights) {
  OMEGA_H_CHECK(std::abs(weights[9] - 0.0) < EPSILON);
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  LO const v0_old_face = old_fv2v[old_face*3 + 0];
  LO const v1_old_face = old_fv2v[old_face*3 + 1];
  LO const v2_old_face = old_fv2v[old_face*3 + 2];
  LO const old_face_e0 = old_fe2e[old_face*3 + 0];
  LO const old_face_e1 = old_fe2e[old_face*3 + 1];
  LO const old_face_e2 = old_fe2e[old_face*3 + 2];
  I8 e0_flip = -1;
  I8 e1_flip = -1;
  I8 e2_flip = -1;
  LO v1 = v1_old_face;
  LO v2 = v2_old_face;
  auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
  auto e0v1 = old_ev2v[old_face_e0*2 + 1];
  auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
  auto e1v1 = old_ev2v[old_face_e1*2 + 1];
  auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
  auto e2v1 = old_ev2v[old_face_e2*2 + 1];
  if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
    e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
  }
  if ((e1v0_old_face == v2) && (e1v1 == v1)) {
    e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
  }
  if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
    e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
  }

  Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
  Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
  Real cz00 = old_vertCtrlPts[v0_old_face*dim + 2];
  Real cx30 = old_vertCtrlPts[v1*dim + 0];
  Real cy30 = old_vertCtrlPts[v1*dim + 1];
  Real cz30 = old_vertCtrlPts[v1*dim + 2];
  Real cx03 = old_vertCtrlPts[v2*dim + 0];
  Real cy03 = old_vertCtrlPts[v2*dim + 1];
  Real cz03 = old_vertCtrlPts[v2*dim + 2];

  auto pts_per_edge = n_edge_pts;
  Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
  Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
  Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
  Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
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

  Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
  Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
  Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
  Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
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

  Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
  Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
  Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
  Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
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

  auto c00 = vector_3(cx00, cy00, cz00);
  auto c10 = vector_3(cx10, cy10, cz10);
  auto c20 = vector_3(cx20, cy20, cz20);
  auto c30 = vector_3(cx30, cy30, cz30);
  auto c21 = vector_3(cx21, cy21, cz21);
  auto c12 = vector_3(cx12, cy12, cz12);
  auto c03 = vector_3(cx03, cy03, cz03);
  auto c02 = vector_3(cx02, cy02, cz02);
  auto c01 = vector_3(cx01, cy01, cz01);

  Vector<3> p11_w;
  for (LO k = 0; k < 3; ++k) {
    p11_w[k] = c00[k]*weights[0] +
      c30[k]*weights[1] +
      c03[k]*weights[2] +
      c10[k]*weights[3] +
      c20[k]*weights[4] +
      c21[k]*weights[5] +
      c12[k]*weights[6] +
      c02[k]*weights[7] +
      c01[k]*weights[8];
  }
  return p11_w;
}

OMEGA_H_DEVICE Vector<2> face_interpToCtrlPt_2d(
    LO const order, LO const newface, LOs new_ev2v, LOs new_fe2e,
    Reals new_vertCtrlPts, Reals new_edgeCtrlPts, Vector<2> p11, LOs new_fv2v) {
  LO const dim=2;
  LO const pts_per_edge = n_internal_ctrlPts(EDGE, order);
  I8 newface_e0_flip = -1;
  I8 newface_e1_flip = -1;
  I8 newface_e2_flip = -1;
  LO newface_v0 = new_fv2v[newface*3 + 0];
  LO newface_v1 = new_fv2v[newface*3 + 1];
  LO newface_v2 = new_fv2v[newface*3 + 2];
  LO newface_e0 = new_fe2e[newface*3 + 0];
  LO newface_e1 = new_fe2e[newface*3 + 1];
  LO newface_e2 = new_fe2e[newface*3 + 2];
  auto newface_e0v0 = new_ev2v[newface_e0*2 + 0];
  auto newface_e0v1 = new_ev2v[newface_e0*2 + 1];
  auto newface_e1v0 = new_ev2v[newface_e1*2 + 0];
  auto newface_e1v1 = new_ev2v[newface_e1*2 + 1];
  auto newface_e2v0 = new_ev2v[newface_e2*2 + 0];
  auto newface_e2v1 = new_ev2v[newface_e2*2 + 1];
  if ((newface_e0v0 == newface_v1) && (newface_e0v1 == newface_v0)) {
    newface_e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e0v0 == newface_v0) && (newface_e0v1 == newface_v1));
  }
  if ((newface_e1v0 == newface_v2) && (newface_e1v1 == newface_v1)) {
    newface_e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e1v0 == newface_v1) && (newface_e1v1 == newface_v2));
  }
  if ((newface_e2v0 == newface_v0) && (newface_e2v1 == newface_v2)) {
    newface_e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e2v0 == newface_v2) && (newface_e2v1 == newface_v0));
  }

  Real newface_cx00 = new_vertCtrlPts[newface_v0*dim + 0];
  Real newface_cy00 = new_vertCtrlPts[newface_v0*dim + 1];
  Real newface_cx30 = new_vertCtrlPts[newface_v1*dim + 0];
  Real newface_cy30 = new_vertCtrlPts[newface_v1*dim + 1];
  Real newface_cx03 = new_vertCtrlPts[newface_v2*dim + 0];
  Real newface_cy03 = new_vertCtrlPts[newface_v2*dim + 1];

  Real newface_cx10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 0];
  Real newface_cy10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 1];
  Real newface_cx20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (newface_e0_flip > 0) {
    swap2(newface_cx10, newface_cx20);
    swap2(newface_cy10, newface_cy20);
  }

  Real newface_cx21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 0];
  Real newface_cy21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 1];
  Real newface_cx12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (newface_e1_flip > 0) {
    swap2(newface_cx21, newface_cx12);
    swap2(newface_cy21, newface_cy12);
  }

  Real newface_cx02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 0];
  Real newface_cy02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 1];
  Real newface_cx01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (newface_e2_flip > 0) {
    swap2(newface_cx01, newface_cx02);
    swap2(newface_cy01, newface_cy02);
  }
  auto newface_c00 = vector_2(newface_cx00, newface_cy00);
  auto newface_c10 = vector_2(newface_cx10, newface_cy10);
  auto newface_c20 = vector_2(newface_cx20, newface_cy20);
  auto newface_c30 = vector_2(newface_cx30, newface_cy30);
  auto newface_c21 = vector_2(newface_cx21, newface_cy21);
  auto newface_c12 = vector_2(newface_cx12, newface_cy12);
  auto newface_c03 = vector_2(newface_cx03, newface_cy03);
  auto newface_c02 = vector_2(newface_cx02, newface_cy02);
  auto newface_c01 = vector_2(newface_cx01, newface_cy01);

  auto xi_11 = xi_11_cube();
  Vector<2> newface_c11_w;
  for (LO k = 0; k < dim; ++k) {
    newface_c11_w[k] = (p11[k] -
        newface_c00[k]*Bij(order, 0, 0, xi_11[0], xi_11[1]) -
        newface_c10[k]*Bij(order, 1, 0, xi_11[0], xi_11[1]) -
        newface_c20[k]*Bij(order, 2, 0, xi_11[0], xi_11[1]) -
        newface_c30[k]*Bij(order, 3, 0, xi_11[0], xi_11[1]) -
        newface_c21[k]*Bij(order, 2, 1, xi_11[0], xi_11[1]) -
        newface_c12[k]*Bij(order, 1, 2, xi_11[0], xi_11[1]) -
        newface_c03[k]*Bij(order, 0, 3, xi_11[0], xi_11[1]) -
        newface_c02[k]*Bij(order, 0, 2, xi_11[0], xi_11[1]) -
        newface_c01[k]*Bij(order, 0, 1, xi_11[0], xi_11[1]))/
      Bij(order, 1, 1, xi_11[0], xi_11[1]);
  }

  return newface_c11_w;
}

OMEGA_H_DEVICE Vector<3> face_interpToCtrlPt_3d(
    LO const order, LO const newface, LOs new_ev2v, LOs new_fe2e,
    Reals new_vertCtrlPts, Reals new_edgeCtrlPts, Vector<3> p11, LOs new_fv2v) {
  LO const dim=3;
  LO const pts_per_edge = n_internal_ctrlPts(EDGE, order);
  I8 newface_e0_flip = -1;
  I8 newface_e1_flip = -1;
  I8 newface_e2_flip = -1;
  LO newface_v0 = new_fv2v[newface*3 + 0];
  LO newface_v1 = new_fv2v[newface*3 + 1];
  LO newface_v2 = new_fv2v[newface*3 + 2];
  LO newface_e0 = new_fe2e[newface*3 + 0];
  LO newface_e1 = new_fe2e[newface*3 + 1];
  LO newface_e2 = new_fe2e[newface*3 + 2];
  auto newface_e0v0 = new_ev2v[newface_e0*2 + 0];
  auto newface_e0v1 = new_ev2v[newface_e0*2 + 1];
  auto newface_e1v0 = new_ev2v[newface_e1*2 + 0];
  auto newface_e1v1 = new_ev2v[newface_e1*2 + 1];
  auto newface_e2v0 = new_ev2v[newface_e2*2 + 0];
  auto newface_e2v1 = new_ev2v[newface_e2*2 + 1];
  if ((newface_e0v0 == newface_v1) && (newface_e0v1 == newface_v0)) {
    newface_e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e0v0 == newface_v0) && (newface_e0v1 == newface_v1));
  }
  if ((newface_e1v0 == newface_v2) && (newface_e1v1 == newface_v1)) {
    newface_e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e1v0 == newface_v1) && (newface_e1v1 == newface_v2));
  }
  if ((newface_e2v0 == newface_v0) && (newface_e2v1 == newface_v2)) {
    newface_e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e2v0 == newface_v2) && (newface_e2v1 == newface_v0));
  }

  Real newface_cx00 = new_vertCtrlPts[newface_v0*dim + 0];
  Real newface_cy00 = new_vertCtrlPts[newface_v0*dim + 1];
  Real newface_cz00 = new_vertCtrlPts[newface_v0*dim + 2];
  Real newface_cx30 = new_vertCtrlPts[newface_v1*dim + 0];
  Real newface_cy30 = new_vertCtrlPts[newface_v1*dim + 1];
  Real newface_cz30 = new_vertCtrlPts[newface_v1*dim + 2];
  Real newface_cx03 = new_vertCtrlPts[newface_v2*dim + 0];
  Real newface_cy03 = new_vertCtrlPts[newface_v2*dim + 1];
  Real newface_cz03 = new_vertCtrlPts[newface_v2*dim + 2];

  Real newface_cx10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 0];
  Real newface_cy10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 1];
  Real newface_cz10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 2];
  Real newface_cx20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real newface_cz20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (newface_e0_flip > 0) {
    swap2(newface_cx10, newface_cx20);
    swap2(newface_cy10, newface_cy20);
    swap2(newface_cz10, newface_cz20);
  }

  Real newface_cx21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 0];
  Real newface_cy21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 1];
  Real newface_cz21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 2];
  Real newface_cx12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real newface_cz12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (newface_e1_flip > 0) {
    swap2(newface_cx21, newface_cx12);
    swap2(newface_cy21, newface_cy12);
    swap2(newface_cz21, newface_cz12);
  }

  Real newface_cx02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 0];
  Real newface_cy02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 1];
  Real newface_cz02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 2];
  Real newface_cx01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real newface_cz01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (newface_e2_flip > 0) {
    swap2(newface_cx01, newface_cx02);
    swap2(newface_cy01, newface_cy02);
    swap2(newface_cz01, newface_cz02);
  }
  auto newface_c00 = vector_3(newface_cx00, newface_cy00, newface_cz00);
  auto newface_c10 = vector_3(newface_cx10, newface_cy10, newface_cz10);
  auto newface_c20 = vector_3(newface_cx20, newface_cy20, newface_cz20);
  auto newface_c30 = vector_3(newface_cx30, newface_cy30, newface_cz30);
  auto newface_c21 = vector_3(newface_cx21, newface_cy21, newface_cz21);
  auto newface_c12 = vector_3(newface_cx12, newface_cy12, newface_cz12);
  auto newface_c03 = vector_3(newface_cx03, newface_cy03, newface_cz03);
  auto newface_c02 = vector_3(newface_cx02, newface_cy02, newface_cz02);
  auto newface_c01 = vector_3(newface_cx01, newface_cy01, newface_cz01);

  auto xi_11 = xi_11_cube();
  Vector<3> newface_c11_w;
  for (LO k = 0; k < 3; ++k) {
    newface_c11_w[k] = (p11[k] -
        newface_c00[k]*Bij(order, 0, 0, xi_11[0], xi_11[1]) -
        newface_c10[k]*Bij(order, 1, 0, xi_11[0], xi_11[1]) -
        newface_c20[k]*Bij(order, 2, 0, xi_11[0], xi_11[1]) -
        newface_c30[k]*Bij(order, 3, 0, xi_11[0], xi_11[1]) -
        newface_c21[k]*Bij(order, 2, 1, xi_11[0], xi_11[1]) -
        newface_c12[k]*Bij(order, 1, 2, xi_11[0], xi_11[1]) -
        newface_c03[k]*Bij(order, 0, 3, xi_11[0], xi_11[1]) -
        newface_c02[k]*Bij(order, 0, 2, xi_11[0], xi_11[1]) -
        newface_c01[k]*Bij(order, 0, 1, xi_11[0], xi_11[1]))/
      Bij(order, 1, 1, xi_11[0], xi_11[1]);

  }

  return newface_c11_w;
}

LOs create_curved_verts_and_edges_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts,
                                     LOs keys2edges);

void create_curved_faces_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                            LOs keys2prods, LOs keys2edges, LOs keys2old_faces,
                            LOs old_verts2new_verts);

LOs create_curved_verts_and_edges_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts,
                                     LOs keys2edges);

void create_curved_verts_and_edges_3d_p2(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts,
                                     LOs keys2edges, bool transfer_cands, Read<I8> oldedge_is_cand);

void create_curved_faces_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                            LOs keys2prods, LOs keys2edges, LOs keys2old_faces,
                            LOs old_verts2new_verts);

LOs checkValidity_2d(Mesh *new_mesh, LOs new_tris, Int const mesh_dim);

const Real askWorstQuality_2d(Mesh *new_mesh, LOs new_tris, Int const mesh_dim);

Reals askQuality_2d(Mesh *new_mesh, LOs new_tris, Int const mesh_dim);

OMEGA_H_INLINE LO computeTriNodeIndex (LO P, LO i, LO j) {
  LO k = P-i-j;
  if (i == P) return 0;
  if (j == P) return 1;
  if (k == P) return 2;
  if (k == 0) return 2+j;
  if (i == 0) return 2+(P-1)+k;
  if (j == 0) return 2+(P-1)*2+i;
  return k*(P-1)-k*(k-1)/2+j+2*P;
}

OMEGA_H_INLINE LO getTriNodeIndex (LO P, LO i, LO j) {
  return computeTriNodeIndex(P,i,j);
}

template <Int mesh_dim>
OMEGA_H_INLINE Real getTriPartialJacobianDet(Few<Real, 200*mesh_dim> nodes,
//OMEGA_H_INLINE Real getTriPartialJacobianDet(Few<Real, 10*mesh_dim> nodes,
  LO P, LO i1, LO j1, LO i2, LO j2) {
  LO p00 = getTriNodeIndex(P,i1+1,j1);
  LO p01 = getTriNodeIndex(P,i1,j1+1);
  LO p10 = getTriNodeIndex(P,i2+1,j2);
  LO p11 = getTriNodeIndex(P,i2,j2);
  Real ret;
  ret = cross(get_vector<2>(nodes, p01) - get_vector<2>(nodes, p00),
              get_vector<2>(nodes, p11) - get_vector<2>(nodes, p10));
  return ret;
}

template <Int mesh_dim>
OMEGA_H_INLINE Real Nijk(Few<Real, 200*mesh_dim> nodes, LO d, LO I, LO J) {
//OMEGA_H_INLINE Real Nijk(Few<Real, 10*mesh_dim> nodes, LO d, LO I, LO J) {
  Real sum = 0.;
  LO CD = trinomial(2*(d-1), I, J);
  for (LO j1 = 0; j1 <= J; ++j1) {
    auto i1start = max2(0, I+J-j1-(d-1));
    auto i1end = min2(I, d-1-j1);
    for (LO i1 = i1start; i1 <= i1end; ++i1){
      sum += trinomial(d-1, i1, j1)*trinomial(d-1, I-i1, J-j1)
        *getTriPartialJacobianDet<mesh_dim>(nodes, d, i1, j1, I-i1, J-j1);
    }
  }
  return sum*d*d/CD;
}

template <Int n, Int mesh_dim>
OMEGA_H_INLINE Few<Real, n> getTriJacDetNodes(
    LO P, Few<Real, 200*mesh_dim> const& elemNodes) {
    //LO P, Few<Real, 10*mesh_dim> const& elemNodes) {
  Few<Real, n> nodes;//n
  for (LO I = 0; I <= 2*(P-1); ++I) {
    for (LO J = 0; J <= 2*(P-1)-I; ++J) {
        OMEGA_H_CHECK(getTriNodeIndex(2*(P-1),I,J) < nodes.size());
        nodes[getTriNodeIndex(2*(P-1),I,J)] = Nijk<mesh_dim>(elemNodes,P,I,J);
    }
  }
  return nodes;
}

template<Int n>
OMEGA_H_INLINE LO checkMinJacDet(Few<Real, n> const& nodes, LO order) {
  // first 3 vertices
  Real minAcceptable = 0.0;
  for (LO i = 0; i < 3; ++i) {
    if (nodes[i] < minAcceptable) {
      return i+2;
    }
  }

  for (int edge = 0; edge < 3; ++edge) {
    for (int i = 0; i < 2*(order-1)-1; ++i) {
      if (nodes[3+edge*(2*(order-1)-1)+i] < minAcceptable) {
        return 8+edge;
      }
    }
  }

  for (int i = 0; i < (2*order-3)*(2*order-4)/2; ++i) {
    if (nodes[6*(order-1)+i] < minAcceptable) {
      return 14;
    }
  }

  return -1;
}

OMEGA_H_DEVICE Few<Real, 60> collect_tet_pts(
    LO const order, LO const tet, LOs ev2v, LOs rv2v,
    Reals vertCtrlPts, Reals edgeCtrlPts, Reals faceCtrlPts,
    LOs re2e, LOs rf2f) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  auto tet_v0 = rv2v[tet*4 + 0];
  auto tet_v1 = rv2v[tet*4 + 1];
  auto tet_v2 = rv2v[tet*4 + 2];
  auto tet_v3 = rv2v[tet*4 + 3];
  auto c000 = get_vector<3>(vertCtrlPts, tet_v0);
  auto c300 = get_vector<3>(vertCtrlPts, tet_v1);
  auto c030 = get_vector<3>(vertCtrlPts, tet_v2);
  auto c003 = get_vector<3>(vertCtrlPts, tet_v3);

  auto tet_e0 = re2e[tet*6 + 0];
  auto tet_e1 = re2e[tet*6 + 1];
  auto tet_e2 = re2e[tet*6 + 2];
  auto tet_e3 = re2e[tet*6 + 3];
  auto tet_e4 = re2e[tet*6 + 4];
  auto tet_e5 = re2e[tet*6 + 5];

  auto tet_f0 = rf2f[tet*4 + 0];
  auto tet_f1 = rf2f[tet*4 + 1];
  auto tet_f2 = rf2f[tet*4 + 2];
  auto tet_f3 = rf2f[tet*4 + 3];
  auto c110 = get_vector<3>(faceCtrlPts, tet_f0);
  auto c101 = get_vector<3>(faceCtrlPts, tet_f1);
  auto c111 = get_vector<3>(faceCtrlPts, tet_f2);
  auto c011 = get_vector<3>(faceCtrlPts, tet_f3);

  LO e0_flip = -1;
  LO e1_flip = -1;
  LO e2_flip = -1;
  LO e3_flip = -1;
  LO e4_flip = -1;
  LO e5_flip = -1;
  auto e0v0 = ev2v[tet_e0*2 + 0];
  auto e0v1 = ev2v[tet_e0*2 + 1];
  auto e1v0 = ev2v[tet_e1*2 + 0];
  auto e1v1 = ev2v[tet_e1*2 + 1];
  auto e2v0 = ev2v[tet_e2*2 + 0];
  auto e2v1 = ev2v[tet_e2*2 + 1];
  auto e3v0 = ev2v[tet_e3*2 + 0];
  auto e3v1 = ev2v[tet_e3*2 + 1];
  auto e4v0 = ev2v[tet_e4*2 + 0];
  auto e4v1 = ev2v[tet_e4*2 + 1];
  auto e5v0 = ev2v[tet_e5*2 + 0];
  auto e5v1 = ev2v[tet_e5*2 + 1];
  e0_flip = edge_is_flip(e0v0, e0v1, tet_v0, tet_v1);
  e1_flip = edge_is_flip(e1v0, e1v1, tet_v1, tet_v2);
  e2_flip = edge_is_flip(e2v0, e2v1, tet_v2, tet_v0);
  e3_flip = edge_is_flip(e3v0, e3v1, tet_v0, tet_v3);
  e4_flip = edge_is_flip(e4v0, e4v1, tet_v1, tet_v3);
  e5_flip = edge_is_flip(e5v0, e5v1, tet_v2, tet_v3);

  auto pts_per_edge = n_edge_pts;
  Real cx100 = edgeCtrlPts[tet_e0*pts_per_edge*dim + 0];
  Real cy100 = edgeCtrlPts[tet_e0*pts_per_edge*dim + 1];
  Real cz100 = edgeCtrlPts[tet_e0*pts_per_edge*dim + 2];
  Real cx200 = edgeCtrlPts[tet_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy200 = edgeCtrlPts[tet_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz200 = edgeCtrlPts[tet_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e0_flip > 0) {
    swap2(cx100, cx200);
    swap2(cy100, cy200);
    swap2(cz100, cz200);
  }
  Real cx210 = edgeCtrlPts[tet_e1*pts_per_edge*dim + 0];
  Real cy210 = edgeCtrlPts[tet_e1*pts_per_edge*dim + 1];
  Real cz210 = edgeCtrlPts[tet_e1*pts_per_edge*dim + 2];
  Real cx120 = edgeCtrlPts[tet_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy120 = edgeCtrlPts[tet_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz120 = edgeCtrlPts[tet_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e1_flip > 0) {
    swap2(cx210, cx120);
    swap2(cy210, cy120);
    swap2(cz210, cz120);
  }
  Real cx020 = edgeCtrlPts[tet_e2*pts_per_edge*dim + 0];
  Real cy020 = edgeCtrlPts[tet_e2*pts_per_edge*dim + 1];
  Real cz020 = edgeCtrlPts[tet_e2*pts_per_edge*dim + 2];
  Real cx010 = edgeCtrlPts[tet_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy010 = edgeCtrlPts[tet_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz010 = edgeCtrlPts[tet_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e2_flip > 0) {
    swap2(cx020, cx010);
    swap2(cy020, cy010);
    swap2(cz020, cz010);
  }
  Real cx001 = edgeCtrlPts[tet_e3*pts_per_edge*dim + 0];
  Real cy001 = edgeCtrlPts[tet_e3*pts_per_edge*dim + 1];
  Real cz001 = edgeCtrlPts[tet_e3*pts_per_edge*dim + 2];
  Real cx002 = edgeCtrlPts[tet_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy002 = edgeCtrlPts[tet_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz002 = edgeCtrlPts[tet_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e3_flip > 0) {
    swap2(cx002, cx001);
    swap2(cy002, cy001);
    swap2(cz002, cz001);
  }
  Real cx201 = edgeCtrlPts[tet_e4*pts_per_edge*dim + 0];
  Real cy201 = edgeCtrlPts[tet_e4*pts_per_edge*dim + 1];
  Real cz201 = edgeCtrlPts[tet_e4*pts_per_edge*dim + 2];
  Real cx102 = edgeCtrlPts[tet_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy102 = edgeCtrlPts[tet_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz102 = edgeCtrlPts[tet_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e4_flip > 0) {
    swap2(cx102, cx201);
    swap2(cy102, cy201);
    swap2(cz102, cz201);
  }
  Real cx021 = edgeCtrlPts[tet_e5*pts_per_edge*dim + 0];
  Real cy021 = edgeCtrlPts[tet_e5*pts_per_edge*dim + 1];
  Real cz021 = edgeCtrlPts[tet_e5*pts_per_edge*dim + 2];
  Real cx012 = edgeCtrlPts[tet_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy012 = edgeCtrlPts[tet_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz012 = edgeCtrlPts[tet_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e5_flip > 0) {
    swap2(cx012, cx021);
    swap2(cy012, cy021);
    swap2(cz012, cz021);
  }
  auto c100 = vector_3(cx100, cy100, cz100);
  auto c200 = vector_3(cx200, cy200, cz200);
  auto c210 = vector_3(cx210, cy210, cz210);
  auto c120 = vector_3(cx120, cy120, cz120);
  auto c020 = vector_3(cx020, cy020, cz020);
  auto c010 = vector_3(cx010, cy010, cz010);
  auto c001 = vector_3(cx001, cy001, cz001);
  auto c002 = vector_3(cx002, cy002, cz002);
  auto c102 = vector_3(cx102, cy102, cz102);
  auto c201 = vector_3(cx201, cy201, cz201);
  auto c012 = vector_3(cx012, cy012, cz012);
  auto c021 = vector_3(cx021, cy021, cz021);

  Few<Real,60> tet_pts;
  for (LO j = 0; j < dim; ++j) {
    tet_pts[0*dim + j] = c000[j];
    tet_pts[1*dim + j] = c300[j];
    tet_pts[2*dim + j] = c030[j];
    tet_pts[3*dim + j] = c003[j];
    tet_pts[4*dim + j] = c100[j];
    tet_pts[5*dim + j] = c200[j];
    tet_pts[6*dim + j] = c210[j];
    tet_pts[7*dim + j] = c120[j];
    tet_pts[8*dim + j] = c020[j];
    tet_pts[9*dim + j] = c010[j];
    tet_pts[10*dim + j] = c001[j];
    tet_pts[11*dim + j] = c002[j];
    tet_pts[12*dim + j] = c201[j];
    tet_pts[13*dim + j] = c102[j];
    tet_pts[14*dim + j] = c012[j];
    tet_pts[15*dim + j] = c021[j];
    tet_pts[16*dim + j] = c110[j];
    tet_pts[17*dim + j] = c101[j];
    tet_pts[18*dim + j] = c111[j];
    tet_pts[19*dim + j] = c011[j];
  }
  return tet_pts;
}

#define OMEGA_H_EXPL_INST_DECL(T)                                              
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

} // namespace Omega_h

#endif
