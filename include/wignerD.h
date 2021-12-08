#pragma once

/*
  Library for computing Wigner D matrices

    D = [D_{m1,m2} ^l(a,b,c)]_{m1,m2=-l}^l

  the with matrix elements

    D_{m1,m2} ^j(a,b,c) = <l m1| exp(I a Jz) exp(I b Jy) exp(I c Jz)| lm2 >
                        = exp[I (a m1 + c m2)] d^l_{m1,m2}(b)

  where d denoted Wigner's (small) d-matrix

    d^l_{m1,m2}(b) = D_{m1,m2}^lj(0,b,0)

  with
    m1, m2, l are integers,
    m1, m2 in [-l,l], l >= 0, Jz = -I d/d phi

  and basis functions are spherical harmonics:

    <theta, phi | l m > = Y_n^l (theta,phi)

    l - azimuthal quantum number/angular momentum quantum number
    m - magnetic quantum number

  normalized as <lm|lm> = 1

  Definition is aligned with Wolfram mathematica and functional  basis

  Author: Martin Horvat, November 2021
*/

#include <cmath>
#include <complex>
#include <algorithm>

#include "matrix.h"

// assuming n >=1
template<int n, class T>
void calc_sincos(const T & x, T s[n], T c[n]){

  #if 0
  for (int i = 0, j = 1; i < n; ++i, ++j) {
    s[i] = std::sin(j*x);
    c[i] = std::cos(j*x);
  }
  #else   // recursive calculations of sin and cos

  s[0] = std::sin(x);
  c[0] = std::cos(x);

  for (int i = 1; i < n; ++i) {
    s[i] = c[i-1]*s[0] + s[i-1]*c[0];
    c[i] = c[i-1]*c[0] - s[i-1]*s[0];
  }
  #endif
}

template <typename T>
void tqli(T *d, T *e, int n, T **z) {

  auto sign_ = [](const T & a, const T &b)
    { return (b >= 0 ? std::abs(a) : -std::abs(a)); };

  int m, l, iter, i, k;

  T s, r, p, g, f, dd, c, b;

  for (i = 1; i < n; ++i) e[i-1] = e[i];
  e[n-1] = 0;

  for (l = 0; l < n; ++l) {

    iter = 0;

    do {

      for (m = l; m < n-1; ++m) {
        dd = std::abs(d[m]) + std::abs(d[m + 1]);
        if ((std::abs(e[m]) + dd) == dd) break;
      }

      if (m != l) {

        if (iter++ == 30) std::cerr << "Too many iterations in tqli\n";

        g =(d[l + 1] - d[l])/(2*e[l]);
        r = std::hypot(g, T(1));
        g = d[m] - d[l] + e[l]/(g + sign_(r, g));
        s = c = T(1);
        p = T(0);

        for (i = m-1; i >= l; --i) {

          f = s*e[i];
          b = c*e[i];
          e[i+1] = (r = std::hypot(f, g));

          if (r == T(0)) {
            d[i+1] -= p;
            e[m] = T(0);
            break;
          }

          s = f/r;
          c = g/r;
          g = d[i+1] - p;
          r = (d[i] - g)*s + 2*c*b;
          d[i+1] = g + (p = s*r);
          g = c*r - b;

          for (k = 0; k < n; ++k) {
            f = z[k][i+1];
            z[k][i+1] = s*z[k][i] + c*f;
            z[k][i] = c*z[k][i] - s*f;
          }
        }

        if (r == T(0) && i >= l) continue;

        d[l] -= p;
        e[l] = g;
        e[m] = T(0);
      }
    } while (m != l);
  }
}



/*
  Calculate Wigner's (small) d-matrix

    d = [d_{m1,m2}^j]_{m1,m2=-l}^l

  at j with matrix elements

      d_{m1,m2}^l = <lm1|exp(I a Jy) |lm2>

  Input:

    a : real    - angle of rotation around y axis
    l : integer - azimuthal quantum number

  Output:
    R[2j+1][2j+1] - Wigner's d matrix is real

*/

template <class T>
void wigner_d(const T & a, int l, T **R) {

  int
    dim = 2*l + 1,
    size = dim*dim;

  if (a == 0) {
    std::fill(R[0], R[0] + size, T(0));
    for (int i = 0; i < dim; ++i) R[i][i] = 1;
    return;
  }

  using F = double;

  switch (l) {

    case 0:
      R[0][0] = 1;
    break;

    case 1:
      {
        F s = std::sin(F(a)), c = std::cos(F(a)), t;

        R[0][0] = R[2][2] = (1 + c)/2;
        R[0][2] = R[2][0] = (1 - c)/2;
        R[1][1] = c;
        t = s/std::sqrt(F(2));
        R[0][1] = R[1][2] = -t;
        R[1][0] = R[2][1] =  t;
      }
    break;

    case 2:
      {
        const int n = 2;
        F s[n], c[n];
        calc_sincos<n>(F(a), s, c);

        F t[9] = {
          (3 - 4*c[0] + c[1])/8., (2*s[0] - s[1])/4.,
          std::sqrt(F(6))*(1 - c[1])/8., (c[0] - c[1])/2.,
          (2*s[0] + s[1])/4., (std::sqrt(F(1.5))*s[1])/2.,
          (1 + 3*c[1])/4.,(c[0] + c[1])/2.,(3 + 4*c[0] + c[1])/8.
        };

       R[4][4] = R[0][0] = t[8], R[4][3] = R[1][0] = t[4],
       R[3][4] = R[0][1] = -t[4],
       R[4][2] = R[2][4] = R[2][0] = R[0][2] = t[2],
       R[4][1] = R[3][0] = t[1], R[1][4] = R[0][3] = -t[1],
       R[4][0] = R[0][4] = t[0], R[3][3] = R[1][1] = t[7],
       R[3][2] = R[2][1] = t[5], R[2][3] = R[1][2] = -t[5],
       R[3][1] = R[1][3] = t[3], R[2][2] = t[6];
      }
    break;

    case 3:
      {
        const int n = 3;
        F s[n], c[n];
        calc_sincos<n>(F(a), s, c);

        F t[16] = {
            (10 - 15*c[0] + 6*c[1] - c[2])/32.,
            (5*std::sqrt(F(6))*s[0] - 4*std::sqrt(F(6))*s[1] + std::sqrt(F(6))*s[2])/32.,
            (2*std::sqrt(F(15)) - std::sqrt(F(15))*c[0] - 2*std::sqrt(F(15))*c[1] + std::sqrt(F(15))*c[2])/32.,
            (5*c[0] - 8*c[1] + 3*c[2])/16.,
            (3*std::sqrt(F(5))*s[0] - std::sqrt(F(5))*s[2])/16.,
            (std::sqrt(F(5))*s[0] + 4*std::sqrt(F(5))*s[1] - 3*std::sqrt(F(5))*s[2])/(16.*std::sqrt(F(2))),
            (2*std::sqrt(F(15)) + std::sqrt(F(15))*c[0] - 2*std::sqrt(F(15))*c[1] - std::sqrt(F(15))*c[2])/32.,
            (std::sqrt(F(30))*c[0] - std::sqrt(F(30))*c[2])/16.,(6 - c[0] + 10*c[1] - 15*c[2])/32.,
            (5*std::sqrt(F(6))*s[0] + 4*std::sqrt(F(6))*s[1] + std::sqrt(F(6))*s[2])/32.,
            -0.0625*(std::sqrt(F(5))*s[0] - 4*std::sqrt(F(5))*s[1] - 3*std::sqrt(F(5))*s[2])/std::sqrt(F(2)),
            (std::sqrt(F(3))*s[0] + 5*std::sqrt(F(3))*s[2])/16.,
            (3*c[0] + 5*c[2])/8.,
            (6 + c[0] + 10*c[1] + 15*c[2])/32.,
            (5*c[0] + 8*c[1] + 3*c[2])/16.,(10 + 15*c[0] + 6*c[1] + c[2])/32.
          };


        R[6][6] = R[0][0] = t[15], R[6][5] = R[1][0] = t[9],
        R[5][6] = R[0][1] = -t[9],R[6][4] = R[4][6] = R[2][0] = R[0][2] = t[6],
        R[6][3] = R[3][0] = t[4],R[3][6] = R[0][3] = -t[4],
        R[6][2] = R[4][0] = R[2][6] = R[0][4] = t[2],R[6][1] = R[5][0] = t[1],
        R[1][6] = R[0][5] = -t[1],R[6][0] = R[0][6] = t[0],R[5][5] = R[1][1] = t[14],
        R[5][4] = R[2][1] = t[10],R[4][5] = R[1][2] = -t[10],
        R[5][3] = R[3][5] = R[3][1] = R[1][3] = t[7],R[5][2] = R[4][1] = t[5],
        R[2][5] = R[1][4] = -t[5],R[5][1] = R[1][5] = t[3],R[4][4] = R[2][2] = t[13],
        R[4][3] = R[3][2] = t[11],R[3][4] = R[2][3] = -t[11],R[4][2] = R[2][4] = t[8],
        R[3][3] = t[12];

      }
    break;
    case 4:
      {
        const int n = 4;
        F s[n], c[n];
        calc_sincos<n>(F(a), s, c);

        F
          t[25] = {(35 - 56*c[0] + 28*c[1] - 8*c[2] + c[3])/128.,
            (14*std::sqrt(F(2))*s[0] - 14*std::sqrt(F(2))*s[1] + 6*std::sqrt(F(2))*s[2] - std::sqrt(F(2))*s[3])/64.,
            (5*std::sqrt(F(7)) - 4*std::sqrt(F(7))*c[0] - 4*std::sqrt(F(7))*c[1] + 4*std::sqrt(F(7))*c[2] - std::sqrt(F(7))*c[3])/64.,
            (7*c[0] - 14*c[1] + 9*c[2] - 2*c[3])/32.,
            (6*std::sqrt(F(14))*s[0] - 2*std::sqrt(F(14))*s[1] - 2*std::sqrt(F(14))*s[2] + std::sqrt(F(14))*s[3])/64.,
            (std::sqrt(F(7))*s[0] + 2*std::sqrt(F(7))*s[1] - 3*std::sqrt(F(7))*s[2] + std::sqrt(F(7))*s[3])/(16.*std::sqrt(F(2))),
            (3*std::sqrt(F(70)) - 4*std::sqrt(F(70))*c[1] + std::sqrt(F(70))*c[3])/128.,
            (3*std::sqrt(F(7))*c[0] - 2*std::sqrt(F(7))*c[1] - 3*std::sqrt(F(7))*c[2] + 2*std::sqrt(F(7))*c[3])/32.,
            (5 - 2*c[0] + 4*c[1] - 14*c[2] + 7*c[3])/32.,
            (6*std::sqrt(F(14))*s[0] + 2*std::sqrt(F(14))*s[1] - 2*std::sqrt(F(14))*s[2] - std::sqrt(F(14))*s[3])/64.,
            (2*std::sqrt(F(35))*s[1] - std::sqrt(F(35))*s[3])/32.,
            (3*std::sqrt(F(2))*s[0] + 2*std::sqrt(F(2))*s[1] + 7*std::sqrt(F(2))*s[2] - 7*std::sqrt(F(2))*s[3])/32.,
            (5*std::sqrt(F(7)) + 4*std::sqrt(F(7))*c[0] - 4*std::sqrt(F(7))*c[1] - 4*std::sqrt(F(7))*c[2] - std::sqrt(F(7))*c[3])/64.,
            (3*std::sqrt(F(7))*c[0] + 2*std::sqrt(F(7))*c[1] - 3*std::sqrt(F(7))*c[2] - 2*std::sqrt(F(7))*c[3])/32.,
            (3*std::sqrt(F(5)) + 4*std::sqrt(F(5))*c[1] - 7*std::sqrt(F(5))*c[3])/(32.*std::sqrt(F(2))),
            (9*c[0] - 2*c[1] + 7*c[2] - 14*c[3])/32.,(9 + 20*c[1] + 35*c[3])/64.,
            (9*c[0] + 2*c[1] + 7*c[2] + 14*c[3])/32.,
            (14*std::sqrt(F(2))*s[0] + 14*std::sqrt(F(2))*s[1] + 6*std::sqrt(F(2))*s[2] + std::sqrt(F(2))*s[3])/64.,
            (5 + 2*c[0] + 4*c[1] + 14*c[2] + 7*c[3])/32.,
            -0.0625*(std::sqrt(F(7))*s[0] - 2*std::sqrt(F(7))*s[1] - 3*std::sqrt(F(7))*s[2] - std::sqrt(F(7))*s[3])/std::sqrt(F(2)),
            (3*std::sqrt(F(2))*s[0] - 2*std::sqrt(F(2))*s[1] + 7*std::sqrt(F(2))*s[2] + 7*std::sqrt(F(2))*s[3])/32.,
            (2*std::sqrt(F(5))*s[1] + 7*std::sqrt(F(5))*s[3])/32.,
            (7*c[0] + 14*c[1] + 9*c[2] + 2*c[3])/32.,
            (35 + 56*c[0] + 28*c[1] + 8*c[2] + c[3])/128.};

         R[8][8] = R[0][0] = t[24],R[8][7] = R[1][0] = t[18],R[7][8] = R[0][1] = -t[18],
         R[8][6] = R[6][8] = R[2][0] = R[0][2] = t[12],R[8][5] = R[3][0] = t[9],
         R[5][8] = R[0][3] = -t[9],R[8][4] = R[4][8] = R[4][0] = R[0][4] = t[6],
         R[8][3] = R[5][0] = t[4],R[3][8] = R[0][5] = -t[4],R[8][2] = R[6][0] = R[2][8] = R[0][6] = t[2],
         R[8][1] = R[7][0] = t[1],R[1][8] = R[0][7] = -t[1],R[8][0] = R[0][8] = t[0],
         R[7][7] = R[1][1] = t[23],R[7][6] = R[2][1] = t[20],R[6][7] = R[1][2] = -t[20],
         R[7][5] = R[5][7] = R[3][1] = R[1][3] = t[13],R[7][4] = R[4][1] = t[10],
         R[4][7] = R[1][4] = -t[10],R[7][3] = R[5][1] = R[3][7] = R[1][5] = t[7],
         R[7][2] = R[6][1] = t[5],R[2][7] = R[1][6] = -t[5],R[7][1] = R[1][7] = t[3],
         R[6][6] = R[2][2] = t[19],R[6][5] = R[3][2] = t[21],R[5][6] = R[2][3] = -t[21],
         R[6][4] = R[4][6] = R[4][2] = R[2][4] = t[14],R[6][3] = R[5][2] = t[11],
         R[3][6] = R[2][5] = -t[11],R[6][2] = R[2][6] = t[8],R[5][5] = R[3][3] = t[17],
         R[5][4] = R[4][3] = t[22],R[4][5] = R[3][4] = -t[22],R[5][3] = R[3][5] = t[15],
         R[4][4] = t[16];
      }
      break;
      case 5:
      {
        const int n = 5;
        F s[n], c[n];
        calc_sincos<n>(F(a), s, c);


        F t[36]= {(126 - 210*c[0] + 120*c[1] - 45*c[2] + 10*c[3] - c[4])/512.,
          (42*std::sqrt(F(10))*s[0] - 48*std::sqrt(F(10))*s[1] + 27*std::sqrt(F(10))*s[2] - 8*std::sqrt(F(10))*s[3] + std::sqrt(F(10))*s[4])/512.,
          (-3*(-14*std::sqrt(F(5)) + 14*std::sqrt(F(5))*c[0] + 8*std::sqrt(F(5))*c[1] - 13*std::sqrt(F(5))*c[2] + 6*std::sqrt(F(5))*c[3] - std::sqrt(F(5))*c[4]))/512.,
          (42*c[0] - 96*c[1] + 81*c[2] - 32*c[3] + 5*c[4])/256.,
          (14*std::sqrt(F(30))*s[0] - 8*std::sqrt(F(30))*s[1] - 3*std::sqrt(F(30))*s[2] + 4*std::sqrt(F(30))*s[3] - std::sqrt(F(30))*s[4])/256.,
          (3*(14*std::sqrt(F(2))*s[0] + 16*std::sqrt(F(2))*s[1] - 39*std::sqrt(F(2))*s[2] + 24*std::sqrt(F(2))*s[3] - 5*std::sqrt(F(2))*s[4]))/512.,
          (6*std::sqrt(F(210)) - 2*std::sqrt(F(210))*c[0] - 8*std::sqrt(F(210))*c[1] + 3*std::sqrt(F(210))*c[2] + 2*std::sqrt(F(210))*c[3] - std::sqrt(F(210))*c[4])/512.,
          (14*std::sqrt(F(3))*c[0] - 16*std::sqrt(F(3))*c[1] - 9*std::sqrt(F(3))*c[2] + 16*std::sqrt(F(3))*c[3] - 5*std::sqrt(F(3))*c[4])/128.,
          (70 - 42*c[0] + 24*c[1] - 169*c[2] + 162*c[3] - 45*c[4])/512.,
          (3*(10*std::sqrt(F(7))*s[0] - 5*std::sqrt(F(7))*s[2] + std::sqrt(F(7))*s[4]))/256.,
          (2*std::sqrt(F(21))*s[0] + 16*std::sqrt(F(21))*s[1] - 9*std::sqrt(F(21))*s[2] - 8*std::sqrt(F(21))*s[3] + 5*std::sqrt(F(21))*s[4])/256.,
          (std::sqrt(F(3))*(14*std::sqrt(F(2))*s[0] + 8*std::sqrt(F(2))*s[1] + 13*std::sqrt(F(2))*s[2] - 36*std::sqrt(F(2))*s[3] + 15*std::sqrt(F(2))*s[4]))/256.,
          (6*std::sqrt(F(210)) + 2*std::sqrt(F(210))*c[0] - 8*std::sqrt(F(210))*c[1] - 3*std::sqrt(F(210))*c[2] + 2*std::sqrt(F(210))*c[3] + std::sqrt(F(210))*c[4])/512.,
          (3*(2*std::sqrt(F(70))*c[0] - 3*std::sqrt(F(70))*c[2] + std::sqrt(F(70))*c[4]))/256.,
          -0.001953125*(std::sqrt(F(3))*(-10*std::sqrt(F(14)) + 2*std::sqrt(F(14))*c[0] - 8*std::sqrt(F(14))*c[1] + 13*std::sqrt(F(14))*c[2] + 18*std::sqrt(F(14))*c[3] - 15*std::sqrt(F(14))*c[4])),
          (14*c[0] - 8*c[1] + 3*c[2] - 24*c[3] + 15*c[4])/64.,
          (14*std::sqrt(F(30))*s[0] + 8*std::sqrt(F(30))*s[1] - 3*std::sqrt(F(30))*s[2] - 4*std::sqrt(F(30))*s[3] - std::sqrt(F(30))*s[4])/256.,
          (-2*std::sqrt(F(21))*s[0] + 16*std::sqrt(F(21))*s[1] + 9*std::sqrt(F(21))*s[2] - 8*std::sqrt(F(21))*s[3] - 5*std::sqrt(F(21))*s[4])/256.,
          (std::sqrt(F(5))*(6*std::sqrt(F(7))*s[0] + 13*std::sqrt(F(7))*s[2] - 9*std::sqrt(F(7))*s[4]))/256.,
          (2*std::sqrt(F(7))*s[0] + 8*std::sqrt(F(7))*s[1] + 3*std::sqrt(F(7))*s[2] + 12*std::sqrt(F(7))*s[3] - 15*std::sqrt(F(7))*s[4])/128.,
          (3*(14*std::sqrt(F(5)) + 14*std::sqrt(F(5))*c[0] - 8*std::sqrt(F(5))*c[1] - 13*std::sqrt(F(5))*c[2] - 6*std::sqrt(F(5))*c[3] - std::sqrt(F(5))*c[4]))/512.,
          (30*c[0] + 35*c[2] + 63*c[4])/128.,
          (30 + 2*c[0] + 56*c[1] + 21*c[2] + 42*c[3] + 105*c[4])/256.,
          (14*std::sqrt(F(3))*c[0] + 16*std::sqrt(F(3))*c[1] - 9*std::sqrt(F(3))*c[2] - 16*std::sqrt(F(3))*c[3] - 5*std::sqrt(F(3))*c[4])/128.,
          (14*c[0] + 8*c[1] + 3*c[2] + 24*c[3] + 15*c[4])/64.,
          (std::sqrt(F(3))*(10*std::sqrt(F(14)) + 2*std::sqrt(F(14))*c[0] + 8*std::sqrt(F(14))*c[1] + 13*std::sqrt(F(14))*c[2] - 18*std::sqrt(F(14))*c[3] - 15*std::sqrt(F(14))*c[4]))/512.,
          (std::sqrt(F(2.5))*(2*std::sqrt(F(21))*c[0] + std::sqrt(F(21))*c[2] - 3*std::sqrt(F(21))*c[4]))/64.,
          (30 - 2*c[0] + 56*c[1] - 21*c[2] + 42*c[3] - 105*c[4])/256.,
          (70 + 42*c[0] + 24*c[1] + 169*c[2] + 162*c[3] + 45*c[4])/512.,
          (42*std::sqrt(F(10))*s[0] + 48*std::sqrt(F(10))*s[1] + 27*std::sqrt(F(10))*s[2] + 8*std::sqrt(F(10))*s[3] + std::sqrt(F(10))*s[4])/512.,
          (-3*(14*std::sqrt(F(2))*s[0] - 16*std::sqrt(F(2))*s[1] - 39*std::sqrt(F(2))*s[2] - 24*std::sqrt(F(2))*s[3] - 5*std::sqrt(F(2))*s[4]))/512.,
          (42*c[0] + 96*c[1] + 81*c[2] + 32*c[3] + 5*c[4])/256.,
          (std::sqrt(F(5))*(2*std::sqrt(F(6))*s[0] + 7*std::sqrt(F(6))*s[2] + 21*std::sqrt(F(6))*s[4]))/256.,
          (-2*std::sqrt(F(7))*s[0] + 8*std::sqrt(F(7))*s[1] - 3*std::sqrt(F(7))*s[2] + 12*std::sqrt(F(7))*s[3] + 15*std::sqrt(F(7))*s[4])/128.,
          (std::sqrt(F(3))*(14*std::sqrt(F(2))*s[0] - 8*std::sqrt(F(2))*s[1] + 13*std::sqrt(F(2))*s[2] + 36*std::sqrt(F(2))*s[3] + 15*std::sqrt(F(2))*s[4]))/256.,
          (126 + 210*c[0] + 120*c[1] + 45*c[2] + 10*c[3] + c[4])/512.};

          R[10][10] = R[0][0] = t[35], R[10][9] = R[1][0] = t[29],
          R[9][10] = R[0][1] = -t[29], R[10][8] = R[8][10] = R[2][0] = R[0][2] = t[20],
          R[10][7] = R[3][0] = t[16], R[7][10] = R[0][3] = -t[16],
          R[10][6] = R[6][10] = R[4][0] = R[0][4] = t[12], R[10][5] = R[5][0] = t[9],
          R[5][10] = R[0][5] = -t[9], R[10][4] = R[6][0] = R[4][10] = R[0][6] = t[6],
          R[10][3] = R[7][0] = t[4], R[3][10] = R[0][7] = -t[4],
          R[10][2] = R[8][0] = R[2][10] = R[0][8] = t[2],
          R[10][1] = R[9][0] = t[1], R[1][10] = R[0][9] = -t[1],
          R[10][0] = R[0][10] = t[0], R[9][9] = R[1][1] = t[31],
          R[9][8] = R[2][1] = t[30], R[8][9] = R[1][2] = -t[30],
          R[9][7] = R[7][9] = R[3][1] = R[1][3] = t[23],
          R[9][6] = R[4][1] = t[17],R[6][9] = R[1][4] = -t[17],
          R[9][5] = R[5][9] = R[5][1] = R[1][5] = t[13],
          R[9][4] = R[6][1] = t[10], R[4][9] = R[1][6] = -t[10],
          R[9][3] = R[7][1] = R[3][9] = R[1][7] = t[7],
          R[9][2] = R[8][1] = t[5], R[2][9] = R[1][8] = -t[5],
          R[9][1] = R[1][9] = t[3], R[8][8] = R[2][2] = t[28],
          R[8][7] = R[3][2] = t[34], R[7][8] = R[2][3] = -t[34],
          R[8][6] = R[6][8] = R[4][2] = R[2][4] = t[25],
          R[8][5] = R[5][2] = t[18], R[5][8] = R[2][5] = -t[18],
          R[8][4] = R[6][2] = R[4][8] = R[2][6] = t[14],
          R[8][3] = R[7][2] = t[11], R[3][8] = R[2][7] = -t[11],
          R[8][2] = R[2][8] = t[8], R[7][7] = R[3][3] = t[24],
          R[7][6] = R[4][3] = t[33], R[6][7] = R[3][4] = -t[33],
          R[7][5] = R[5][7] = R[5][3] = R[3][5] = t[26],
          R[7][4] = R[6][3] = t[19], R[4][7] = R[3][6] = -t[19],
          R[7][3] = R[3][7] = t[15], R[6][6] = R[4][4] = t[22],
          R[6][5] = R[5][4] = t[32],R[5][6] = R[4][5] = -t[32],
          R[6][4] = R[4][6] = t[27], R[5][5] = t[21];
      }
      break;
    default:

      //
      // compute A = exp(i a Jx): Jx is sym. tridigonal, Jx = O diag(d) O^T,
      //           = O exp(I d diag(a)) O^T
      //

      F **O = matrix <F> (dim, dim),
        *d = new F [dim],
        *e = new F [dim];

      std::fill(d, d + dim, F(0));
      for (int m = 0, n = dim; m < dim; ++m, --n) e[m] = std::sqrt(F(m*n))/2;

      std::fill(O[0], O[0] + size, F(0));
      for (int m = 0; m < dim; ++m) O[m][m] = F(1);

      tqli(d, e, dim, O);

      F t;
      std::complex<F> **c = matrix<std::complex<F>>(2,dim);

      // calculate propagation

      for (int i = 0, k = 4 + (-l % 4); i < dim; ++i, ++k) {
        t = F(a)*d[i];
        c[0][i] = std::complex<F> (std::cos(t), std::sin(t));

        if (k >= 4) k = 0;

        switch (k) {
          case 0: c[1][i] = F(1); break;
          case 1: c[1][i] = std::complex<F> (T(0),F(1)); break;
          case 2: c[1][i] = F(-1); break;
          case 3: c[1][i] = std::complex<F> (T(0),F(-1));
        }
      }

      delete [] e;
      delete [] d;

      //
      // compute exp(i a Jy) =  exp(i pi/2 Jz)  A exp(-i pi/2 Jz)
      //

      std::complex<F> sum;

      for (int m = 0; m < dim; ++m)
        for (int n = 0; n <= m; ++n){

          sum = 0;
          for (int k = 0; k < dim; ++k) sum += O[m][k]*O[n][k]*c[0][k];

          R[m][n] = std::real(sum*std::conj(c[1][m])*c[1][n]);
          if (m != n) R[n][m] = ((m - n) % 2 ? -R[m][n] : R[m][n]);
        }

      free_matrix(O);
      free_matrix(c);
  }
}


/*
  Calculate

    v[l + m] = exp(I m x)      m = -j, ..., j

  Input:
    x - argument
    l - range of m

  Output:
    v = [exp(-I l x), ..., exp(I l x)]
*/
template <class T>
void calc_expI(const T &x, int l, std::complex<T> *v){

  using C = std::complex<T>;

  C one(T(1));

  if (x == 0)
    std::fill(v, v + 2*l + 1, one);
  else {
  #if 0
    T s, c, v[l] = one;
    C *vm = v + l - 1, *vp = v + l + 1;
    for (int i = 1; i <= l; ++i, --vm, ++vp)
      *vm = std::conj*vp = C(std::cos(i*x), std::sin(i*x)));
  #else  // recursive calculations of expI(ix)
    C *u = v + l, e = C(std::cos(x), std::sin(x));
    *u = one;
    for (int i = 1; i <= l; ++i) u[-i] = conj(u[i] = e*u[i-1]);

  #endif
  }
}

/*
  Calculate Wigner's D-matrix

    D = [D_{m1,m2} ^j(a,b,c)]_{m1,m2=-j}^j

  at j.

  Input:

    a - angle of rotation around z axis
    b - angle of rotation around y axis
    c - angle of rotation around z axis

    j - azimuthal quantum number

  Output:

    D[2j+1][2j+1] - Wigner's D matrix
*/

template <class T>
void wigner_D(const T &a, const T &b,  const T& c,
              int l, std::complex<T> **D) {

  using C = std::complex<T>;

  int dim = 2*l + 1;

  //
  // calculate small d matrix due to Jy rotation
  //
  T **d = matrix<T> (dim, dim);

  wigner_d(b, l,  d);

  //
  // calculate phasors due to Jz rotation
  //
  C *va = new C[dim],
    *vc = new C[dim];

  calc_expI(a, l, va);
  calc_expI(c, l, vc);

  //
  // calculate D matrix
  //
  for (int m = 0; m < dim; ++m)
    for (int n = 0; n < dim; ++n)
      D[m][n] = va[m]*d[m][n]*vc[n];

  //
  // cleanup
  //
  free_matrix(d);
  delete [] va;
  delete [] vc;
}
