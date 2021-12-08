/*
  Testing Wigner D matrices.

  Author: Martin Horvat, November 2021
*/

#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <complex>
#include <limits>
#include <string>

#include <matrix.h>
#include <wignerD.h>
#include <print.h>

int main(){

  using myreal = double;

  print_title(std::cout, "Wigner d - matrix, j=2, a=1, res1.txt");
  {
    int
      j = 2,
      dim = 2*j + 1;

    myreal
      a = 1;

    auto d = matrix <myreal> (dim, dim);

    wigner_d(a, j, d);

    print_mat(std::cout, dim, dim, d);

    free_matrix(d);
  }


  print_title(std::cout, "Wigner d - matrix, j=3, a=1, res2.txt");
  {
    int
      j = 3,
      dim = 2*j + 1;

    myreal
      a = 1;

    auto d = matrix <myreal> (dim, dim);

    wigner_d(a, j, d);

    print_mat(std::cout, dim, dim, d);

    free_matrix(d);
  }


  print_title(std::cout, "Wigner D - matrix, j=2, (a,b,c) = (1,2,3), res3.txt");
  {
    int
      j = 2,
      dim = 2*j + 1;

    myreal
      a = 1, b = 2, c = 3;

    auto **D = matrix <std::complex<myreal>>(dim, dim);

    wigner_D(a, b, c, j, D);

    print_mat(std::cout, dim, dim, D);

    free_matrix(D);
  }


  print_title(std::cout, "Wigner D - matrix, j=10, (a,b,c) = (1,2,3), res4.txt");
  {
    int
      j = 10,
      dim = 2*j + 1;

    myreal
      a = 1, b = 2, c = 3;

    auto **D = matrix <std::complex<myreal>>(dim, dim);

    wigner_D(a, b, c, j, D);

    print_Cmat(std::cout, dim, dim, D);

    free_matrix(D);
  }

  return 0;
}
