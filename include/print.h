#pragma once

/*
  Auxillary routines for printing of results

  Author: Martin Horvat, December 2021
*/

#include <iostream>
#include <cmath>
#include <complex>
#include <string>

template <typename T>
std::string W(const std::complex<T> & x) {

  int len = std::numeric_limits<T>::max_digits10;
  char buf[256], format[256];

  sprintf(format, "(%%%d.%de %%%d.%de)", len+8,len, len+8,len);
  sprintf(buf, format, x.real(), x.imag());

  return std::string(buf);
}

template <typename T>
std::string W(const T & x) {

  int len = std::numeric_limits<T>::max_digits10;
  char buf[256], format[256];

  sprintf(format, "%%%d.%de", len+8,len);
  sprintf(buf, format, x);

  return std::string(buf);
}

template <class T, int n, int m>
void print_mat(std::ostream & o, T A[n][m]){
  int m1 = m - 1;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      o << W(A[i][j]) << (j < m1 ? ' ' : '\n');
}

template <class T>
void print_mat(std::ostream & o, int n, int m, T **A){
  int m1 = m - 1;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      o << W(A[i][j]) << (j < m1 ? ' ' : '\n');
}

template <class T>
void print_Cmat(std::ostream & o, int n, int m, std::complex<T> **A){
  int m1 = m - 1;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      o << W(A[i][j].real()) << ' ' << W(A[i][j].imag()) << (j < m1 ? ' ' : '\n');
}

void print_title(std::ostream & o, std::string title){
  o << std::string(10, '=') <<  " " << title << " "<< std::string(70-2 - title.size(),'=') << '\n';
}
