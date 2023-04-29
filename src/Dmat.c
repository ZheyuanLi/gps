#include <Rinternals.h>
#include "base.h"
void Diff (int n, int k, double *x, double *dx) {
  double *xi = x, *yi = x + k, *xn = x + n, *dxi = dx, alpha, tmp;
  if (k == 1) {
    while (yi < xn) {
      tmp = yi[0] - yi[-1];
      *dxi++ = tmp;
      yi++;
    }
  } else {
    alpha = 1.0 / k;
    while (yi < xn) {
      tmp = (*yi++) - (*xi++);
      tmp *= alpha;
      *dxi++ = tmp;
    }
  }
}
SEXP C_Diff (SEXP x, SEXP k, SEXP n, SEXP xi) {
  if (!isReal(x)) error("'x' is not in double-precision mode!");
  int i = asInteger(xi), l = length(x);
  if (i < 1 || i > l) error("'xi' is out of bound!");
  double *subx = REAL(x) + i - 1;
  int N = asInteger(n), K = asInteger(k);
  if (N > l - i + 1) error("n <= length(x) - xi + 1 required!");
  if (K <= 0 || K >= N) error("1 <= k <= n - 1 required!");
  SEXP dx = PROTECT(allocVector(REALSXP, N - K));
  Diff(N, K, subx, REAL(dx));
  UNPROTECT(1);
  return dx;
}
void ComputeLD (double *xt, int k, int d, double *ld) {
  int m = d - 1, p = k - d, i;
  double *dx1, *dx2;
  for (i = 1; i <= m; i++) {
    dx1 = ld + (i - 1) * p; dx2 = dx1 + i;
    while (dx1 < dx2) *dx1++ = 0.0;
    Diff(k - 2 * i, d - i, xt + i, dx1);
  }
}
SEXP C_ComputeLD (SEXP xt, SEXP d) {
  if (!isReal(xt)) error("'xt' is not in double-precision mode!");
  int K = length(xt), D = asInteger(d);
  SEXP ld = PROTECT(allocMatrix(REALSXP, K - D, D - 1));
  ComputeLD(REAL(xt), K, D, REAL(ld));
  UNPROTECT(1);
  return ld;
}
void NullVec (double *ld, int p, int m, double *h) {
  double *dx, *hp = h + p, *hi, c; int j, skip;
  skip = (m - 1); ZeroVec(skip, h);
  hi = h + skip; while (hi < hp) *hi++ = 1.0;
  j = m - 1;
  while (j--) {
    dx = ld + j * p + skip; hi = h + skip; c = 0.0;
    while (hi < hp) {
      c += (*dx) * (*hi);
      *hi = c; hi++; dx++;
    }
  }
  hi = h + skip; c = 0.0;
  while (hi < hp) {c += (*hi) * (*hi); hi++;}
  hi = h + skip; c = 1.0 / sqrt(c);
  while (hi < hp) *hi++ *= c;
}
void NullGD (double *ld, int p, int m, double *H) {
  int i; double *h = H;
  for (i = 1; i <= m; i++, h += p) NullVec(ld, p, i, h);
}
SEXP C_NullGD (SEXP ld, SEXP m) {
  int P = nrows(ld), M = asInteger(m);
  SEXP H = PROTECT(allocMatrix(REALSXP, P, M));
  NullGD(REAL(ld), P, M, REAL(H));
  UNPROTECT(1);
  return H;
}
