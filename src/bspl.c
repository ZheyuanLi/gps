#include <stdlib.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "base.h"
int IsMonoInc (int n, double *x) {
  int flag = 1; double *xi, *xn = x + n - 1;
  for (xi = x; xi < xn; xi++) {
    if (xi[1] <= xi[0]) {flag = 0; break;}
  }
  return flag;
}
SEXP C_IsMonoInc (SEXP x, SEXP n, SEXP xi) {
  if (!isReal(x)) error("'x' is not in double-precision mode!");
  int i = asInteger(xi), l = length(x);
  if (i < 1 || i > l) error("'xi' is out of bound!");
  double *subx = REAL(x) + i - 1;
  int N = asInteger(n);
  if (N > l - i + 1) error("n <= length(x) - xi + 1 required!");
  int flag = IsMonoInc(N, subx);
  return ScalarInteger(flag);
}
void MakeGrid (double *b, int k, int n, double *x, int rmdup) {
  int ni = k - 1;
  double *ptrb = b, *bend = b + ni, b0, b1;
  double step0 = 1.0 / (n - 1), step;
  int n0 = n - rmdup;
  if (rmdup) x[0] = b[0];
  double *ptrx = x + rmdup, *pend = ptrx + (n0 - 1);
  for (; ptrb < bend; ptrb++, pend += n0) {
    b0 = ptrb[0]; b1 = ptrb[1];
    step = (b1 - b0) * step0;
    if (rmdup) b0 += step;
    for (; ptrx < pend; b0 += step, ptrx++) *ptrx = b0;
    if (b0 < b1) *ptrx++ = b0;
    else {
      if (b1 > 0) step = b1;
      else if (b1 < 0) step = -b1;
      else step = 1.0;
      *ptrx++ = b1 - 2.22e-16 * step;
    }
  }
}
SEXP C_MakeGrid (SEXP b, SEXP n, SEXP rmdup) {
  int K = length(b), N = asInteger(n), RMDUP = asInteger(rmdup);
  SEXP x = PROTECT(allocVector(REALSXP, (K - 1) * (N - RMDUP) + RMDUP));
  MakeGrid(REAL(b), K, N, REAL(x), RMDUP);
  UNPROTECT(1);
  return x;
}
void SmallAtA (int n, double alpha, double *A, double *X) {
  double *A0i, *A0j = A, *A0n = A + n * n, *Aki, *Akj, *Anj;
  double *Xjj = X, *Xij, *Xji;
  double c;
  while (A0j < A0n) {
    Anj = A0j + n;
    Akj = A0j; c = 0.0;
    while (Akj < Anj) {
      c += Akj[0] * Akj[0];
      Akj++;
    }
    c *= alpha; *Xjj = c;
    A0i = A0j + n; Xij = Xjj + 1; Xji = Xjj + n;
    while (A0i < A0n) {
      Aki = A0i; Akj = A0j; c = 0.0;
      while (Akj < Anj) {
        c += Aki[0] * Akj[0];
        Aki++; Akj++;
      }
      c *= alpha; *Xij = c; *Xji = c;
      A0i += n; Xij++; Xji += n;
    }
    A0j += n; Xjj += n + 1;
  }
}
void SmallLtA (int n, double *L, double *A, double *X) {
  double *Xij = X, *Xnn = X + n * n;
  double *Lii, *Lki;
  double *Aij = A, *Akj, *Anj = A;
  double c;
  while (Xij < Xnn) {
    Lii = L; Anj += n;
    while (Aij < Anj) {
      Lki = Lii; Akj = Aij; c = 0.0;
      while (Akj < Anj) {
        c += Lki[0] * Akj[0];
        Lki++; Akj++;
      }
      *Xij++ = c; Lii += n + 1; Aij++;
    }
  }
}
SEXP C_SbarBlocks (SEXP xd, SEXP W, SEXP B) {
  int ord = nrows(W);
  int k1 = length(xd) - 1;
  double *s = REAL(xd);
  double *b = s + k1;
  double *L = REAL(W); int blocksize;
  F77_CALL(dpotf2)("l", &ord, L, &ord, &blocksize FCONE);
  double *Bj = REAL(B);
  blocksize = ord * ord;
  double alpha;
  double *X = malloc(blocksize * sizeof(double));
  SEXP S = PROTECT(alloc3DArray(REALSXP, ord, ord, k1));
  double *Sj = REAL(S);
  while (s < b) {
    SmallLtA(ord, L, Bj, X);
    alpha = 0.5 * (s[1] - s[0]);
    SmallAtA(ord, alpha, X, Sj);
    Bj += blocksize; Sj += blocksize; s++;
  }
  free(X);
  UNPROTECT(1);
  return S;
}
static inline void Block2LTB (int n, double *A, double *L) {
  double *Ajj = A, *Aij, *Anj = A, *Ann = A + n * n;
  double *L0j = L, *Lij;
  while (Ajj < Ann) {
    Aij = Ajj; Lij = L0j; Anj += n;
    while (Aij < Anj) *Lij++ += *Aij++;
    Ajj += n + 1; L0j += n;
  }
}
SEXP C_SbarLTB (SEXP S) {
  SEXP Dim = getAttrib(S, R_DimSymbol);
  int *dim = INTEGER(Dim);
  int ord = dim[0];
  int k1 = dim[2];
  int n = k1 + ord - 1;
  SEXP LTB = PROTECT(allocMatrix(REALSXP, ord, n));
  double *L = REAL(LTB), *Lj = L; int blocksize = length(LTB);
  ZeroVec(blocksize, L);
  blocksize = ord * ord;
  double *Sj = REAL(S), *Send = Sj + blocksize * k1;
  while (Sj < Send) {
    Block2LTB(ord, Sj, Lj);
    Sj += blocksize; Lj += ord;
  }
  UNPROTECT(1);
  return LTB;
}
double xtAx (int n, double *A, double *x) {
  double c = 0.0, alpha;
  double *xj = x, *xi, *xn = x + n;
  double *Ajj = A, *Aij;
  while (xj < xn) {
    alpha = xj[0];
    c += Ajj[0] * alpha * alpha;
    Aij = Ajj + 1; xi = xj + 1; alpha += alpha;
    while (xi < xn) {
      c += Aij[0] * xi[0] * alpha;
      Aij++; xi++;
    }
    Ajj += n + 1; xj++;
  }
  return c;
}
SEXP C_btSb (SEXP S, SEXP b) {
  SEXP Dim = getAttrib(S, R_DimSymbol);
  int *dim = INTEGER(Dim);
  int ord = dim[0];
  int k1 = dim[2];
  int n = k1 + ord - 1;
  if (length(b) != n) error("Incorrect number of coefficients!");
  SEXP c = PROTECT(allocVector(REALSXP, k1));
  double *cj = REAL(c), *cend = cj + k1;
  int blocksize = ord * ord;
  double *Sj = REAL(S);
  double *bj = REAL(b);
  while (cj < cend) {
    *cj = xtAx(ord, Sj, bj);
    Sj += blocksize; cj++; bj++;
  }
  UNPROTECT(1);
  return c;
}
