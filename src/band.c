#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "base.h"
void Csparse2LTB (int b1, int n, int k, double *x, double *L) {
  int b = b1 - 1, r = k - n;
  double *start, *end1, *end2, *end3, *ptrL, *ptrx = x;
  start = L; end1 = L + b; end2 = L;
  while (start < end1) {
    ptrL = start;
    while (ptrL <= end2) {
      *ptrL = *ptrx++; ptrL += b;
    }
    start++; end2 += b1;
  }
  start = end1; end1 += (n - b) * b1;
  while (start < end1) {
    ptrL = start;
    while (ptrL <= end2) {
      *ptrL = *ptrx++; ptrL += b;
    }
    start += b1; end2 += b1;
  }
  end3 = end2 + r;
  while (end2 < end3) {
    ptrL = start;
    while (ptrL < end2) {
      *ptrL = *ptrx++; ptrL += b;
    }
    start += b1; end2++;
  }
  end3 = end2 + (b - r);
  while (end2 < end3) {
    ptrL = start;
    while (ptrL < end2) {
      *ptrL = 0.0; ptrL += b;
    }
    start += b1; end2++;
  }
}
SEXP C_Csparse2LTB (SEXP b1, SEXP n, SEXP k, SEXP x) {
  int B1 = asInteger(b1), N = asInteger(n), K = asInteger(k);
  SEXP L = PROTECT(allocMatrix(REALSXP, B1, N));
  Csparse2LTB(B1, N, K, REAL(x), REAL(L));
  UNPROTECT(1);
  return L;
}
void LTB2Dense (int b1, int n, int k, double *L, double *A) {
  double *Lcut = L + b1 * (k - b1);
  double *Ajj = A, *Aij = A, *Akj, *Akn;
  double *Lij = L, *Lbj = Lij + b1;
  while (Lij < Lcut) {
    while (Aij < Ajj) *Aij++ = 0.0;
    while (Lij < Lbj) *Aij++ = *Lij++;
    Ajj += k + 1; Lbj += b1;
  }
  Akj = Ajj + b1; Akn = A + k * n;
  while (Akj <= Akn) {
    while (Aij < Ajj) *Aij++ = 0.0;
    Lij = Lcut; while (Aij < Akj) *Aij++ = *Lij++;
    Ajj += k + 1; Akj += k; Lcut += b1;
  }
}
SEXP C_LTB2Dense (SEXP L, SEXP k) {
  int B1 = nrows(L), N = ncols(L), K = asInteger(k);
  if (K < N || K >= N + B1) error("'k' is out of bound!");
  SEXP A = PROTECT(allocMatrix(REALSXP, K, N));
  LTB2Dense(B1, N, K, REAL(L), REAL(A));
  UNPROTECT(1);
  return A;
}
SEXP C_SolveLTB (SEXP transA, SEXP A, SEXP y, SEXP overwrite) {
  int ione = 1, n = ncols(A), b1 = nrows(A), bw = b1 - 1, k;
  char nt = 'n'; if (asInteger(transA)) nt = 't';
  double *ptrA = REAL(A);
  if (isMatrix(y)) {
    if (nrows(y) != n) error("nrow(y) == ncol(A) expected!");
    k = ncols(y);
  } else {
    if (length(y) != n) error("length(y) == ncol(A) expected!");
    k = 1;
  }
  SEXP x = y; double *ptry = REAL(y), *ptrx = ptry;
  int MakeCopy = 1 - asInteger(overwrite);
  if (MakeCopy) {
    x = PROTECT(allocVector(REALSXP, n * k));
    ptrx = REAL(x); VecCopy(n * k, ptry, ptrx);
    if (k > 1) setAttrib(x, R_DimSymbol, getAttrib(y, R_DimSymbol));
  }
  double *xend = ptrx + n * k;
  while (ptrx < xend) {
    F77_CALL(dtbsv)("l", &nt, "n", &n, &bw, ptrA, &b1, ptrx, &ione FCONE FCONE FCONE);
    ptrx += n;
  }
  if (MakeCopy) UNPROTECT(1);
  return x;
}
SEXP C_LPBTRF (SEXP A, SEXP overwrite) {
  int n = ncols(A), b1 = nrows(A), bw = b1 - 1;
  SEXP X = A; double *ptrA = REAL(A), *ptrX = ptrA;
  int MakeCopy = 1 - asInteger(overwrite);
  if (MakeCopy) {
    X = PROTECT(allocMatrix(REALSXP, b1, n));
    ptrX = REAL(X); VecCopy(b1 * n, ptrA, ptrX);
  }
  int info; F77_CALL(dpbtrf)("l", &n, &bw, ptrX, &b1, &info FCONE);
  ZeroAntiLowerTri(b1, ptrX + (n - b1) * b1, b1);
  if (MakeCopy) UNPROTECT(1);
  if (info) error("The leading minor of order %d is not positive definite!", info);
  return X;
}
void Dx (int n, int b1, double *Dt, double *x, double *y) {
  double *Dt0j = Dt, *xj = x, *yj = y, *yn = y + n, *xi, *di, a;
  while (yj < yn) {
    di = Dt0j; xi = xj; Dt0j += b1;
    a = 0.0; while (di < Dt0j) a += (*di++) * (*xi++);
    *yj++ = a; xj++;
  }
}
void DX (int n, int b1, int k, double *Dt, double *X, int LDX, double *Y, int LDY) {
  double *X0j = X, *Y0j = Y, *X0k = X + k * LDX;
  while (X0j < X0k) {
    Dx(n, b1, Dt, X0j, Y0j);
    X0j += LDX; Y0j += LDY;
  }
}
void Dtx (int n, int b1, double *Dt, double *x, double *y) {
  double *Dt0j = Dt, *xj = x, *xn = x + n, *yj = y, *yi, *di, a;
  ZeroVec(n + b1 - 1, y);
  while (xj < xn) {
    di = Dt0j; yi = yj; Dt0j += b1;
    a = *xj; while (di < Dt0j) *yi++ += a * (*di++);
    yj++; xj++;
  }
}
void DtX (int n, int b1, int k, double *Dt, double *X, int LDX, double *Y, int LDY) {
  double *X0j = X, *Y0j = Y, *X0k = X + k * LDX;
  while (X0j < X0k) {
    Dtx(n, b1, Dt, X0j, Y0j);
    X0j += LDX; Y0j += LDY;
  }
}
void DDt (int n, int b1, double *Dt, double *X) {
  double *dj = Dt, *di, *db, *ptrX = X, *bar, *x, *y, a;
  bar = Dt + (n - b1) * b1;
  while (dj < bar) {
    db = dj + b1;
    a = 0.0; y = dj; while (y < db) {a += (*y) * (*y); y++;}
    *ptrX++ = a; di = db; dj++;
    while (dj < db) {
      a = 0.0; x = di; y = dj; while (y < db) a += (*x++) * (*y++);
      *ptrX++ = a; di += b1; dj++;
    }
    dj = db;
  }
  bar += b1 * b1;
  while (dj < bar) {
    db = dj + b1;
    a = 0.0; y = dj; while (y < db) {a += (*y) * (*y); y++;}
    *ptrX++ = a; di = db; dj++;
    while (di < bar) {
      a = 0.0; x = di; y = dj; while (y < db) a += (*x++) * (*y++);
      *ptrX++ = a; di += b1; dj++;
    }
    while (dj < db) {*ptrX++ = 0.0; dj++;}
    dj = db;
  }
}
SEXP C_DDt (SEXP Dt) {
  int b1 = nrows(Dt), n = ncols(Dt);
  SEXP X = PROTECT(allocMatrix(REALSXP, b1, n));
  DDt(n, b1, REAL(Dt), REAL(X));
  UNPROTECT(1);
  return X;
}
void DtD (int n, int b1, double *Dt, double *X) {
  double *Dt0j = Dt, *Dt0n = Dt + b1 * n, *X0j = X, *Y0j, *Yij, a, *di, *dj;
  int m = n + b1 - 1;
  ZeroVec(b1 * m, X);
  while (Dt0j < Dt0n) {
    dj = Dt0j; Dt0j += b1; Y0j = X0j;
    while (dj < Dt0j) {
      a = *dj; Yij = Y0j;
      di = dj; while (di < Dt0j) *Yij++ += a * (*di++);
      Y0j += b1; dj++;
    }
    X0j += b1;
  }
}
SEXP C_DtD (SEXP Dt) {
  int b1 = nrows(Dt), n = ncols(Dt);
  SEXP X = PROTECT(allocMatrix(REALSXP, b1, n + b1 - 1));
  DtD(n, b1, REAL(Dt), REAL(X));
  UNPROTECT(1);
  return X;
}
