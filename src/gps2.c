#include <stdlib.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "base.h"
SEXP C_CheckSupp (SEXP nx, SEXP ord) {
  int k = length(nx) - 1, *x = INTEGER(nx), d = asInteger(ord);
  int flag = 0, *xk = x + k, *xl, *xr, *xi, n;
  if (x[0] == 0 || xk[0] == 0) flag = 1;
  if (k > d) {
    for (xl = x + 1, xr = xl + d; xr < xk; xl++, xr++) {
      n = 0; for (xi = xl; xi < xr; xi++) n += xi[0];
      if (n == 0) {flag = 1; break;}
    }
  }
  return ScalarInteger(flag);
}
void LTB2Dense (int b1, int n, int k, double *L, double *A);
void Dx (int n, int b1, double *Dt, double *x, double *y);
void Dtx (int n, int b1, double *Dt, double *x, double *y);
void LTB2Dense (int b1, int n, int k, double *L, double *A);
void FormE (int d, int p, double *L, int b1D, int q, double *Dt, double *E) {
  int bwL = d - 1, ione = 1, n = p;
  double *Ejj = E, *Epq = E + p * q, *L0j = L;
  LTB2Dense(b1D, q, p, Dt, E);
  while (Ejj < Epq) {
    F77_CALL(dtbsv)("l", "n", "n", &n, &bwL, L0j, &d, Ejj, &ione);
    L0j += d; Ejj += p + 1; n--;
  }
}
SEXP C_FormE (SEXP L, SEXP Dt) {
  int d = nrows(L), p = ncols(L), b1D = nrows(Dt), q = ncols(Dt);
  SEXP E = PROTECT(allocMatrix(REALSXP, p, q));
  FormE(d, p, REAL(L), b1D, q, REAL(Dt), REAL(E));
  UNPROTECT(1);
  return E;
}
double MeanEigen (int p, int q, double *E) {
  double w = 0.0, *Ejj = E, *Eij, *Epj = E + p, *Epq = E + p * q;
  while (Ejj < Epq) {
    for (Eij = Ejj; Eij < Epj; Eij++) w += Eij[0] * Eij[0];
    Ejj += p + 1; Epj += p;
  }
  w /= (double)q;
  return w;
}
SEXP C_MeanEigen (SEXP E) {
  double w = MeanEigen(nrows(E), ncols(E), REAL(E));
  return ScalarReal(w);
}
void WoodburyA (int p, int q, double *E, double *F, double *G) {
  int m = p - q, ione = 1;
  double done = 1.0, dzero = 0.0, *E2 = E + q, *ptr, *end;
  if (m > 1) {
    MatTrans(m, q, E2, p, F, q);
    F77_CALL(dtrsm)("l", "l", "t", "n", &q, &m, &done, E, &p, F, &q);
    F77_CALL(dsyrk)("l", "t", &m, &q, &done, F, &q, &dzero, G, &m);
    end = G + m * m; for (ptr = G; ptr < end; ptr += m + 1) *ptr += 1.0;
    F77_CALL(dpotf2)("l", &m, G, &m, &ione);
    F77_CALL(dtrsm)("l", "l", "n", "n", &q, &m, &done, E, &p, F, &q);
  } else {
    end = F + q; for (ptr = F; ptr < end; ptr++, E2 += p) *ptr = *E2;
    F77_CALL(dtrsv)("l", "t", "n", &q, E, &p, F, &ione);
    done = DOT(q, F, F); done += 1.0; *G = sqrt(done);
    F77_CALL(dtrsv)("l", "n", "n", &q, E, &p, F, &ione);
  }
}
void SolveA (int p, int q, double *E, double *F, double *G, double *v, double *u) {
  int m = p - q, ione = 1;
  double done = 1.0, dzero = 0.0, *w1 = u + q, *w2 = w1 + q;
  VecCopy(q, v, u);
  F77_CALL(dtrsv)("l", "t", "n", &q, E, &p, u, &ione);
  F77_CALL(dtrsv)("l", "n", "n", &q, E, &p, u, &ione);
  if (m > 1) {
    F77_CALL(dgemv)("t", &q, &m, &done, F, &q, v, &ione, &dzero, w2, &ione);
    F77_CALL(dtrsv)("l", "n", "n", &m, G, &m, w2, &ione);
    F77_CALL(dtrsv)("l", "t", "n", &m, G, &m, w2, &ione);
    F77_CALL(dgemv)("n", &q, &m, &done, F, &q, w2, &ione, &dzero, w1, &ione);
    AXPY(q, -1.0, w1, u);
  } else {
    dzero = G[0] * G[0];
    done = -1.0 * DOT(q, F, v) / dzero;
    AXPY(q, done, F, u);
  }
}
int MinEigen (int p, int q, double *E, double *w, double tol) {
  int m = p - q, k, i; double norm, w0 = 0.0, w1 = 0.0;
  int sizeF = q * m, sizeG = m * m, sizeu = q;
  if (m > 1) sizeu += q + m;
  double *work = malloc((sizeF + sizeG + q + sizeu) * sizeof(double));
  double *F = work, *G = F + sizeF, *v = G + sizeG, *u = v + q;
  WoodburyA(p, q, E, F, G);
  while (v < u) *v++ = 1.0; v -= q;
  SolveA(p, q, E, F, G, v, u);
  k = 0;
  while (1) {
    norm = DOT(q, u, u);
    norm = sqrt(norm);
    SCALnew(q, 1.0 / norm, u, v);
    SolveA(p, q, E, F, G, v, u);
    w1 = DOT(q, u, v);
    if (fabs(w1 - w0) < w0 * tol) break;
    w0 = w1; k++;
  }
  free(work);
  *w = 1.0 / w1;
  return k;
}
SEXP C_MinEigen (SEXP E, SEXP tol) {
  int p = nrows(E), q = ncols(E), k;
  SEXP w = PROTECT(allocVector(REALSXP, 1));
  k = MinEigen(p, q, REAL(E), REAL(w), asReal(tol));
  setAttrib(w, R_NamesSymbol, ScalarInteger(k));
  UNPROTECT(1);
  return w;
}
int MaxEigen (int d, int p, double *L, int b1D, int q, double *Dt,
              double *w, double tol) {
  double norm, w0 = 0.0, w1 = 0.0;
  int m = p - q, bwL = d - 1, ione = 1, k;
  int sizeu = p; if (b1D == d) sizeu += bwL - m;
  double *work = malloc((q + sizeu) * sizeof(double));
  double *v = work, *u = v + q;
  while (v < u) *v++ = 1.0; v = work;
  Dtx(q, b1D, Dt, v, u);
  F77_CALL(dtbsv)("l", "n", "n", &p, &bwL, L, &d, u, &ione);
  k = 0;
  while (1) {
    F77_CALL(dtbsv)("l", "t", "n", &p, &bwL, L, &d, u, &ione);
    Dx(q, b1D, Dt, u, v);
    norm = DOT(q, v, v);
    norm = sqrt(norm);
    SCAL(q, 1.0 / norm, v);
    Dtx(q, b1D, Dt, v, u);
    F77_CALL(dtbsv)("l", "n", "n", &p, &bwL, L, &d, u, &ione);
    w1 = DOT(p, u, u);
    if (fabs(w1 - w0) < w0 * tol) break;
    w0 = w1; k++;
  }
  free(work);
  *w = w1;
  return k;
}
SEXP C_MaxEigen (SEXP L, SEXP Dt, SEXP tol) {
  int d = nrows(L), p = ncols(L), b1D = nrows(Dt), q = ncols(Dt), k;
  SEXP w = PROTECT(allocVector(REALSXP, 1));
  k = MaxEigen(d, p, REAL(L), b1D, q, REAL(Dt), REAL(w), asReal(tol));
  setAttrib(w, R_NamesSymbol, ScalarInteger(k));
  UNPROTECT(1);
  return w;
}
int ApproxEigen (int q, double *w, double min, double max, double mean, double tol) {
  double *work = malloc(3 * q * sizeof(double));
  double *B2 = work, *alpha1 = B2 + q, *alpha2 = alpha1 + q, *end = alpha1;
  int k; double beta0, beta1, g = (double)q, dg, Newton, thresh;
  double logMin = log(min), logMax = log(max), sum = g * mean;
  double h = 1.0 / (g + 1.0), p = h, xmax = log(g), r = 0.5 / xmax, x01, x10;
  while (B2 < end) {
    x01 = (log(1.0 / p - 1.0) + xmax) * r;
    x10 = 1.0 - x01;
    g = 2 * x01 * x10;
    dg = exp(x10 * x10 * logMin + x01 * x01 * logMax);
    *B2++ = g; *alpha1++ = dg; *alpha2++ = g * dg;
    p += h;
  }
  B2 = work; alpha1 = end; alpha2 = end + q;
  g = -sum;
  while (B2 < end) {
    g += alpha1[0] * exp(B2[0] * beta0);
    B2++; alpha1++;
  }
  B2 = work; alpha1 = end;
  beta0 = 0.0;
  k = 0;
  while (1) {
#ifdef trace
    printf("k = %d, beta = %.2f, g = %.2e\n", k, beta0, g);
#endif
    dg = 0.0;
    while (B2 < end) {
      dg += alpha2[0] * exp(B2[0] * beta0);
      B2++; alpha2++;
    }
    B2 = work; alpha2 = alpha1 + q;
    Newton = g / dg;
    if (Newton > 5) {
      Newton = 5;
    } else if (Newton < -5) {
      Newton = -5;
    }
    if (fabs(Newton) < fabs(beta0) * tol) break;
    thresh = fabs(g);
    while (1) {
#ifdef trace
      printf(" -> g' = %.2e, Newton = %.2e\n", dg, Newton);
#endif
      beta1 = beta0 - Newton;
      g = -sum;
      while (B2 < end) {
        g += alpha1[0] * exp(B2[0] * beta1);
        B2++; alpha1++;
      }
      B2 = work; alpha1 = end;
      if (fabs(g) < thresh) break;
      Newton *= 0.5;
    }
    beta0 = beta1; k++;
  }
  while (B2 < end) {
    *w++ = alpha1[0] * exp(B2[0] * beta0);
    alpha1++; B2++;
    }
  free(work);
  return k;
}
SEXP C_ApproxEigen (SEXP q, SEXP min, SEXP max, SEXP mean, SEXP tol) {
  int Q = asInteger(q);
  SEXP w = PROTECT(allocVector(REALSXP, Q));
  ApproxEigen(Q, REAL(w), asReal(min), asReal(max), asReal(mean), asReal(tol));
  UNPROTECT(1);
  return w;
}
SEXP C_LAUUM (SEXP E) {
  int p = nrows(E), q = ncols(E), m;
  double *E1 = REAL(E), *E2 = E1 + q;
  SEXP A = PROTECT(allocMatrix(REALSXP, q, q));
  double done = 1.0, *ptrA = REAL(A);
  MatCopy(q, q, E1, p, ptrA, q);
  F77_CALL(dlauum)("l", &q, ptrA, &q, &m);
  m = p - q;
  F77_CALL(dsyrk)("l", "t", &q, &m, &done, E2, &p, &done, ptrA, &q);
  Mat2Sym(0, q, ptrA, q);
  UNPROTECT(1);
  return A;
}
void Rho2EDF (int q, double *w, int n, double *rho, double *edf) {
  double *wi, *wq = w + q;
  double *rhoj, *end = rho + n;
  double *edfj = edf, df, lambda;
  for (rhoj = rho; rhoj < end; rhoj++) {
    lambda = exp(rhoj[0]); df = 0.0;
    for (wi = w; wi < wq; wi++) {
      df += 1.0 / (1.0 + lambda * wi[0]);
    }
    *edfj++ = df;
  }
}
SEXP C_Rho2EDF (SEXP w, SEXP rho) {
  int q = length(w), n = length(rho);
  SEXP edf = PROTECT(allocVector(REALSXP, n));
  Rho2EDF(q, REAL(w), n, REAL(rho), REAL(edf));
  UNPROTECT(1);
  return edf;
}
int EDF2Rho_unstable (int q, double *w, double edf, double *init, double tol) {
  double *wi, *wq = w + q; int k;
  double rho = init[0], lambda, x, g, dg, Newton;
  if (edf <= 0 || edf >= q) error("'edf' is out of bound!");
  for (k = 0; k < 100; k++) {
    lambda = exp(rho); g = -edf; dg = 0.0;
    for (wi = w; wi < wq; wi++) {
      x = 1.0 + lambda * wi[0];
      g += 1.0 / x;
      dg += lambda * wi[0] / (x * x);
    }
    Newton = g / dg;
#ifdef trace
    printf("k = %d, rho = %.2f, g = %.2e, g' = %.2e, Newton = %.2e\n",
           k, rho, g, dg, Newton);
#endif
    if (fabs(Newton) < fabs(rho) * tol) break;
    rho += Newton;
  }
  *init = rho;
  return k;
}
int EDF2Rho (int q, double *w, double edf, double *init,
             double MaxNewton, double tol) {
  double *wi, *wq = w + q; int k;
  double rho0 = init[0], lambda, rho1;
  double alpha, g, dg, Newton, thresh;
  if (edf <= 0 || edf >= q) error("'edf' is out of bound!");
  lambda = exp(rho0);
  for (g = -edf, wi = w; wi < wq; wi++) {
    alpha = 1.0 + lambda * wi[0];
    g += 1.0 / alpha;
  }
  k = 0;
  while (1) {
#ifdef trace
    printf("k = %d, rho = %.2f, g = %.2e\n", k, rho0, g);
#endif
    for (dg = 0.0, wi = w; wi < wq; wi++) {
      alpha = 1.0 + lambda * wi[0];
      dg += lambda * wi[0] / (alpha * alpha);
    }
    Newton = g / dg;
    if (Newton > MaxNewton) {
      Newton = MaxNewton;
    } else if (Newton < -MaxNewton) {
      Newton = -MaxNewton;
    }
    if (fabs(Newton) < fabs(rho0) * tol) break;
    thresh = fabs(g);
    while (1) {
#ifdef trace
      printf(" -> g' = %.2e, Newton = %.2e\n", dg, Newton);
#endif
      rho1 = rho0 + Newton;
      lambda = exp(rho1);
      for (g = -edf, wi = w; wi < wq; wi++) {
        alpha = 1.0 + lambda * wi[0];
        g += 1.0 / alpha;
      }
      if (fabs(g) < thresh) break;
      Newton *= 0.5;
    }
    rho0 = rho1; k++;
  }
  *init = rho0;
  return k;
}
SEXP C_EDF2Rho (SEXP w, SEXP edf, SEXP init, SEXP MaxNewton, SEXP tol) {
  int q = length(w), k;
  SEXP rho = PROTECT(allocVector(REALSXP, 1));
  double *RHO = REAL(rho); RHO[0] = asReal(init);
  k = EDF2Rho(q, REAL(w), asReal(edf), RHO, asReal(MaxNewton), asReal(tol));
  setAttrib(rho, R_NamesSymbol, ScalarInteger(k));
  UNPROTECT(1);
  return rho;
}
void FormK (int p, double *Z, int d, double *S, int b1S, double lambda, double *K) {
  int bwZ = d - 1, info; double *Sij = S, *Kij = K;
  double *Zij = Z, *Zpd = Z + p * d, *end1 = Z + b1S, *end2 = Z + d;
  while (Zij < Zpd) {
    while (Zij < end1) (*Kij++) = (*Zij++) + lambda * (*Sij++);
    while (Zij < end2) (*Kij++) = (*Zij++);
    end1 += d; end2 += d;
  }
  F77_CALL(dpbtrf)("l", &p, &bwZ, K, &d, &info);
  if (info) error("Normal matrix is rank-deficient!");
}
void SolvePLS (int d, int p, double *K, double *z, int r, double *Ct, double *beta,
               double *V, double *G, double *E, double *u, double *beta_c) {
  int bwZ = d - 1, ione = 1, info;
  double dzero = 0.0, done = 1.0, *V0j, *V0r, *Ct0j, tmp;
  VecCopy(p, z, beta);
  F77_CALL(dtbsv)("l", "n", "n", &p, &bwZ, K, &d, beta, &ione);
  V0j = V; V0r = V + p * r; Ct0j = Ct;
  while (V0j < V0r) {
    VecCopy(p, Ct0j, V0j);
    F77_CALL(dtbsv)("l", "n", "n", &p, &bwZ, K, &d, V0j, &ione);
    V0j += p; Ct0j += p;
  }
  if (r > 1) {
    F77_CALL(dsyrk)("l", "t", &r, &p, &done, V, &p, &dzero, G, &r);
    F77_CALL(dpotrf)("l", &r, G, &r, &info);
    if (info) error("Constraints are not linearly independent!");
    MatTrans(p, r, V, p, E, r);
    F77_CALL(dtrsm)("l", "l", "n", "n", &r, &p, &done, G, &r, E, &r);
    F77_CALL(dgemv)("n", &r, &p, &done, E, &r, beta, &ione, &dzero, u, &ione);
    F77_CALL(dgemv)("t", &r, &p, &done, E, &r, u, &ione, &dzero, beta_c, &ione);
    F77_CALL(dtbsv)("l", "t", "n", &p, &bwZ, K, &d, beta_c, &ione);
  } else if (r == 1) {
    tmp = DOT(p, V, V); tmp = sqrt(tmp); *G = tmp;
    SCALnew(p, 1.0 / tmp, V, E);
    tmp = DOT(p, E, beta);
    SCALnew(p, tmp, E, beta_c);
    F77_CALL(dtbsv)("l", "t", "n", &p, &bwZ, K, &d, beta_c, &ione);
  }
  F77_CALL(dtbsv)("l", "t", "n", &p, &bwZ, K, &d, beta, &ione);
  if (r) AXPY(p, -1.0, beta_c, beta);
}
double PLS2EDF (int d, int p, double *K, double *L, int r, double *E,
                double *B0, double *B0_c) {
  int bwK = d - 1, ione = 1, pi;
  double edf = 0.0, *Kii, *B0ii, done = 1.0;
  LTB2Dense(d, p, p, L, B0);
  pi = p; Kii = K; B0ii = B0;
  while (pi) {
    F77_CALL(dtbsv)("l", "n", "n", &pi, &bwK, Kii, &d, B0ii, &ione);
    edf += DOT(pi, B0ii, B0ii);
    Kii += d; B0ii += p + 1; pi--;
  }
  if (r > 1) {
    MatCopy(r, p, E, r, B0_c, r);
    F77_CALL(dtrmm)("r", "l", "n", "n", &r, &p, &done, B0, &p, B0_c, &r);
    edf -= DOT(r * p, B0_c, B0_c);
  } else if (r == 1) {
    VecCopy(p, E, B0_c);
    F77_CALL(dtrmv)("l", "t", "n", &p, B0, &p, B0_c, &ione);
    edf -= DOT(p, B0_c, B0_c);
  }
  return edf;
}
SEXP C_GridPLS (SEXP Z, SEXP L, SEXP S, SEXP z, SEXP Ct, SEXP rho) {
  int d = nrows(Z), p = ncols(Z), b1S = nrows(S), r = ncols(Ct), l = length(rho);
  double *ptrZ = REAL(Z), *ptrL = REAL(L), *ptrS = REAL(S), *ptrz = REAL(z);
  double *ptrCt = REAL(Ct), *rhoj = REAL(rho), *rhol = rhoj + l, lambda;
  SEXP beta = PROTECT(allocMatrix(REALSXP, p, l));
  SEXP edf = PROTECT(allocVector(REALSXP, l));
  double *betaj = REAL(beta), *edfj = REAL(edf);
  int sizeK = d * p, sizeV = p * r, sizeG = r * r, sizeB0 = p * p;
  int sizework = sizeK + 3 * sizeV + sizeG + r + p + sizeB0;
  double *work = malloc(sizework * sizeof(double));
  double *K = work, *V = K + sizeK, *G = V + sizeV, *E = G + sizeG;
  double *u = E + sizeV, *beta_c = u + r, *B0 = beta_c + p, *B0_c = B0 + sizeB0;
  while (rhoj < rhol) {
    lambda = exp(rhoj[0]);
    FormK(p, ptrZ, d, ptrS, b1S, lambda, K);
    SolvePLS(d, p, K, ptrz, r, ptrCt, betaj, V, G, E, u, beta_c);
    edfj[0] = PLS2EDF(d, p, K, ptrL, r, E, B0, B0_c);
    rhoj++; betaj += p; edfj++;
  }
  free(work);
  SEXP LIST = PROTECT(allocVector(VECSXP, 3));
  SEXP NAMES = PROTECT(allocVector(STRSXP, 3));
  SET_VECTOR_ELT(LIST, 0, rho);
  SET_STRING_ELT(NAMES, 0, mkChar("rho"));
  SET_VECTOR_ELT(LIST, 1, beta);
  SET_STRING_ELT(NAMES, 1, mkChar("beta"));
  SET_VECTOR_ELT(LIST, 2, edf);
  SET_STRING_ELT(NAMES, 2, mkChar("edf"));
  setAttrib(LIST, R_NamesSymbol, NAMES);
  UNPROTECT(4);
  return LIST;
}
double GCVscore (int n, int d, int p, double *L, double *f, double minRSS,
                 double *beta, double edf, double *fHat) {
  int bwL = d - 1, ione = 1; double N = (double)n, ResiEDF, RSS, gcv;
  VecCopy(p, beta, fHat);
  F77_CALL(dtbmv)("l", "t", "n", &p, &bwL, L, &d, fHat, &ione);
  ResiEDF = N - edf;
  RSS = L2Loss(p, f, fHat) + minRSS;
  gcv = (N * RSS) / (ResiEDF * ResiEDF);
  return gcv;
}
SEXP C_GridGCV (SEXP n, SEXP L, SEXP f, SEXP minRSS, SEXP beta, SEXP edf) {
  int d = nrows(L), p = ncols(L), N = asInteger(n), l = length(edf);
  double *ptrL = REAL(L), *ptrf = REAL(f), RSSmin = asReal(minRSS);
  double *betaj = REAL(beta), *edfj = REAL(edf), *edfl = edfj + l;
  SEXP gcv = PROTECT(allocVector(REALSXP, l));
  double *gcvj = REAL(gcv);
  double *fHat = malloc(p * sizeof(double));
  while (edfj < edfl) {
    *gcvj = GCVscore(N, d, p, ptrL, ptrf, RSSmin, betaj, *edfj, fHat);
    betaj += p; edfj++; gcvj++;
  }
  free(fHat);
  UNPROTECT(1);
  return gcv;
}
