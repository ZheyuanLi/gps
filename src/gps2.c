#include <stdlib.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "base.h"
#define MaxIterEigen 1000
#define MaxIterNewton 200
#include <R_ext/Random.h>
void RandomVec (int n, double *x) {
  double *xn = x + n, *xi = x;
  GetRNGstate();
  while (xi < xn) *xi++ = unif_rand();
  PutRNGstate();
}
SEXP C_CheckSupp (SEXP nx, SEXP ord) {
  int k = length(nx) - 1, *x = INTEGER(nx), d = asInteger(ord);
  int flag = 0, *xk = x + k, *xl, *xr, *xi, n;
  if (*x == 0 || *xk == 0) flag = 1;
  if (k > d) {
    for (xl = x + 1, xr = xl + d; xr < xk; xl++, xr++) {
      n = 0; for (xi = xl; xi < xr; xi++) n += *xi;
      if (n == 0) {flag = 1; break;}
    }
  }
  return ScalarInteger(flag);
}
void LTB2Dense (int b1, int n, int k, double *L, double *A);
void Dx (int n, int b1, double *Dt, double *x, double *y);
void Dtx (int n, int b1, double *Dt, double *x, double *y);
void FormE (int d, int p, double *L, int b1D, int q, double *Dt, double *E) {
  int bwL = d - 1, ione = 1, n = p;
  double *Ejj = E, *Epq = E + p * q, *L0j = L;
  LTB2Dense(b1D, q, p, Dt, E);
  while (Ejj < Epq) {
    F77_CALL(dtbsv)("l", "n", "n", &n, &bwL, L0j, &d, Ejj, &ione FCONE FCONE FCONE);
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
double MeanDR (int p, int q, double *E) {
  double w = 0.0, *Ejj = E, *Eij, *Epj = E + p, *Epq = E + p * q;
  while (Ejj < Epq) {
    for (Eij = Ejj; Eij < Epj; Eij++) w += Eij[0] * Eij[0];
    Ejj += p + 1; Epj += p;
  }
  w /= (double)q;
  return w;
}
SEXP C_MeanDR (SEXP E) {
  double w = MeanDR(nrows(E), ncols(E), REAL(E));
  return ScalarReal(w);
}
void WoodburyA (int p, int q, double *E, double *F, double *G) {
  int m = p - q, ione = 1;
  double done = 1.0, dzero = 0.0, *E2 = E + q, *ptr, *end;
  if (m > 1) {
    MatTrans(m, q, E2, p, F, q);
    F77_CALL(dtrsm)("l", "l", "t", "n", &q, &m, &done, E, &p, F, &q
                      FCONE FCONE FCONE FCONE);
    F77_CALL(dsyrk)("l", "t", &m, &q, &done, F, &q, &dzero, G, &m FCONE FCONE);
    end = G + m * m; for (ptr = G; ptr < end; ptr += m + 1) *ptr += 1.0;
    F77_CALL(dpotf2)("l", &m, G, &m, &ione FCONE);
    F77_CALL(dtrsm)("l", "l", "n", "n", &q, &m, &done, E, &p, F, &q
                      FCONE FCONE FCONE FCONE);
  } else {
    end = F + q; for (ptr = F; ptr < end; ptr++, E2 += p) *ptr = *E2;
    F77_CALL(dtrsv)("l", "t", "n", &q, E, &p, F, &ione FCONE FCONE FCONE);
    done = DOT(q, F, F); done += 1.0; *G = sqrt(done);
    F77_CALL(dtrsv)("l", "n", "n", &q, E, &p, F, &ione FCONE FCONE FCONE);
  }
}
void SolveA (int p, int q, double *E, double *F, double *G, double *v, double *u) {
  int m = p - q, ione = 1;
  double done = 1.0, dzero = 0.0, *w1 = u + q, *w2 = w1 + q;
  VecCopy(q, v, u);
  F77_CALL(dtrsv)("l", "t", "n", &q, E, &p, u, &ione FCONE FCONE FCONE);
  F77_CALL(dtrsv)("l", "n", "n", &q, E, &p, u, &ione FCONE FCONE FCONE);
  if (m > 1) {
    F77_CALL(dgemv)("t", &q, &m, &done, F, &q, v, &ione, &dzero, w2, &ione FCONE);
    F77_CALL(dtrsv)("l", "n", "n", &m, G, &m, w2, &ione FCONE FCONE FCONE);
    F77_CALL(dtrsv)("l", "t", "n", &m, G, &m, w2, &ione FCONE FCONE FCONE);
    F77_CALL(dgemv)("n", &q, &m, &done, F, &q, w2, &ione, &dzero, w1, &ione FCONE);
    AXPY(q, -1.0, w1, u);
  } else {
    dzero = G[0] * G[0];
    done = -1.0 * DOT(q, F, v) / dzero;
    AXPY(q, done, F, u);
  }
}
int MinDR (int p, int q, double *E, double *w, double tol) {
  int m = p - q, k; double norm, w0 = 0.0, w1 = 0.0;
  int sizeF = q * m, sizeG = m * m, sizeu = q;
  if (m > 1) sizeu += q + m;
  double *work = malloc((sizeF + sizeG + q + sizeu) * sizeof(double));
  double *F = work, *G = F + sizeF, *v = G + sizeG, *u = v + q;
  WoodburyA(p, q, E, F, G);
  RandomVec(q, v);
  SolveA(p, q, E, F, G, v, u);
#ifdef trace
  printf("Finding the smallest Demmler-Reinsch eigenvalue:\n");
#endif
  k = 0;
  while (k < MaxIterEigen) {
    norm = DOT(q, u, u);
    norm = sqrt(norm);
    SCALnew(q, 1.0 / norm, u, v);
    SolveA(p, q, E, F, G, v, u);
    w1 = DOT(q, u, v);
#ifdef trace
    printf("-> k = %d, eigenvalue = %.4e\n", k, w1);
#endif
    if (w1 < 0) break;
    if (fabs(w1 - w0) < w0 * tol) break;
    w0 = w1; k++;
  }
#ifdef trace
  putchar('\n');
#endif
  free(work);
  *w = 1.0 / w1;
  return k;
}
SEXP C_MinDR (SEXP E, SEXP tol) {
  int p = nrows(E), q = ncols(E), k;
  SEXP w = PROTECT(allocVector(REALSXP, 1));
  double *ptrw = REAL(w);
  k = MinDR(p, q, REAL(E), ptrw, asReal(tol));
  if (k == MaxIterEigen) {
    warning("Unable to find the smallest eigenvalue in %d iterations!", MaxIterEigen);
    *ptrw = 0.0;
  }
  UNPROTECT(1);
  return w;
}
int MaxDR (int d, int p, double *L, int b1D, int q, double *Dt,
           double *w, double tol) {
  double norm, w0 = 0.0, w1 = 0.0;
  int m = p - q, bwL = d - 1, ione = 1, k;
  int sizeu = p; if (b1D == d) sizeu += bwL - m;
  double *work = malloc((q + sizeu) * sizeof(double));
  double *v = work, *u = v + q;
  RandomVec(q, v);
  Dtx(q, b1D, Dt, v, u);
  F77_CALL(dtbsv)("l", "n", "n", &p, &bwL, L, &d, u, &ione FCONE FCONE FCONE);
#ifdef trace
  printf("Finding the largest Demmler-Reinsch eigenvalue:\n");
#endif
  k = 0;
  while (k < MaxIterEigen) {
    F77_CALL(dtbsv)("l", "t", "n", &p, &bwL, L, &d, u, &ione FCONE FCONE FCONE);
    Dx(q, b1D, Dt, u, v);
    norm = DOT(q, v, v);
    norm = sqrt(norm);
    SCAL(q, 1.0 / norm, v);
    Dtx(q, b1D, Dt, v, u);
    F77_CALL(dtbsv)("l", "n", "n", &p, &bwL, L, &d, u, &ione FCONE FCONE FCONE);
    w1 = DOT(p, u, u);
#ifdef trace
    printf("-> k = %d, eigenvalue = %.4e\n", k, w1);
#endif
    if (fabs(w1 - w0) < w0 * tol) break;
    w0 = w1; k++;
  }
#ifdef trace
  putchar('\n');
#endif
  free(work);
  *w = w1;
  return k;
}
SEXP C_MaxDR (SEXP L, SEXP Dt, SEXP tol) {
  int d = nrows(L), p = ncols(L), b1D = nrows(Dt), q = ncols(Dt), k;
  SEXP w = PROTECT(allocVector(REALSXP, 1));
  k = MaxDR(d, p, REAL(L), b1D, q, REAL(Dt), REAL(w), asReal(tol));
  if (k == MaxIterEigen) {
    warning("Unable to find the largest eigenvalue in %d iterations!", MaxIterEigen);
  }
  UNPROTECT(1);
  return w;
}
void Q1ApproxDR (int q, double min, double max, double gamma,
                 double *theta, double *h, double *eta,
                 double *alpha_l, double *alpha_r) {
  double a = log(min), b = log(max);
  double Q = (double)q, log_q = log(Q), log_q1 = log(Q + 1.0);
  double t_j = (gamma - 1.0) * log_q1;
  double t_1 = log_q + t_j, t_q = t_j - gamma * log_q, r = 1.0 / (t_1 - t_q);
  double p_j = 1.0 / (Q + 1.0), step = p_j, z_j, dtmp1, dtmp2;
  double *theta_j = theta, *h_j = h, *h_q = h + q, *eta_j = eta;
  while (h_j < h_q) {
    t_j = log(1.0 - p_j) - gamma * log(p_j);
    z_j = (t_j - t_q) * r;
    dtmp1 = z_j * z_j - z_j; h_j[0] = dtmp1;
    dtmp2 = exp(a + (b - a) * z_j); theta_j[0] = dtmp2;
    eta_j[0] = dtmp1 * dtmp2;
    p_j += step; h_j++; theta_j++; eta_j++;
  }
  *alpha_l = 0; *alpha_r = b - a;
}
void Q2ApproxDR (int q, double min, double max, double gamma,
                 double *theta, double *h, double *eta,
                 double *alpha_l, double *alpha_r) {
  double a = log(min), b = log(max);
  double Q = (double)q, log_q = log(Q), log_q1 = log(Q + 1.0);
  double t_j = (gamma - 1.0) * log_q1;
  double t_1 = log_q + t_j, t_q = t_j - gamma * log_q, r = 1.0 / (t_1 - t_q);
  double p_j = 1.0 / (Q + 1.0), step = p_j, z_j, y_j, B0, B1, B2, B3, dtmp1, dtmp2;
  double *theta_j = theta, *h_j = h, *h_q = h + q, *eta_j = eta;
  while (h_j < h_q) {
    t_j = log(1.0 - p_j) - gamma * log(p_j);
    z_j = (t_j - t_q) * r;
    y_j = 1.0 - z_j;
    B0 = y_j * y_j * y_j;
    B1 = 3.0 * z_j * y_j * y_j;
    B2 = 3.0 * z_j * z_j * y_j;
    B3 = z_j * z_j * z_j;
    dtmp1 = B1 - B2; h_j[0] = dtmp1;
    dtmp2 = exp((B0 + B2) * a + (B2 + B3) * b); theta_j[0] = dtmp2;
    eta_j[0] = dtmp1 * dtmp2;
    p_j += step; h_j++; theta_j++; eta_j++;
  }
  *alpha_l = a; *alpha_r = (2.0 * a + b) / 3.0;
}
int RootApproxDR (int q, double *w, double *theta, double *h, double *eta,
                  double mean, double alpha_l, double alpha_r, double tol) {
  double *theta_j, *eta_j, *h_j, *h_q = h + q;
  int k; double alpha0, alpha1, g, dg, MaxNewton, Newton, thresh;
  double sum = q * mean, g_l = -sum, g_r = -sum;
  theta_j = theta; h_j = h;
  while (h_j < h_q) {
    g_l += theta_j[0] * exp(h_j[0] * alpha_l);
    g_r += theta_j[0] * exp(h_j[0] * alpha_r);
    theta_j++; h_j++;
  }
  if (g_l * g_r > 0.0) {
#ifdef trace
    printf("No 'alpha' is found!\n\n");
#endif
    return 0;
  }
  alpha0 = 0.5 * (alpha_l + alpha_r);
  g = -sum; theta_j = theta; h_j = h;
  while (h_j < h_q) {
    g += theta_j[0] * exp(h_j[0] * alpha0);
    theta_j++; h_j++;
  }
  MaxNewton = 0.25 * (alpha_r - alpha_l);
  k = 0;
  while (k < MaxIterNewton) {
#ifdef trace
    printf("-> k = %d, alpha = %.4f, g = %.4e\n", k, alpha0, g);
#endif
    dg = 0.0; eta_j = eta; h_j = h;
    while (h_j < h_q) {
      dg += eta_j[0] * exp(h_j[0] * alpha0);
      eta_j++; h_j++;
    }
    Newton = g / dg;
    if (Newton > MaxNewton) {
      Newton = MaxNewton;
    } else if (Newton < -MaxNewton) {
      Newton = -MaxNewton;
    }
    if (fabs(Newton) < fabs(alpha0) * tol) break;
    thresh = fabs(g);
    while (1) {
#ifdef trace
      printf("--> g' = %.4e, Newton = %.4e\n", dg, Newton);
#endif
      alpha1 = alpha0 - Newton;
      g = -sum; theta_j = theta; h_j = h;
      while (h_j < h_q) {
        g += theta_j[0] * exp(h_j[0] * alpha1);
        theta_j++; h_j++;
      }
      if (fabs(g) < thresh) break;
      Newton *= 0.5;
    }
    alpha0 = alpha1; k++;
  }
#ifdef trace
  putchar('\n');
#endif
  theta_j = theta; h_j = h; eta_j = w;
  while (h_j < h_q) {
    eta_j[0] += theta_j[0] * exp(h_j[0] * alpha0);
    theta_j++; h_j++; eta_j++;
  }
  return k;
}
SEXP C_RootApproxDR (SEXP theta, SEXP h, SEXP eta, SEXP mean, SEXP interval) {
  int k, q = length(theta); double *bounds = REAL(interval);
  SEXP w = PROTECT(allocVector(REALSXP, q));
  double *ptrw = REAL(w); ZeroVec(q, ptrw);
  k = RootApproxDR(q, ptrw, REAL(theta), REAL(h), REAL(eta),
                   asReal(mean), bounds[0], bounds[1], 1e-6);
  UNPROTECT(1);
  return w;
}
int ApproxDR (int q, double *w, double min, double max, double mean, double tol) {
  double *work = malloc(3 * q * sizeof(double));
  double *h = work, *theta = h + q, *eta = theta + q;
  ZeroVec(q, w);
  double alpha_l, alpha_r;
  double gamma = 0.0; int k, success1, success2, success;
  for (success1 = 0; gamma < 1.01; gamma += 0.05) {
#ifdef trace
    printf("Trying z(%.2f) with convex Q1(z, alpha).\n", gamma);
#endif
    Q1ApproxDR(q, min, max, gamma, theta, h, eta, &alpha_l, &alpha_r);
    k = RootApproxDR(q, w, theta, h, eta, mean, alpha_l, alpha_r, tol);
    if (k) success1++; else if (success1) break;
  }
  for (success2 = 0; gamma < 1.01; gamma += 0.05) {
#ifdef trace
    printf("Trying z(%.2f) with 'S'-shaped Q2(z, alpha).\n", gamma);
#endif
    Q2ApproxDR(q, min, max, gamma, theta, h, eta, &alpha_l, &alpha_r);
    k = RootApproxDR(q, w, theta, h, eta, mean, alpha_l, alpha_r, tol);
    if (k) success2++; else if (success2) break;
  }
  success = success1 + success2;
  free(work);
  return success;
}
SEXP C_ApproxDR (SEXP q, SEXP min, SEXP max, SEXP mean, SEXP tol) {
  int Q = asInteger(q), N;
  SEXP w = PROTECT(allocVector(REALSXP, Q));
  double *ptrw = REAL(w), *w_q = ptrw + Q, r;
  N = ApproxDR(Q, REAL(w), asReal(min), asReal(max), asReal(mean), asReal(tol));
  if (N) {
    r = 1.0 / (double)N;
    while (ptrw < w_q) *ptrw++ *= r;
  } else {
    while (ptrw < w_q) *ptrw++ = NA_REAL;
    warning("Unable to approximate Demmler-Reinsch eigenvalues!");
  }
  UNPROTECT(1);
  return w;
}
SEXP C_LAUUM (SEXP E) {
  int p = nrows(E), q = ncols(E), m;
  double *E1 = REAL(E), *E2 = E1 + q;
  SEXP A = PROTECT(allocMatrix(REALSXP, q, q));
  double done = 1.0, *ptrA = REAL(A);
  MatCopy(q, q, E1, p, ptrA, q);
  F77_CALL(dlauum)("l", &q, ptrA, &q, &m FCONE);
  m = p - q;
  F77_CALL(dsyrk)("l", "t", &q, &m, &done, E2, &p, &done, ptrA, &q FCONE FCONE);
  Mat2Sym(0, q, ptrA, q);
  UNPROTECT(1);
  return A;
}
void Rho2REDF (int q, double *w, int n, double *rho, double *redf) {
  double *wi, *wq = w + q;
  double *rhoj, *end = rho + n;
  double *redfj = redf, dtmp, lambda;
  for (rhoj = rho; rhoj < end; rhoj++) {
    lambda = exp(*rhoj); dtmp = 0.0;
    for (wi = w; wi < wq; wi++) {
      dtmp += 1.0 / (1.0 + lambda * wi[0]);
    }
    *redfj++ = dtmp;
  }
}
SEXP C_Rho2REDF (SEXP w, SEXP rho) {
  int q = length(w), n = length(rho);
  SEXP redf = PROTECT(allocVector(REALSXP, n));
  Rho2REDF(q, REAL(w), n, REAL(rho), REAL(redf));
  UNPROTECT(1);
  return redf;
}
int REDF2Rho (int q, double *w, double redf, double *rho,
              double MaxNewton, double tol) {
  double *wi, *wq = w + q; int k;
  double rho0 = *rho, lambda, rho1;
  double alpha, g, dg, Newton, thresh;
  if (redf <= 0 || redf >= q) error("'redf' is out of bound!");
  lambda = exp(rho0);
  for (g = -redf, wi = w; wi < wq; wi++) {
    alpha = 1.0 + lambda * wi[0];
    g += 1.0 / alpha;
  }
#ifdef trace
  printf("Finding 'rho' that matches the target 'redf':\n");
#endif
  k = 0;
  while (k < MaxIterNewton) {
#ifdef trace
    printf("-> k = %d, rho = %.4f, g = %.4e\n", k, rho0, g);
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
      printf("--> g' = %.4e, Newton = %.4e\n", dg, Newton);
#endif
      rho1 = rho0 + Newton;
      lambda = exp(rho1);
      for (g = -redf, wi = w; wi < wq; wi++) {
        alpha = 1.0 + lambda * wi[0];
        g += 1.0 / alpha;
      }
      if (fabs(g) < thresh) break;
      Newton *= 0.5;
    }
    rho0 = rho1; k++;
  }
#ifdef trace
  putchar('\n');
#endif
  *rho = rho0;
  return k;
}
SEXP C_REDF2Rho (SEXP w, SEXP redf, SEXP rho0, SEXP MaxNewton, SEXP tol) {
  int q = length(w), k;
  SEXP rho = PROTECT(allocVector(REALSXP, 1));
  double *RHO = REAL(rho); *RHO = asReal(rho0);
  k = REDF2Rho(q, REAL(w), asReal(redf), RHO, asReal(MaxNewton), asReal(tol));
  if (k == MaxIterNewton) {
    warning("Unable to find 'rho' to match 'redf' in %d iterations!", MaxIterNewton);
  }
  UNPROTECT(1);
  return rho;
}
int FormK (int p, double *Z, int d, double *S, int b1S, double lambda, double *K) {
  int bwZ = d - 1, info; double *Sij = S, *Kij = K;
  double *Zij = Z, *Zpd = Z + p * d, *end1 = Z + b1S, *end2 = Z + d;
  while (Zij < Zpd) {
    while (Zij < end1) (*Kij++) = (*Zij++) + lambda * (*Sij++);
    while (Zij < end2) (*Kij++) = (*Zij++);
    end1 += d; end2 += d;
  }
  F77_CALL(dpbtrf)("l", &p, &bwZ, K, &d, &info FCONE);
  return info;
}
void SolvePWLS (int d, int p, double *K, double *z, int r, double *Ct, double *beta,
                double *V, double *G, double *E, double *u, double *beta_c) {
  int bwZ = d - 1, ione = 1, info;
  double dzero = 0.0, done = 1.0, *V0j, *V0r, *Ct0j, tmp;
  VecCopy(p, z, beta);
  F77_CALL(dtbsv)("l", "n", "n", &p, &bwZ, K, &d, beta, &ione FCONE FCONE FCONE);
  V0j = V; V0r = V + p * r; Ct0j = Ct;
  while (V0j < V0r) {
    VecCopy(p, Ct0j, V0j);
    F77_CALL(dtbsv)("l", "n", "n", &p, &bwZ, K, &d, V0j, &ione FCONE FCONE FCONE);
    V0j += p; Ct0j += p;
  }
  if (r > 1) {
    F77_CALL(dsyrk)("l", "t", &r, &p, &done, V, &p, &dzero, G, &r FCONE FCONE);
    F77_CALL(dpotrf)("l", &r, G, &r, &info FCONE);
    if (info) error("Constraints are not linearly independent!");
    MatTrans(p, r, V, p, E, r);
    F77_CALL(dtrsm)("l", "l", "n", "n", &r, &p, &done, G, &r, E, &r
                      FCONE FCONE FCONE FCONE);
    F77_CALL(dgemv)("n", &r, &p, &done, E, &r, beta, &ione, &dzero, u, &ione FCONE);
    F77_CALL(dgemv)("t", &r, &p, &done, E, &r, u, &ione, &dzero, beta_c, &ione FCONE);
    F77_CALL(dtbsv)("l", "t", "n", &p, &bwZ, K, &d, beta_c, &ione FCONE FCONE FCONE);
  } else if (r == 1) {
    tmp = DOT(p, V, V); tmp = sqrt(tmp); *G = tmp;
    SCALnew(p, 1.0 / tmp, V, E);
    tmp = DOT(p, E, beta);
    SCALnew(p, tmp, E, beta_c);
    F77_CALL(dtbsv)("l", "t", "n", &p, &bwZ, K, &d, beta_c, &ione FCONE FCONE FCONE);
  }
  F77_CALL(dtbsv)("l", "t", "n", &p, &bwZ, K, &d, beta, &ione FCONE FCONE FCONE);
  if (r) AXPY(p, -1.0, beta_c, beta);
}
double PWLS2EDF (int d, int p, double *K, double *L, int r, double *E,
                 double *B0, double *B0_c) {
  int bwK = d - 1, ione = 1, pi;
  double edf = 0.0, *Kii, *B0ii, done = 1.0;
  LTB2Dense(d, p, p, L, B0);
  pi = p; Kii = K; B0ii = B0;
  while (pi) {
    F77_CALL(dtbsv)("l", "n", "n", &pi, &bwK, Kii, &d, B0ii, &ione FCONE FCONE FCONE);
    edf += DOT(pi, B0ii, B0ii);
    Kii += d; B0ii += p + 1; pi--;
  }
  if (r > 1) {
    MatCopy(r, p, E, r, B0_c, r);
    F77_CALL(dtrmm)("r", "l", "n", "n", &r, &p, &done, B0, &p, B0_c, &r
                      FCONE FCONE FCONE FCONE);
    edf -= DOT(r * p, B0_c, B0_c);
  } else if (r == 1) {
    VecCopy(p, E, B0_c);
    F77_CALL(dtrmv)("l", "t", "n", &p, B0, &p, B0_c, &ione FCONE FCONE FCONE);
    edf -= DOT(p, B0_c, B0_c);
  }
  return edf;
}
SEXP C_GridPWLS (SEXP Z, SEXP L, SEXP S, SEXP z, SEXP Ct, SEXP rho) {
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
  int info = 0;
  while (rhoj < rhol) {
    lambda = exp(*rhoj);
    info = FormK(p, ptrZ, d, ptrS, b1S, lambda, K);
    if (info) {
      for (info = 0; info < p; info++) betaj[info] = NA_REAL;
      *edfj = NA_REAL;
    } else {
      SolvePWLS(d, p, K, ptrz, r, ptrCt, betaj, V, G, E, u, beta_c);
      *edfj = PWLS2EDF(d, p, K, ptrL, r, E, B0, B0_c);
    }
    rhoj++; betaj += p; edfj++;
  }
  if (info) warning("Penalized least squares is not solvable for some 'rho' values!");
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
  F77_CALL(dtbmv)("l", "t", "n", &p, &bwL, L, &d, fHat, &ione FCONE FCONE FCONE);
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
    if (ISNA(*edfj)) {
      *gcvj = NA_REAL;
    } else {
      *gcvj = GCVscore(N, d, p, ptrL, ptrf, RSSmin, betaj, *edfj, fHat);
    }
    betaj += p; edfj++; gcvj++;
  }
  free(fHat);
  UNPROTECT(1);
  return gcv;
}
