static inline double DOT (int n, double *x, double *y) {
  double c = 0.0, *xi = x, *yi = y, *xn = x + n;
  if (x == y) {
    while (xi < xn) {c += (*xi) * (*xi); xi++;}
  } else {
    while (xi < xn) c += (*xi++) * (*yi++);
  }
  return c;
}
static inline void AXPY (int n, double alpha, double *x, double *y) {
  double *xi = x, *xn = x + n, *yi = y;
  if (alpha == 1.0) {while (xi < xn) (*yi++) += (*xi++);}
  else if (alpha == -1.0) {while (xi < xn) (*yi++) -= (*xi++);}
  else {while (xi < xn) (*yi++) += alpha * (*xi++);}
}
static inline void SCAL (int n, double alpha, double *x) {
  double *xi = x, *xn = x + n;
  if (alpha != 1.0) while (xi < xn) (*xi++) *= alpha;
}
static inline void SCALnew (int n, double alpha, double *x, double *y) {
  double *xi = x, *xn = x + n, *yi = y;
  while (xi < xn) *yi++ = alpha * (*xi++);
}
static inline double L2Loss (int n, double *x, double *y) {
  double *xi = x, *yi = y, *xn = x + n, c = 0.0, a;
  while (xi < xn) {a = (*xi++); a -= (*yi++); c += a * a;}
  return c;
}
static inline void ZeroVec (int n, double *x) {
  double *xi = x, *xn = x + n;
  if (n > 0) xi[0] = 0;
  for (xi += n & 1; xi < xn; xi += 2) {
    xi[0] = 0; xi[1] = 0;
  }
}
static inline void ZeroAntiLowerTri (int n, double *A, int LDA) {
  double *Aij, *Anj = A + n, *AD = A + (n - 1), *Ann = A + LDA * n;
  while (Anj <= Ann) {
    Aij = AD + 1;
    while (Aij < Anj) *Aij++ = 0.0;
    AD += LDA - 1; Anj += LDA;
  }
}
static inline void VecCopy (int n, double *x1, double *x2) {
  int i;
  if (n > 0) x2[0] = x1[0];
  for (i = n & 1; i < n; i += 2) {
    x2[i] = x1[i]; x2[i + 1] = x1[i + 1];
  }
}
static inline void MatCopy (int m, int n, double *A, int LDA, double *B, int LDB) {
  double *A0j = A, *B0j = B, *A0n = A + n * LDA;
  for (; A0j < A0n; A0j += LDA, B0j += LDB) VecCopy(m, A0j, B0j);
}
static inline void MatTrans (int m, int n, double *A, int LDA, double *B, int LDB) {
  double *A0j = A, *Amj = A0j + m, *B0j = B, *B0n = B + n, *Aij, *Bij;
  while (B0j < B0n) {
    Aij = A0j; Bij = B0j;
    while (Aij < Amj) {*Bij = *Aij; Bij += LDB; Aij++;}
    A0j += LDA; Amj += LDA; B0j++;
  }
}
static inline void Mat2Sym (int u2l, int n, double *A, int LDA) {
  double *Ajj = A, *Anj = A + n, *Ann = A + LDA * n, *Aij, *Aji;
  while (Ajj < Ann) {
    Aij = Ajj + 1; Aji = Ajj + LDA;
    while (Aij < Anj) {
      if (u2l) *Aij = *Aji; else *Aji = *Aij;
      Aij++; Aji += LDA;
    }
    Ajj += LDA + 1; Anj += LDA;
  }
}
