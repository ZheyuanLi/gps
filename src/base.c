#include <Rinternals.h>
#include "base.h"
SEXP C_SetAttr (SEXP x, SEXP Name, SEXP Value) {
  setAttrib(x, install(CHAR(asChar(Name))), Value);
  return x;
}
SEXP C_SetDim (SEXP x, SEXP Value) {
  setAttrib(x, R_DimSymbol, Value);
  return x;
}
SEXP C_VecDot (SEXP x, SEXP y) {
  int n = length(x);
  if (length(y) != n) error("length(x) == length(y) expected!");
  double c = DOT(n, REAL(x), REAL(y));
  return ScalarReal(c);
}
SEXP C_VecScal (SEXP alpha, SEXP x, SEXP overwrite) {
  int n = length(x);
  SEXP y = x; double *ptrx = REAL(x), *ptry = ptrx;
  int MakeCopy = 1 - asInteger(overwrite);
  if (MakeCopy) {
    y = PROTECT(allocVector(REALSXP, n));
    ptry = REAL(y); VecCopy(n, ptrx, ptry);
  }
  SCAL(n, asReal(alpha), ptry);
  if (MakeCopy) UNPROTECT(1);
  return y;
}
