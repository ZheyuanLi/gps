#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
extern SEXP C_ApproxEigen(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_btSb(SEXP, SEXP);
extern SEXP C_CheckSupp(SEXP, SEXP);
extern SEXP C_ComputeLD(SEXP, SEXP);
extern SEXP C_Csparse2LTB(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_DDt(SEXP);
extern SEXP C_Diff(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_DtD(SEXP);
extern SEXP C_EDF2Rho(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_FormE(SEXP, SEXP);
extern SEXP C_GridGCV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_GridPLS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_IsAscending(SEXP, SEXP, SEXP);
extern SEXP C_LAUUM(SEXP);
extern SEXP C_LPBTRF(SEXP, SEXP);
extern SEXP C_MakeGrid(SEXP, SEXP, SEXP);
extern SEXP C_MaxEigen(SEXP, SEXP, SEXP);
extern SEXP C_MeanEigen(SEXP);
extern SEXP C_MinEigen(SEXP, SEXP);
extern SEXP C_NullD(SEXP, SEXP);
extern SEXP C_Rho2EDF(SEXP, SEXP);
extern SEXP C_SbarBlocks(SEXP, SEXP, SEXP);
extern SEXP C_SbarLTB(SEXP, SEXP);
extern SEXP C_SetAttr(SEXP, SEXP, SEXP);
extern SEXP C_SetDim(SEXP, SEXP);
extern SEXP C_SolveLTB(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_VecDot(SEXP, SEXP);
extern SEXP C_VecScal(SEXP, SEXP, SEXP);
static const R_CallMethodDef CallEntries[] = {
    {"C_ApproxEigen", (DL_FUNC) &C_ApproxEigen, 5},
    {"C_btSb", (DL_FUNC) &C_btSb, 2},
    {"C_CheckSupp", (DL_FUNC) &C_CheckSupp, 2},
    {"C_ComputeLD", (DL_FUNC) &C_ComputeLD, 2},
    {"C_Csparse2LTB", (DL_FUNC) &C_Csparse2LTB, 4},
    {"C_DDt", (DL_FUNC) &C_DDt, 1},
    {"C_Diff", (DL_FUNC) &C_Diff, 4},
    {"C_DtD", (DL_FUNC) &C_DtD, 1},
    {"C_EDF2Rho", (DL_FUNC) &C_EDF2Rho, 5},
    {"C_FormE", (DL_FUNC) &C_FormE, 2},
    {"C_GridGCV", (DL_FUNC) &C_GridGCV, 6},
    {"C_GridPLS", (DL_FUNC) &C_GridPLS, 6},
    {"C_IsAscending", (DL_FUNC) &C_IsAscending, 3},
    {"C_LAUUM", (DL_FUNC) &C_LAUUM, 1},
    {"C_LPBTRF", (DL_FUNC) &C_LPBTRF, 2},
    {"C_MakeGrid", (DL_FUNC) &C_MakeGrid, 3},
    {"C_MaxEigen", (DL_FUNC) &C_MaxEigen, 3},
    {"C_MeanEigen", (DL_FUNC) &C_MeanEigen, 1},
    {"C_MinEigen", (DL_FUNC) &C_MinEigen, 2},
    {"C_NullD", (DL_FUNC) &C_NullD, 2},
    {"C_Rho2EDF", (DL_FUNC) &C_Rho2EDF, 2},
    {"C_SbarBlocks", (DL_FUNC) &C_SbarBlocks, 3},
    {"C_SbarLTB", (DL_FUNC) &C_SbarLTB, 2},
    {"C_SetAttr", (DL_FUNC) &C_SetAttr, 3},
    {"C_SetDim", (DL_FUNC) &C_SetDim, 2},
    {"C_SolveLTB", (DL_FUNC) &C_SolveLTB, 4},
    {"C_VecDot", (DL_FUNC) &C_VecDot, 2},
    {"C_VecScal", (DL_FUNC) &C_VecScal, 3},
    {NULL, NULL, 0}
};
void R_init_gps(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
