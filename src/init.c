#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
extern SEXP C_ApproxDR(void *, void *, void *, void *, void *);
extern SEXP C_btSb(void *, void *);
extern SEXP C_CheckSupp(void *, void *);
extern SEXP C_ComputeLD(void *, void *);
extern SEXP C_Csparse2LTB(void *, void *, void *, void *);
extern SEXP C_DDt(void *);
extern SEXP C_Diff(void *, void *, void *, void *);
extern SEXP C_DtD(void *);
extern SEXP C_FormE(void *, void *);
extern SEXP C_GridGCV(void *, void *, void *, void *, void *, void *);
extern SEXP C_GridPWLS(void *, void *, void *, void *, void *, void *);
extern SEXP C_IsMonoInc(void *, void *, void *);
extern SEXP C_LAUUM(void *);
extern SEXP C_LPBTRF(void *, void *);
extern SEXP C_LTB2Dense(void *, void *);
extern SEXP C_MakeGrid(void *, void *, void *);
extern SEXP C_MaxDR(void *, void *, void *);
extern SEXP C_MeanDR(void *);
extern SEXP C_MinDR(void *, void *);
extern SEXP C_NullGD(void *, void *);
extern SEXP C_REDF2Rho(void *, void *, void *, void *, void *);
extern SEXP C_Rho2REDF(void *, void *);
extern SEXP C_RootApproxDR(void *, void *, void *, void *, void *);
extern SEXP C_SbarBlocks(void *, void *, void *);
extern SEXP C_SbarLTB(void *);
extern SEXP C_SetDim(void *, void *);
extern SEXP C_SolveLTB(void *, void *, void *, void *);
extern SEXP C_VecDot(void *, void *);
extern SEXP C_VecScal(void *, void *, void *);
static const R_CallMethodDef CallEntries[] = {
    {"C_ApproxDR", (DL_FUNC) &C_ApproxDR, 5},
    {"C_btSb", (DL_FUNC) &C_btSb, 2},
    {"C_CheckSupp", (DL_FUNC) &C_CheckSupp, 2},
    {"C_ComputeLD", (DL_FUNC) &C_ComputeLD, 2},
    {"C_Csparse2LTB", (DL_FUNC) &C_Csparse2LTB, 4},
    {"C_DDt", (DL_FUNC) &C_DDt, 1},
    {"C_Diff", (DL_FUNC) &C_Diff, 4},
    {"C_DtD", (DL_FUNC) &C_DtD, 1},
    {"C_FormE", (DL_FUNC) &C_FormE, 2},
    {"C_GridGCV", (DL_FUNC) &C_GridGCV, 6},
    {"C_GridPWLS", (DL_FUNC) &C_GridPWLS, 6},
    {"C_IsMonoInc", (DL_FUNC) &C_IsMonoInc, 3},
    {"C_LAUUM", (DL_FUNC) &C_LAUUM, 1},
    {"C_LPBTRF", (DL_FUNC) &C_LPBTRF, 2},
    {"C_LTB2Dense", (DL_FUNC) &C_LTB2Dense, 2},
    {"C_MakeGrid", (DL_FUNC) &C_MakeGrid, 3},
    {"C_MaxDR", (DL_FUNC) &C_MaxDR, 3},
    {"C_MeanDR", (DL_FUNC) &C_MeanDR, 1},
    {"C_MinDR", (DL_FUNC) &C_MinDR, 2},
    {"C_NullGD", (DL_FUNC) &C_NullGD, 2},
    {"C_REDF2Rho", (DL_FUNC) &C_REDF2Rho, 5},
    {"C_Rho2REDF", (DL_FUNC) &C_Rho2REDF, 2},
    {"C_RootApproxDR", (DL_FUNC) &C_RootApproxDR, 5},
    {"C_SbarBlocks", (DL_FUNC) &C_SbarBlocks, 3},
    {"C_SbarLTB", (DL_FUNC) &C_SbarLTB, 1},
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
