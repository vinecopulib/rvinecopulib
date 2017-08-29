#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _rvinecopulib_bicop_cdf_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_bicop_check_cpp(SEXP);
extern SEXP _rvinecopulib_bicop_hfunc1_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_bicop_hfunc2_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_bicop_hinv1_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_bicop_hinv2_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_bicop_loglik_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_bicop_par_to_tau_cpp(SEXP);
extern SEXP _rvinecopulib_bicop_pdf_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_bicop_select_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rvinecopulib_bicop_tau_to_par_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_rvine_matrix_check_cpp(SEXP);
extern SEXP _rvinecopulib_vinecop_cdf_cpp(SEXP, SEXP, SEXP);
extern SEXP _rvinecopulib_vinecop_check_cpp(SEXP);
extern SEXP _rvinecopulib_vinecop_inverse_rosenblatt_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_vinecop_loglik_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_vinecop_pdf_cpp(SEXP, SEXP);
extern SEXP _rvinecopulib_vinecop_select_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_rvinecopulib_bicop_cdf_cpp",                  (DL_FUNC) &_rvinecopulib_bicop_cdf_cpp,                   2},
    {"_rvinecopulib_bicop_check_cpp",                (DL_FUNC) &_rvinecopulib_bicop_check_cpp,                 1},
    {"_rvinecopulib_bicop_hfunc1_cpp",               (DL_FUNC) &_rvinecopulib_bicop_hfunc1_cpp,                2},
    {"_rvinecopulib_bicop_hfunc2_cpp",               (DL_FUNC) &_rvinecopulib_bicop_hfunc2_cpp,                2},
    {"_rvinecopulib_bicop_hinv1_cpp",                (DL_FUNC) &_rvinecopulib_bicop_hinv1_cpp,                 2},
    {"_rvinecopulib_bicop_hinv2_cpp",                (DL_FUNC) &_rvinecopulib_bicop_hinv2_cpp,                 2},
    {"_rvinecopulib_bicop_loglik_cpp",               (DL_FUNC) &_rvinecopulib_bicop_loglik_cpp,                2},
    {"_rvinecopulib_bicop_par_to_tau_cpp",           (DL_FUNC) &_rvinecopulib_bicop_par_to_tau_cpp,            1},
    {"_rvinecopulib_bicop_pdf_cpp",                  (DL_FUNC) &_rvinecopulib_bicop_pdf_cpp,                   2},
    {"_rvinecopulib_bicop_select_cpp",               (DL_FUNC) &_rvinecopulib_bicop_select_cpp,                7},
    {"_rvinecopulib_bicop_tau_to_par_cpp",           (DL_FUNC) &_rvinecopulib_bicop_tau_to_par_cpp,            2},
    {"_rvinecopulib_rvine_matrix_check_cpp",         (DL_FUNC) &_rvinecopulib_rvine_matrix_check_cpp,          1},
    {"_rvinecopulib_vinecop_cdf_cpp",                (DL_FUNC) &_rvinecopulib_vinecop_cdf_cpp,                 3},
    {"_rvinecopulib_vinecop_check_cpp",              (DL_FUNC) &_rvinecopulib_vinecop_check_cpp,               1},
    {"_rvinecopulib_vinecop_inverse_rosenblatt_cpp", (DL_FUNC) &_rvinecopulib_vinecop_inverse_rosenblatt_cpp,  2},
    {"_rvinecopulib_vinecop_loglik_cpp",             (DL_FUNC) &_rvinecopulib_vinecop_loglik_cpp,              2},
    {"_rvinecopulib_vinecop_pdf_cpp",                (DL_FUNC) &_rvinecopulib_vinecop_pdf_cpp,                 2},
    {"_rvinecopulib_vinecop_select_cpp",             (DL_FUNC) &_rvinecopulib_vinecop_select_cpp,             14},
    {NULL, NULL, 0}
};

void R_init_rvinecopulib(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
