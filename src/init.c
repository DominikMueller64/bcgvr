#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP embvr_dhbv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP embvr_embv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"embvr_dhbv", (DL_FUNC) &embvr_dhbv,  7},
  {"embvr_embv", (DL_FUNC) &embvr_embv, 10},
  {NULL, NULL, 0}
};

void R_init_embvr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
