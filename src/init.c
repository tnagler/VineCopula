#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

 /* FIXME:
Check these declarations against the C/Fortran source code.
*/

 /* .C calls */
 extern void ADtest(void *, void *, void *);
 extern void archCDF(void *, void *, void *, void *, void *, void *);
 extern void ChatZj(void *, void *, void *, void *, void *, void *);
 extern void CumDist(void *, void *, void *, void *);
 extern void CvMtest(void *, void *, void *);
 extern void d2Tawn(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_mod(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_nu_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_nu_v_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_nu_v_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_par_v_mod(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_par_v_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_rho_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_rho_nu_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_rho_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_rho_v_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_v_mod(void *, void *, void *, void *, void *, void *);
 extern void diff2hfunc_v_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2lPDF_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2lPDF_rho_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2lPDF_rho_tCopula(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_mod(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_nu_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_nu_u_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_nu_u_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_nu_v_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_nu_v_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_par_u_mod(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_par_u_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_par_v_mod(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_par_v_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_rho_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_rho_nu_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_rho_tCopula(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_rho_u_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_rho_v_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_u_mod(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_u_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_v_mod(void *, void *, void *, void *, void *, void *);
 extern void diff2PDF_v_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diffhfunc_mod(void *, void *, void *, void *, void *, void *);
 extern void diffhfunc_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diffhfunc_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diffhfunc_nu_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diffhfunc_rho_tCopula(void *, void *, void *, void *, void *, void *);
 extern void diffhfunc_v_mod(void *, void *, void *, void *, void *, void *);
 extern void diffhfunc_v_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void difflPDF_mod(void *, void *, void *, void *, void *, void *);
 extern void difflPDF_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void difflPDF_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void difflPDF_nu_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void difflPDF_rho_tCopula(void *, void *, void *, void *, void *, void *);
 extern void diffPDF_mod(void *, void *, void *, void *, void *, void *);
 extern void diffPDF_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diffPDF_nu_tCopula_new(void *, void *, void *, void *, void *, void *);
 extern void diffPDF_nu_tCopula_new_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diffPDF_rho_tCopula(void *, void *, void *, void *, void *, void *);
 extern void diffPDF_u_mod(void *, void *, void *, void *, void *, void *);
 extern void diffPDF_u_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void diffPDF_v_mod(void *, void *, void *, void *, void *, void *);
 extern void diffPDF_v_mod_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void getRVM(void *, void *, void *);
 extern void gofECP(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void gofECP2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void gofECP2_pvalue(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void gofECP_pvalue(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void gofPIT_AD(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void gofPIT_AD_pvalue(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void hesse(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void Hfunc1(void *, void *, void *, void *, void *, void *, void *);
 extern void Hfunc1_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void Hfunc2(void *, void *, void *, void *, void *, void *, void *);
 extern void Hfunc2_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void Hinv1(void *, void *, void *, void *, void *, void *, void *);
 extern void Hinv1_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void Hinv2(void *, void *, void *, void *, void *, void *, void *);
 extern void Hinv2_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void KStest(void *, void *, void *);
 extern void ktau(void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void ktau_matrix(void *, void *, void *, void *);
 extern void LL_mod2(void *, void *, void *, void *, void *, void *, void *);
 extern void LL_mod_seperate(void *, void *, void *, void *, void *, void *, void *);
 extern void pcc(void *, void *, void *, void *, void *, void *, void *);
 extern void PDF_seperate(void *, void *, void *, void *, void *, void *, void *);
 extern void PDF_seperate_vec(void *, void *, void *, void *, void *, void *, void *);
 extern void RvinePIT(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void SimulateRVine(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void Tawn2(void *, void *, void *, void *, void *, void *);
 extern void TawnC(void *, void *, void *, void *, void *, void *, void *);
 extern void VineLogLikRvine(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void VineLogLikRvine2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void VineLogLikRvineGradient(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
 extern void White(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

 static const R_CMethodDef CEntries[] = {
     {"ADtest",                            (DL_FUNC) &ADtest,                             3},
     {"archCDF",                           (DL_FUNC) &archCDF,                            6},
     {"ChatZj",                            (DL_FUNC) &ChatZj,                             6},
     {"CumDist",                           (DL_FUNC) &CumDist,                            4},
     {"CvMtest",                           (DL_FUNC) &CvMtest,                            3},
     {"d2Tawn",                            (DL_FUNC) &d2Tawn,                             6},
     {"diff2hfunc_mod",                    (DL_FUNC) &diff2hfunc_mod,                     6},
     {"diff2hfunc_mod_vec",                (DL_FUNC) &diff2hfunc_mod_vec,                 7},
     {"diff2hfunc_nu_tCopula_new",         (DL_FUNC) &diff2hfunc_nu_tCopula_new,          6},
     {"diff2hfunc_nu_tCopula_new_vec",     (DL_FUNC) &diff2hfunc_nu_tCopula_new_vec,      7},
     {"diff2hfunc_nu_v_tCopula_new",       (DL_FUNC) &diff2hfunc_nu_v_tCopula_new,        6},
     {"diff2hfunc_nu_v_tCopula_new_vec",   (DL_FUNC) &diff2hfunc_nu_v_tCopula_new_vec,    7},
     {"diff2hfunc_par_v_mod",              (DL_FUNC) &diff2hfunc_par_v_mod,               6},
     {"diff2hfunc_par_v_mod_vec",          (DL_FUNC) &diff2hfunc_par_v_mod_vec,           7},
     {"diff2hfunc_rho_nu_tCopula_new",     (DL_FUNC) &diff2hfunc_rho_nu_tCopula_new,      6},
     {"diff2hfunc_rho_nu_tCopula_new_vec", (DL_FUNC) &diff2hfunc_rho_nu_tCopula_new_vec,  7},
     {"diff2hfunc_rho_tCopula_new",        (DL_FUNC) &diff2hfunc_rho_tCopula_new,         6},
     {"diff2hfunc_rho_v_tCopula_new",      (DL_FUNC) &diff2hfunc_rho_v_tCopula_new,       6},
     {"diff2hfunc_v_mod",                  (DL_FUNC) &diff2hfunc_v_mod,                   6},
     {"diff2hfunc_v_mod_vec",              (DL_FUNC) &diff2hfunc_v_mod_vec,               7},
     {"diff2lPDF_nu_tCopula_new",          (DL_FUNC) &diff2lPDF_nu_tCopula_new,           6},
     {"diff2lPDF_rho_nu_tCopula_new",      (DL_FUNC) &diff2lPDF_rho_nu_tCopula_new,       6},
     {"diff2lPDF_rho_tCopula",             (DL_FUNC) &diff2lPDF_rho_tCopula,              6},
     {"diff2PDF_mod",                      (DL_FUNC) &diff2PDF_mod,                       6},
     {"diff2PDF_mod_vec",                  (DL_FUNC) &diff2PDF_mod_vec,                   7},
     {"diff2PDF_nu_tCopula_new",           (DL_FUNC) &diff2PDF_nu_tCopula_new,            6},
     {"diff2PDF_nu_tCopula_new_vec",       (DL_FUNC) &diff2PDF_nu_tCopula_new_vec,        7},
     {"diff2PDF_nu_u_tCopula_new",         (DL_FUNC) &diff2PDF_nu_u_tCopula_new,          6},
     {"diff2PDF_nu_u_tCopula_new_vec",     (DL_FUNC) &diff2PDF_nu_u_tCopula_new_vec,      7},
     {"diff2PDF_nu_v_tCopula_new",         (DL_FUNC) &diff2PDF_nu_v_tCopula_new,          6},
     {"diff2PDF_nu_v_tCopula_new_vec",     (DL_FUNC) &diff2PDF_nu_v_tCopula_new_vec,      7},
     {"diff2PDF_par_u_mod",                (DL_FUNC) &diff2PDF_par_u_mod,                 6},
     {"diff2PDF_par_u_mod_vec",            (DL_FUNC) &diff2PDF_par_u_mod_vec,             7},
     {"diff2PDF_par_v_mod",                (DL_FUNC) &diff2PDF_par_v_mod,                 6},
     {"diff2PDF_par_v_mod_vec",            (DL_FUNC) &diff2PDF_par_v_mod_vec,             7},
     {"diff2PDF_rho_nu_tCopula_new",       (DL_FUNC) &diff2PDF_rho_nu_tCopula_new,        6},
     {"diff2PDF_rho_nu_tCopula_new_vec",   (DL_FUNC) &diff2PDF_rho_nu_tCopula_new_vec,    7},
     {"diff2PDF_rho_tCopula",              (DL_FUNC) &diff2PDF_rho_tCopula,               6},
     {"diff2PDF_rho_u_tCopula_new",        (DL_FUNC) &diff2PDF_rho_u_tCopula_new,         6},
     {"diff2PDF_rho_v_tCopula_new",        (DL_FUNC) &diff2PDF_rho_v_tCopula_new,         6},
     {"diff2PDF_u_mod",                    (DL_FUNC) &diff2PDF_u_mod,                     6},
     {"diff2PDF_u_mod_vec",                (DL_FUNC) &diff2PDF_u_mod_vec,                 7},
     {"diff2PDF_v_mod",                    (DL_FUNC) &diff2PDF_v_mod,                     6},
     {"diff2PDF_v_mod_vec",                (DL_FUNC) &diff2PDF_v_mod_vec,                 7},
     {"diffhfunc_mod",                     (DL_FUNC) &diffhfunc_mod,                      6},
     {"diffhfunc_mod_vec",                 (DL_FUNC) &diffhfunc_mod_vec,                  7},
     {"diffhfunc_nu_tCopula_new",          (DL_FUNC) &diffhfunc_nu_tCopula_new,           6},
     {"diffhfunc_nu_tCopula_new_vec",      (DL_FUNC) &diffhfunc_nu_tCopula_new_vec,       7},
     {"diffhfunc_rho_tCopula",             (DL_FUNC) &diffhfunc_rho_tCopula,              6},
     {"diffhfunc_v_mod",                   (DL_FUNC) &diffhfunc_v_mod,                    6},
     {"diffhfunc_v_mod_vec",               (DL_FUNC) &diffhfunc_v_mod_vec,                7},
     {"difflPDF_mod",                      (DL_FUNC) &difflPDF_mod,                       6},
     {"difflPDF_mod_vec",                  (DL_FUNC) &difflPDF_mod_vec,                   7},
     {"difflPDF_nu_tCopula_new",           (DL_FUNC) &difflPDF_nu_tCopula_new,            6},
     {"difflPDF_nu_tCopula_new_vec",       (DL_FUNC) &difflPDF_nu_tCopula_new_vec,        7},
     {"difflPDF_rho_tCopula",              (DL_FUNC) &difflPDF_rho_tCopula,               6},
     {"diffPDF_mod",                       (DL_FUNC) &diffPDF_mod,                        6},
     {"diffPDF_mod_vec",                   (DL_FUNC) &diffPDF_mod_vec,                    7},
     {"diffPDF_nu_tCopula_new",            (DL_FUNC) &diffPDF_nu_tCopula_new,             6},
     {"diffPDF_nu_tCopula_new_vec",        (DL_FUNC) &diffPDF_nu_tCopula_new_vec,         7},
     {"diffPDF_rho_tCopula",               (DL_FUNC) &diffPDF_rho_tCopula,                6},
     {"diffPDF_u_mod",                     (DL_FUNC) &diffPDF_u_mod,                      6},
     {"diffPDF_u_mod_vec",                 (DL_FUNC) &diffPDF_u_mod_vec,                  7},
     {"diffPDF_v_mod",                     (DL_FUNC) &diffPDF_v_mod,                      6},
     {"diffPDF_v_mod_vec",                 (DL_FUNC) &diffPDF_v_mod_vec,                  7},
     {"getRVM",                            (DL_FUNC) &getRVM,                             3},
     {"gofECP",                            (DL_FUNC) &gofECP,                            11},
     {"gofECP2",                           (DL_FUNC) &gofECP2,                           15},
     {"gofECP2_pvalue",                    (DL_FUNC) &gofECP2_pvalue,                    17},
     {"gofECP_pvalue",                     (DL_FUNC) &gofECP_pvalue,                     13},
     {"gofPIT_AD",                         (DL_FUNC) &gofPIT_AD,                         18},
     {"gofPIT_AD_pvalue",                  (DL_FUNC) &gofPIT_AD_pvalue,                  19},
     {"hesse",                             (DL_FUNC) &hesse,                             14},
     {"Hfunc1",                            (DL_FUNC) &Hfunc1,                             7},
     {"Hfunc1_vec",                        (DL_FUNC) &Hfunc1_vec,                         7},
     {"Hfunc2",                            (DL_FUNC) &Hfunc2,                             7},
     {"Hfunc2_vec",                        (DL_FUNC) &Hfunc2_vec,                         7},
     {"Hinv1",                             (DL_FUNC) &Hinv1,                              7},
     {"Hinv1_vec",                         (DL_FUNC) &Hinv1_vec,                          7},
     {"Hinv2",                             (DL_FUNC) &Hinv2,                              7},
     {"Hinv2_vec",                         (DL_FUNC) &Hinv2_vec,                          7},
     {"KStest",                            (DL_FUNC) &KStest,                             3},
     {"ktau",                              (DL_FUNC) &ktau,                               9},
     {"ktau_matrix",                       (DL_FUNC) &ktau_matrix,                        4},
     {"LL_mod2",                           (DL_FUNC) &LL_mod2,                            7},
     {"LL_mod_seperate",                   (DL_FUNC) &LL_mod_seperate,                    7},
     {"pcc",                               (DL_FUNC) &pcc,                                7},
     {"PDF_seperate",                      (DL_FUNC) &PDF_seperate,                       7},
     {"PDF_seperate_vec",                  (DL_FUNC) &PDF_seperate_vec,                   7},
     {"RvinePIT",                          (DL_FUNC) &RvinePIT,                          14},
     {"SimulateRVine",                     (DL_FUNC) &SimulateRVine,                     11},
     {"Tawn2",                             (DL_FUNC) &Tawn2,                              6},
     {"TawnC",                             (DL_FUNC) &TawnC,                              7},
     {"VineLogLikRvine",                   (DL_FUNC) &VineLogLikRvine,                   16},
     {"VineLogLikRvine2",                  (DL_FUNC) &VineLogLikRvine2,                  11},
     {"VineLogLikRvineGradient",           (DL_FUNC) &VineLogLikRvineGradient,           15},
     {"White",                             (DL_FUNC) &White,                             12},
     {NULL, NULL, 0}
 };

 void R_init_VineCopula(DllInfo *dll)
 {
     R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
     R_useDynamicSymbols(dll, FALSE);
     R_RegisterCCallable("VineCopula", "archCDF",                           (DL_FUNC) &archCDF);
     R_RegisterCCallable("VineCopula", "diff2hfunc_mod",                    (DL_FUNC) &diff2hfunc_mod);
     R_RegisterCCallable("VineCopula", "diff2hfunc_mod_vec",                (DL_FUNC) &diff2hfunc_mod_vec);
     R_RegisterCCallable("VineCopula", "diff2hfunc_nu_tCopula_new",         (DL_FUNC) &diff2hfunc_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2hfunc_nu_tCopula_new_vec",     (DL_FUNC) &diff2hfunc_nu_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diff2hfunc_nu_v_tCopula_new",       (DL_FUNC) &diff2hfunc_nu_v_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2hfunc_nu_v_tCopula_new_vec",   (DL_FUNC) &diff2hfunc_nu_v_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diff2hfunc_par_v_mod",              (DL_FUNC) &diff2hfunc_par_v_mod);
     R_RegisterCCallable("VineCopula", "diff2hfunc_par_v_mod_vec",          (DL_FUNC) &diff2hfunc_par_v_mod_vec);
     R_RegisterCCallable("VineCopula", "diff2hfunc_rho_nu_tCopula_new",     (DL_FUNC) &diff2hfunc_rho_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2hfunc_rho_nu_tCopula_new_vec", (DL_FUNC) &diff2hfunc_rho_nu_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diff2hfunc_rho_tCopula_new",        (DL_FUNC) &diff2hfunc_rho_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2hfunc_rho_v_tCopula_new",      (DL_FUNC) &diff2hfunc_rho_v_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2hfunc_v_mod",                  (DL_FUNC) &diff2hfunc_v_mod);
     R_RegisterCCallable("VineCopula", "diff2hfunc_v_mod_vec",              (DL_FUNC) &diff2hfunc_v_mod_vec);
     R_RegisterCCallable("VineCopula", "diff2lPDF_nu_tCopula_new",          (DL_FUNC) &diff2lPDF_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2lPDF_rho_nu_tCopula_new",      (DL_FUNC) &diff2lPDF_rho_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2lPDF_rho_tCopula",             (DL_FUNC) &diff2lPDF_rho_tCopula);
     R_RegisterCCallable("VineCopula", "diff2PDF_mod",                      (DL_FUNC) &diff2PDF_mod);
     R_RegisterCCallable("VineCopula", "diff2PDF_mod_vec",                  (DL_FUNC) &diff2PDF_mod_vec);
     R_RegisterCCallable("VineCopula", "diff2PDF_nu_tCopula_new",           (DL_FUNC) &diff2PDF_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2PDF_nu_tCopula_new_vec",       (DL_FUNC) &diff2PDF_nu_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diff2PDF_nu_u_tCopula_new",         (DL_FUNC) &diff2PDF_nu_u_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2PDF_nu_u_tCopula_new_vec",     (DL_FUNC) &diff2PDF_nu_u_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diff2PDF_nu_v_tCopula_new",         (DL_FUNC) &diff2PDF_nu_v_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2PDF_nu_v_tCopula_new_vec",     (DL_FUNC) &diff2PDF_nu_v_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diff2PDF_par_u_mod",                (DL_FUNC) &diff2PDF_par_u_mod);
     R_RegisterCCallable("VineCopula", "diff2PDF_par_u_mod_vec",            (DL_FUNC) &diff2PDF_par_u_mod_vec);
     R_RegisterCCallable("VineCopula", "diff2PDF_par_v_mod",                (DL_FUNC) &diff2PDF_par_v_mod);
     R_RegisterCCallable("VineCopula", "diff2PDF_par_v_mod_vec",            (DL_FUNC) &diff2PDF_par_v_mod_vec);
     R_RegisterCCallable("VineCopula", "diff2PDF_rho_nu_tCopula_new",       (DL_FUNC) &diff2PDF_rho_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2PDF_rho_nu_tCopula_new_vec",   (DL_FUNC) &diff2PDF_rho_nu_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diff2PDF_rho_tCopula",              (DL_FUNC) &diff2PDF_rho_tCopula);
     R_RegisterCCallable("VineCopula", "diff2PDF_rho_u_tCopula_new",        (DL_FUNC) &diff2PDF_rho_u_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2PDF_rho_v_tCopula_new",        (DL_FUNC) &diff2PDF_rho_v_tCopula_new);
     R_RegisterCCallable("VineCopula", "diff2PDF_u_mod",                    (DL_FUNC) &diff2PDF_u_mod);
     R_RegisterCCallable("VineCopula", "diff2PDF_u_mod_vec",                (DL_FUNC) &diff2PDF_u_mod_vec);
     R_RegisterCCallable("VineCopula", "diff2PDF_v_mod",                    (DL_FUNC) &diff2PDF_v_mod);
     R_RegisterCCallable("VineCopula", "diff2PDF_v_mod_vec",                (DL_FUNC) &diff2PDF_v_mod_vec);
     R_RegisterCCallable("VineCopula", "diffhfunc_mod",                     (DL_FUNC) &diffhfunc_mod);
     R_RegisterCCallable("VineCopula", "diffhfunc_mod_vec",                 (DL_FUNC) &diffhfunc_mod_vec);
     R_RegisterCCallable("VineCopula", "diffhfunc_nu_tCopula_new",          (DL_FUNC) &diffhfunc_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "diffhfunc_nu_tCopula_new_vec",      (DL_FUNC) &diffhfunc_nu_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diffhfunc_rho_tCopula",             (DL_FUNC) &diffhfunc_rho_tCopula);
     R_RegisterCCallable("VineCopula", "diffhfunc_v_mod",                   (DL_FUNC) &diffhfunc_v_mod);
     R_RegisterCCallable("VineCopula", "diffhfunc_v_mod_vec",               (DL_FUNC) &diffhfunc_v_mod_vec);
     R_RegisterCCallable("VineCopula", "difflPDF_mod",                      (DL_FUNC) &difflPDF_mod);
     R_RegisterCCallable("VineCopula", "difflPDF_mod_vec",                  (DL_FUNC) &difflPDF_mod_vec);
     R_RegisterCCallable("VineCopula", "difflPDF_nu_tCopula_new",           (DL_FUNC) &difflPDF_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "difflPDF_nu_tCopula_new_vec",       (DL_FUNC) &difflPDF_nu_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "difflPDF_rho_tCopula",              (DL_FUNC) &difflPDF_rho_tCopula);
     R_RegisterCCallable("VineCopula", "diffPDF_mod",                       (DL_FUNC) &diffPDF_mod);
     R_RegisterCCallable("VineCopula", "diffPDF_mod_vec",                   (DL_FUNC) &diffPDF_mod_vec);
     R_RegisterCCallable("VineCopula", "diffPDF_nu_tCopula_new",            (DL_FUNC) &diffPDF_nu_tCopula_new);
     R_RegisterCCallable("VineCopula", "diffPDF_nu_tCopula_new_vec",        (DL_FUNC) &diffPDF_nu_tCopula_new_vec);
     R_RegisterCCallable("VineCopula", "diffPDF_rho_tCopula",               (DL_FUNC) &diffPDF_rho_tCopula);
     R_RegisterCCallable("VineCopula", "diffPDF_u_mod",                     (DL_FUNC) &diffPDF_u_mod);
     R_RegisterCCallable("VineCopula", "diffPDF_u_mod_vec",                 (DL_FUNC) &diffPDF_u_mod_vec);
     R_RegisterCCallable("VineCopula", "diffPDF_v_mod",                     (DL_FUNC) &diffPDF_v_mod);
     R_RegisterCCallable("VineCopula", "diffPDF_v_mod_vec",                 (DL_FUNC) &diffPDF_v_mod_vec);
     R_RegisterCCallable("VineCopula", "Hfunc1",                            (DL_FUNC) &Hfunc1);
     R_RegisterCCallable("VineCopula", "Hfunc1_vec",                        (DL_FUNC) &Hfunc1_vec);
     R_RegisterCCallable("VineCopula", "Hfunc2",                            (DL_FUNC) &Hfunc2);
     R_RegisterCCallable("VineCopula", "Hfunc2_vec",                        (DL_FUNC) &Hfunc2_vec);
     R_RegisterCCallable("VineCopula", "Hinv1",                             (DL_FUNC) &Hinv1);
     R_RegisterCCallable("VineCopula", "Hinv1_vec",                         (DL_FUNC) &Hinv1_vec);
     R_RegisterCCallable("VineCopula", "Hinv2",                             (DL_FUNC) &Hinv2);
     R_RegisterCCallable("VineCopula", "Hinv2_vec",                         (DL_FUNC) &Hinv2_vec);
     R_RegisterCCallable("VineCopula", "ktau",                              (DL_FUNC) &ktau);
     R_RegisterCCallable("VineCopula", "ktau_matrix",                       (DL_FUNC) &ktau_matrix);
     R_RegisterCCallable("VineCopula", "LL_mod2",                           (DL_FUNC) &LL_mod2);
     R_RegisterCCallable("VineCopula", "LL_mod_seperate",                   (DL_FUNC) &LL_mod_seperate);
     R_RegisterCCallable("VineCopula", "pcc",                               (DL_FUNC) &pcc);
     R_RegisterCCallable("VineCopula", "PDF_seperate",                      (DL_FUNC) &PDF_seperate);
     R_RegisterCCallable("VineCopula", "PDF_seperate_vec",                  (DL_FUNC) &PDF_seperate_vec);
     R_RegisterCCallable("VineCopula", "RvinePIT",                          (DL_FUNC) &RvinePIT);
     R_RegisterCCallable("VineCopula", "SimulateRVine",                     (DL_FUNC) &SimulateRVine);
     R_RegisterCCallable("VineCopula", "Tawn2",                             (DL_FUNC) &Tawn2);
     R_RegisterCCallable("VineCopula", "TawnC",                             (DL_FUNC) &TawnC);
     R_RegisterCCallable("VineCopula", "VineLogLikRvine",                   (DL_FUNC) &VineLogLikRvine);
     R_RegisterCCallable("VineCopula", "VineLogLikRvine2",                  (DL_FUNC) &VineLogLikRvine2);
     R_RegisterCCallable("VineCopula", "VineLogLikRvineGradient",           (DL_FUNC) &VineLogLikRvineGradient);
}
