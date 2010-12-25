#ifndef AOS_LIB_STFUN
#define AOS_LIB_STFUN
typedef struct stfun_t stfun_t;
void stfun_init(stfun_t *A, long nx, long ny, double *amp);
void stfun_push(stfun_t *A, dmat *opd);
dmat *stfun_finalize(stfun_t *A);
#endif