#ifndef __C_INTERFACE_H__
#define __C_INTERFACE_H__

#ifdef __cplusplus
extern "C"
{
#endif

int estimate_circle_taubin(double x[], double y[], int length, double *out_x, double *out_y, double *out_r);

int estimate_circle_lm(double x[], double y[], int length, double rxy_init[3], double rxy_out[3]);

int estimate_circle(double x[], double y[], int length, double *out_x, double *out_y, double *out_r);

int estimate_circle_lm_trace(double x[], double y[], int length, double rxy_init[3], int *out_numiters, double **out_x, double **out_y, double **out_r, double **out_grad_norm, double **out_mse);

#ifdef __cplusplus
}
#endif

#endif // __C_INTERFACE_H__

