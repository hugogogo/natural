#ifndef _COORD_H
#define _COORD_H

// matrix functions
void ComputeColumnMeans(double *, int, int, double *);

// Lasso related functions
static inline double SoftThreshold(const double c, const double t) {
  /* Soft thresholds c by t. */
  return (c < -t) ? (c + t) : ((c > t) ?  (c - t) : 0.0);
}

double ComputeObjective(const double *, const int, const int, const cs *, const double);
void cd_update(const double *, const double *, const int, const int, const int, const double, const int,
    double, double *, cs *, int *, int *, int *,
    int *, int *, const double);

// active set functions
int IsInArray(const int *, const int, const int);
int AddToArray(int *, int*, const int, const int);

// Sparse matrix related functions
int cs_GetElement(cs *, int, int, double *);
int cs_Predict(double *, int, int, cs *, double *);

#endif
