/*
Based on coord7 from Jacob's code.
Organic lasso
does not standardize matrix or response, i.e., no estimate of intercept
*/

#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "csparse.h"
#include "olasso.h"

// This function
// arguments:
//    x: nobs -by- ndim design matrix, not standardized
//    y:    an nobs-vector of response
//    lambda: an nlambda-dim vector storing the list of lambdas used
//    gamma_i, gamma_p, gamma_x: tuples of sparse matrix representation of the resulting estimates of gamma
//                               which is a ndim -by- nlambda matrix, with each column representing an estimate of gamma for one lambda
//
//    It solves:
//    min_\beta { 0.5/n ||y - X\beta||^2 + \lambda ||\beta||_1^2}
//
void olasso_c(double *x, double *y, int *nobs, int *ndim,
    int *nlambda, double *lambda,
    int *beta_i, int *beta_p, double *beta_x, int *maxnz, double *thresh){

  // sample size n
  const int n = *nobs;
  // dimension p
  const int p = *ndim;
  // number of tuning parameters
  const int nlam = *nlambda;
  // maximum number of nonzero elements
  const int nz = *maxnz;
  // iterators, more specifically,
  // i used for iterating row, j used for iterating column, l used for iterating lambda
  // ll is used for iterating nonzero elements in the sparse matrix representation
  int i, j, l, ll;
  // constants used in active sets strategy
  // cycle_count is the number of cycles of coordinates updates for each lambda
  // nelt is constant used in warm-starting strategy
  // number of nonzero elements in previous column (previous lambda)
  int cycle_count, last_changed_cycle, all_var_update_mode, nelt;
  double tol = *thresh;
  double temp, nulldev;

  // residuals
  double *r = malloc(n * sizeof(double));
  // l2norm of each columns
  double *sd = malloc(p * sizeof(double));

  // temporarily store the value of the l1-norm of beta
  double l1norm;

  // active sets
  int *act = malloc(p * sizeof(int));
  int act_length = 0;

  // let sp_beta's i,p,x point memory that's been allocated by R.
  // sp_beta is used to store the result of the third step
  cs *sp_beta = cs_calloc(1, sizeof(cs)); // let sp_beta's i,p,x point memory that's been allocated by R.
  sp_beta->m = p;
  sp_beta->n = nlam;
  sp_beta->nzmax = nz;
  sp_beta->nz = -1; // indicates it's in compressed column form
  sp_beta->i = beta_i;
  sp_beta->p = beta_p;
  sp_beta->x = beta_x;

  // initialize residual
  for (i = 0; i < n; i++){
    r[i] = y[i];
  }

  // Does not standardize design matrix
  // just store ||x_j||^2 / n into sd
  for (j = 0; j < p; j++){
    temp = 0.0;
    for (i = 0; i < n; i++)
      temp += pow(x[i + n * j], 2);
    sd[j] = temp / n;
  }

  // Compute null deviance:
  nulldev = 0.0;
  // Modification: previously it was y[i], but really should be centered y if we model intercept
  for (i = 0; i < n; i++)
    nulldev += pow(r[i], 2);

  // update tolerance after data standardization
  tol = tol / (sqrt(n / nulldev));
  // Rprintf("Tolerance: %e", tol);

  act_length = 0;
  last_changed_cycle = 0;
  all_var_update_mode = 1;
  cycle_count = 0;

  for (l = 1; l < nlam; l++) {
    // begin with warm start from previous lambda -- copy column l - 1 to column l:
    nelt = sp_beta->p[l] - sp_beta->p[l - 1]; // number of non-zero elements in col l - 1
    l1norm = 0.0;
    for (ll = sp_beta->p[l - 1]; ll < sp_beta->p[l]; ll++) {
      // (loop over the elements in column l - 1)
      sp_beta->x[nelt + ll] = sp_beta->x[ll];
      sp_beta->i[nelt + ll] = sp_beta->i[ll];
      l1norm += fabs(sp_beta->x[ll]);
    }
    for (ll = l + 1; ll < nlam + 1; ll++)
      sp_beta->p[ll] += nelt;

    // Main coordinate descent loop:
    cd_update(x, sd, n, p, l, lambda[l], nlam,
        l1norm, r, sp_beta, act, &act_length,
        &all_var_update_mode, &last_changed_cycle, &cycle_count, tol);
  }

  // de-allocate memory:
  free(r);
  free(sd);
  free(act);
  // note: didn't dynamically allocate i, p, or x... so don't call cs_spfree.
  cs_free(sp_beta);
}

void cd_update(const double *x, const double *sd, const int nobs, const int ndim, const int l, const double lam, const int ncol,
    double l1norm, double *r, cs *sp_beta, int *act, int *act_length,
    int *all_var_update_mode, int *last_changed_cycle, int *cycle_count, const double tol) {
//  This function takes a nobs-by-ndim matrix x, a nobs-vector r, and
//  computes the lasso with tuning parameter lambda from a warm start sp_beta
//  and with active set initialized to act.
//
//  The function updates sp_beta (passing the result by reference).
//
//  Args:
//    x: nobs-by-ndim standardized design matrix with column mean zero and column standard-deviation 1.
//    l: the column number in the resulting solution matrix
//    lam: tuning parameter
//    ncol: total number of columns in the solution matrix
//    r: residuals: initialized by residuals from previous tuning parameter

  int i, j, k, ll;
  double cxr, oldth, newth, delta;
  while (1) {
    // if we need to soft-threshold all variables
    if (*all_var_update_mode) {
      for (j = 0; j < ndim; j++) {
        // NAIVE UPDATE: compute x_j^Tr
        cxr = 0;
        for (i = 0; i < nobs; i++)
          cxr += x[i + nobs * j] * r[i];

        ll = cs_GetElement(sp_beta, j, l, &oldth);
        // l1norm = ||beta_{-j}||_1
        l1norm -= fabs(oldth);
        // Soft-Threshold
        newth = SoftThreshold(sd[j] * oldth + cxr / nobs, 2 * lam * l1norm) / (2 * lam + sd[j]);
        // update l1norm
        l1norm += fabs(newth);

        if ((ll == -1) && (newth != 0)) {
          // note: this works because l is the last column with nonzero entries
          sp_beta->i[sp_beta->p[l + 1]] = j;
          sp_beta->x[sp_beta->p[l + 1]] = newth;
          for (ll = l + 1; ll < ncol + 1; ll++)
            sp_beta->p[ll]++;
        }
        else if (ll >= 0){
          sp_beta->x[ll] = newth;
        }

        delta = newth - oldth;

        if (delta != 0) {
          // parameter just changed
          if ((oldth == 0) & !IsInArray(act, *act_length, j)) {
            // add main effect to ever-active set
            AddToArray(act, act_length, ndim, j);
          }
          // update residual:
          for (i = 0; i < nobs; i++)
            r[i] -= delta * x[i + nobs * j];
          if (fabs(delta) > tol)
            *last_changed_cycle = *cycle_count;
        }
      }
    }
    // otherwise, we only need to update those in active set
    else {
      for (k = 0; k < *act_length; k++) {
        j = act[k];
        // NAIVE UPDATE: compute x_j^Tr
        cxr = 0;
        for (i = 0; i < nobs; i++)
          cxr += x[i + nobs * j] * r[i];

        ll = cs_GetElement(sp_beta, j, l, &oldth);
        // l1norm = ||beta_{-j}||_1
        l1norm -= fabs(oldth);
        // Soft-Threshold
        newth = SoftThreshold(oldth * sd[j] + cxr / nobs, 2 * lam * l1norm) / (2 * lam + sd[j]);
        // update l1norm
        l1norm += fabs(newth);

        if ((ll == -1) && (newth != 0)) {
          // note: this works because l is the last column with nonzero entries
          sp_beta->i[sp_beta->p[l + 1]] = j;
          sp_beta->x[sp_beta->p[l + 1]] = newth;
          for (ll = l + 1; ll < ncol + 1; ll++)
            sp_beta->p[ll]++;
        }
        else if (ll >= 0){
          sp_beta->x[ll] = newth;
        }

        delta = newth - oldth;

        if (delta != 0) {
          // parameter just changed
          if ((oldth == 0) & !IsInArray(act, *act_length, j)) {
            // add main effect to ever-active set
            AddToArray(act, act_length, ndim, j);
          }
          // update residual:
          for (i = 0; i < nobs; i++)
            r[i] -= delta * x[i + nobs * j];
          if (fabs(delta) > tol)
            *last_changed_cycle = *cycle_count;
        }
      }
    }

    // one iteration of coordinates updates finished
    if ((*cycle_count > 0) && (*last_changed_cycle < *cycle_count)) {
      if (*all_var_update_mode) {
        *all_var_update_mode = 0;
        // Rprintf("Total number of cycles: %d \n", *cycle_count);
        break;
      }
      else {
        // converged on active set... now do pass over all variables.
        *all_var_update_mode = 1;
      }
    } else {
      // we haven't yet converged
      if (*all_var_update_mode) {
        *all_var_update_mode = 0; // just did a pass over all vars - switch back to everactive
      }
    }
    (*cycle_count)++;
  }
}

void ComputeColumnMeans(double *m, int nrow, int ncol, double *cmeans) {
  int i, j, jj;
  for (j = 0; j < ncol; j++) {
    cmeans[j] = 0;
    jj = nrow * j;
    for (i = 0; i < nrow; i++) {
      cmeans[j] += m[i + jj] / nrow;
    }
  }
}

// Active set functions
int IsInArray(const int *arr, const int len, const int j) {
  int i;
  for (i = 0; i < len; i++) {
    if (arr[i] == j)
      return 1;
  }
  return 0;
}

int AddToArray(int *arr, int *len, const int maxlen, const int j) {
  if (*len >= maxlen) {
    Rprintf("could not add element.  Array full.\n");
    return 0;
  }
  else {
    arr[*len] = j;
    (*len)++;
    return 1;
  }
}

int cs_GetElement(cs *A, int row, int col, double *element) {
  /* Returns by reference the (first occurence of) the element A[row, col].
     If A is in triplet form, returns the
     index ll such that A->i[ll] == row and A->p[ll] == col.  If A is in compressed column form
     returns ll such that A->i[ll] == row and A->p[col] <= ll < A->p[col + 1].

     Returns -1 if element is not in A (in which case *element = 0). */
  int ll;
  if (row >= A->m || col >= A->n) {
    Rprintf("Error: row or col exceeds A->m and A->n!\n");
    *element = 0;
    return -1;
  }
  if (A->nz == -1) {
    // compressed column form
    for (ll = A->p[col]; ll < A->p[col + 1]; ll++)
      if (A->i[ll] == row) {
        *element = A->x[ll];
        return ll;
      }
    *element = 0;
    return -1;
  }
  // triplet form
  for (ll = 0; ll < A->nz; ll++)
    if ((A->i[ll] == row) && (A->p[ll] == col)) {
      *element = A->x[ll];
      return ll;
    }
  *element = 0;
  return -1;
}

int cs_Predict(double *B, int nn, int np, cs *A, double *C) {
  // Computes C = C + B_aug * A^T, where A is a sparse matrix and B is a dense design matrix
  // B_aug is the augmented design matrix
  // Args:
  // A: sparse matrix in compressed column format with nq columns.
  // B: a dense nn-by-np matrix of standardized design matrix.

  if (A->n != np || A->nz != -1)
    return 0;
  int row, col, l;
  for (col = 0; col < np; col++){
    // main effects
    for (l = A->p[col]; l < A->p[col + 1]; l++)
      for (row = 0; row < nn; row++)
        C[row + nn * A->i[l]] += A->x[l] * B[row + nn * col];
  }
  return 1;
}

double ComputeObjective(const double *r, const int n, const int l, const cs *vec, const double lambda) {
  // Computes ||r||^2 /(n) + 2 * lambda * ||th||_1^2
  int i, ll;
  double obj = 0;
  double norm = 0;
  for (i = 0; i < n; i++)
    obj += pow(r[i], 2);
  obj /= n;
  // (loop over the elements in column l)
  for (ll = vec->p[l]; ll < vec->p[l + 1]; ll++)
    norm += fabs(vec->x[ll]);
  obj += 2 * lambda * norm * norm;
  return obj;
}

int main(int argc, char** argv)
{
  return 1;
}
