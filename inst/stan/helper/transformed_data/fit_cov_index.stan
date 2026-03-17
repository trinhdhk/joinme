/**
 * @brief Number of free elements in Q_idm x Q_idm Cholesky factor L_i.
 * Uses Stan integer division operator %/% (parenthesize to avoid precedence issues).
 *
 * When indep_idmarker_cov == 1, only the diagonal is modeled (M_cov = Q_idm).
 * Otherwise, include the full lower triangle (M_cov = Q_idm*(Q_idm+1)/2).
 */
int M_cov = (indep_idmarker_cov == 1) ? Q_idm
            : ((Q_idm * (Q_idm + 1)) %/% 2);

/* number of off-diagonal correlation association features */
int M_corr = (Q_idm >= 2) ? ((Q_idm * (Q_idm - 1)) %/% 2) : 0;

/* number of vcov association features (lower-triangular entries of L) */
int M_vcov = M_cov;

/* mapping m -> (r,c) in the lower triangle of L_i */
array[M_cov] int r_idx; // row index for element m
array[M_cov] int c_idx; // column index for element m

/* Diagonal-only path: m indexes (r=r, c=r) */
if (indep_idmarker_cov == 1) {
  for (m in 1 : M_cov) {
    r_idx[m] = m;
    c_idx[m] = m;
  }
} else {
  /* Full lower triangle path: row-major within each r */
  int idx = 1;
  for (r in 1 : Q_idm)
    for (c in 1 : r) {
      r_idx[idx] = r;
      c_idx[idx] = c;
      idx += 1;
    }
}
