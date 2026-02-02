/**
 * @brief Number of free elements in Q_idm×Q_idm Cholesky factor L_i.
 * Uses Stan integer division operator %/% (parenthesize to avoid precedence issues).
 */
int M_cov = (indep_idmarker_cov == 1) ? Q_idm
            : ((Q_idm * (Q_idm + 1)) %/% 2);

// mapping m -> (r,c)
array[M_cov] int r_idx;
array[M_cov] int c_idx;

if (indep_idmarker_cov == 1) {
  for (m in 1 : M_cov) {
    r_idx[m] = m;
    c_idx[m] = m;
  }
} else {
  int idx = 1;
  for (r in 1 : Q_idm)
    for (c in 1 : r) {
      r_idx[idx] = r;
      c_idx[idx] = c;
      idx += 1;
    }
}
