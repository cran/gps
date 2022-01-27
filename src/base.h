static inline void ZeroVec (int n, double *x) {
  double *xi = x, *xn = x + n;
  if (n > 0) xi[0] = 0;
  for (xi += n & 1; xi < xn; xi += 2) {
    xi[0] = 0; xi[1] = 0;
  }
}
