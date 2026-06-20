/*
 * scatter2grid.i — IDW interpolation from scattered/hex 2D data to Cartesian grid
 *
 * func scatter2grid(xi, yi, zi, nx, ny, power=, nneighbors=, xout=, yout=)
 *
 * Inputs:
 *   xi, yi      - 1D arrays of input point coordinates
 *   zi          - 1D array of values at those points
 *   nx, ny      - output grid dimensions
 *
 * Keywords:
 *   power=      - IDW distance weighting exponent (default: 2.0)
 *   nneighbors= - number of nearest neighbours to use (default: 7)
 *   xout=       - if provided, named variable receives the output x axis (1D, length nx)
 *   yout=       - if provided, named variable receives the output y axis (1D, length ny)
 *
 * Returns:
 *   2D array [nx, ny] of interpolated values on a regular Cartesian grid
 *   Grid extent is the bounding box of (xi, yi) with a small margin.
 *
 * Notes:
 *   - For a hexagonal grid, nneighbors=7 (centre + 6 hex neighbours) is natural.
 *   - power=2 is standard IDW; lower values (1) give smoother results,
 *     higher values (3-4) make the interpolant more local/faithful to samples.
 *   - Exact hits (distance < eps) return the sample value directly.
 * 
 * Note: Coded by Claude.
 */
func scatter2grid(xi, yi, zi, nx, ny, &xout, &yout, power=, nneighbors=)
{
  if (is_void(power))      power      = 2.0;
  if (is_void(nneighbors)) nneighbors = 7;

  npts = numberof(xi);
  if (numberof(yi) != npts || numberof(zi) != npts)
    error, "xi, yi, zi must have the same number of elements";
  if (nneighbors > npts) nneighbors = npts;

  // Build output grid with a small margin (5% of range) so edge points aren't clipped
  xmin = min(xi); xmax = max(xi);
  ymin = min(yi); ymax = max(yi);
  xmargin = (xmax - xmin) * 0.05;
  ymargin = (ymax - ymin) * 0.05;
  xout = span(xmin - xmargin, xmax + xmargin, nx);
  yout = span(ymin - ymargin, ymax + ymargin, ny);

  // Output grid as 2D coordinate arrays [nx, ny]
  // xg varies along first dim, yg along second
  xg = xout(, -:1:ny);   // [nx, ny]
  yg = yout(-:1:nx, );   // [nx, ny]

  result = array(0.0, nx, ny);

  eps = 1e-14 * max(xmax - xmin, ymax - ymin);

  // Loop over output grid points; inner work is vectorised over input scatter points
  // For large grids, looping over output pixels is unavoidable without a kd-tree,
  // but the inner distance+sort is array ops so it's reasonably fast.
  for (j = 1; j <= ny; j++) {
    for (i = 1; i <= nx; i++) {
      dx = xi - xg(i, j);
      dy = yi - yg(i, j);
      d2 = dx*dx + dy*dy;

      // Find nneighbors nearest: partial sort via index sort, take first nneighbors
      idx = sort(d2);
      idx = idx(1:nneighbors);

      d2k = d2(idx);
      zk  = zi(idx);

      // Exact hit guard
      if (d2k(1) < eps*eps) {
        result(i, j) = zk(1);
        continue;
      }

      w = 1.0 / d2k^(power * 0.5);
      result(i, j) = sum(w * zk) / sum(w);
    }
  }

  xout = xg;
  yout = yg;
  
  return result;
}
