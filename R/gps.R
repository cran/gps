ExampleBS <- function (d = 4, uniform = TRUE, clamped = FALSE) {
  if (d < 2 || d > 4) stop("d = 2, 3, 4 required!", call. = FALSE)
  if (uniform) {
    lbnd <- as.double(seq.int(to = -1, length.out = 4, by = 2))
    rbnd <- as.double(seq.int(from = 1, length.out = 4, by = 2))
    clamped <- FALSE
  } else if (clamped) {
    lbnd <- rep.int(-1, 4)
    rbnd <- c(1, 2, 2.5, 4)
  } else {
    lbnd <- c(-3, -2.25, -1.75, -1)
    rbnd <- c(1, 2, 2.5, 4)
  }
  xt <- c(lbnd, rbnd)
  if (clamped) {
    x <- MakeGrid(xt[4:8], n = 21)
  } else {
    x <- MakeGrid(xt, n = 21)
  }
  Bsparse <- splines::splineDesign(xt, x, d, outer.ok = TRUE, sparse = TRUE)
  B <- Zero2NA(Bsparse)
  list(xt = xt, clamped = clamped, d = d, x = x, B = B)
}

PlotExampleBS <- function (object) {
  xt <- object$xt
  d <- object$d
  x <- object$x
  B <- object$B
  p <- ncol(B)
  if (object$clamped) {
    ip <- c(rep.int(NA_real_, 4 - d), apply(B[, (4 - d + 1):p], 2, which.max))
  } else {
    ip <- apply(B, 2, which.max)
  }
  xp <- x[ip]
  Bp <- B[cbind(ip, 1:p)]
  ymax <- max(Bp, na.rm = TRUE) * 1.2
  ylim <- c(0, ymax)
  nonzero <- seq.int(5 - d, length.out = d)
  col.Bspl <- rep.int(8, p); col.Bspl[nonzero] <- seq_len(d)
  op <- par(xpd = TRUE)
  on.exit(par(op))
  matplot(x, B, type = "l", lty = 1, lwd = 2, col = col.Bspl, ylim = ylim,
          xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n", ann = FALSE)
  segments(xt[1], 0, xt[8], 0)
  polygon(c(-1, 1, 1, -1), c(0, 0, ymax, ymax), col = gray(0.3, 0.2), border = NA)
  col.knots <- rep.int(8, 8); col.knots[seq.int(5 - d, 4 + d)] <- 1
  if (object$clamped) {
    stepsize <- 0.2
    ystack <- seq.int(to = stepsize, length.out = 3, by = -stepsize)
    points.default(xt[4:8], numeric(5), pch = 19, col = col.knots[4:8], cex = 1.5)
    points.default(xt[1:3], ystack, pch = 19, col = col.knots[1:3], cex = 1.5)
  } else {
    points.default(xt, numeric(8), pch = 19, col = col.knots, cex = 1.5)
  }
  lab.Bspl <- vector("expression", p)
  for (j in 1:p) {
    lab.Bspl[j] <- eval(substitute(expression(B[list(j, d)]), list(j = j, d = d)))
  }
  text.default(xp, Bp, labels = lab.Bspl, col = col.Bspl, adj = c(0.2, 0), cex = 2)
  lab.knots <- vector("expression", 8)
  for (j in 1:8) {
    lab.knots[j] <- eval(substitute(expression(t[j]), list(j = j)))
  }
  if (object$clamped) {
    text.default(xt[1:3], ystack, labels = lab.knots[1:3], col = col.knots[1:3],
                 cex = 2, pos = 4)
    text.default(xt[4:8], numeric(5), labels = lab.knots[4:8], col = col.knots[4:8],
                 cex = 2, pos = 1)
  } else {
    text.default(xt, numeric(8), labels = lab.knots, col = col.knots,
                 cex = 2, pos = 1)
  }
  usr <- par("usr")
  text.default(usr[2], usr[4], labels = c(paste0("order: ", d)),
               cex = 2, adj = c(1, 1))
}

DemoBS <- function (uniform = TRUE, clamped = FALSE) {
  if (uniform && clamped) {
    warning("uniform = TRUE implies clamped = FALSE!")
  }
  spl3 <- ExampleBS(d = 4, uniform, clamped)
  spl2 <- ExampleBS(d = 3, uniform, clamped)
  spl1 <- ExampleBS(d = 2, uniform, clamped)
  if (uniform) {
    caption <- "uniform B-splines"
  } else if (clamped) {
    caption <- "non-uniform & clamped B-splines"
  } else {
    caption <- "non-uniform B-splines"
  }
  op <- par(mfrow = c(3, 1), mar = c(2, 0.75, 0, 0.75), oma = c(0, 0, 2, 0))
  on.exit(par(op))
  PlotExampleBS(spl3)
  PlotExampleBS(spl2)
  PlotExampleBS(spl1)
  title(caption, outer = TRUE, cex.main = 1.5)
}

DemoKnots <- function (aligned = TRUE) {
  xt1 <- seq.int(0, 1, length.out = 12)
  xt2 <- c(0, 0.08, 0.2, 0.27, 0.31, 0.42, 0.54, 0.59, 0.69, 0.82, 0.85, 1)
  xt3 <- rep.int(xt2[4:9], c(4, 1, 1, 1, 1, 4))
  xg1 <- MakeGrid(xt1, n = 21)
  xg2 <- MakeGrid(xt2, n = 21)
  xg3 <- MakeGrid(xt3[4:9], n = 21)
  B1sparse <- splines::splineDesign(xt1, xg1, outer.ok = TRUE, sparse = TRUE)
  B2sparse <- splines::splineDesign(xt2, xg2, outer.ok = TRUE, sparse = TRUE)
  B3sparse <- splines::splineDesign(xt3, xg3, sparse = TRUE)
  B1max <- max(B1sparse@x)
  B2max <- max(B2sparse@x)
  B3max <- max(B3sparse@x)
  B1 <- Zero2NA(B1sparse)
  B2 <- Zero2NA(B2sparse)
  B3 <- Zero2NA(B3sparse)
  Bspl.col <- c(1:6, 8, 7)
  op <- par(mfrow = c(3, 1), mar = c(0.5, 0.5, 1.5, 0.5), xpd = TRUE)
  on.exit(par(op))
  matplot(xg1, B1, ann = FALSE, axes = FALSE, type = "l", col = Bspl.col, bty = "n",
          lty = 1, lwd = 2, xlim = c(0, 1), xaxs = "i", yaxs = "i")
  abline(h = 0, col = 8)
  points.default(xt1, numeric(12), pch = 19, cex = 1.5)
  polygon(xt1[c(4, 9, 9, 4)], c(0, 0, B1max, B1max), col = gray(0.3, 0.2), border = NA)
  title("(a) uniform cubic B-splines", cex.main = 1.5)
  matplot(xg2, B2, ann = FALSE, axes = FALSE, type = "l", col = Bspl.col, bty = "n",
          lty = 1, lwd = 2, xlim = c(0, 1), xaxs = "i", yaxs = "i")
  abline(h = 0, col = 8)
  points.default(xt2, numeric(12), pch = 19, cex = 1.5)
  polygon(xt2[c(4, 9, 9, 4)], c(0, 0, B2max, B2max), col = gray(0.3, 0.2), border = NA)
  title("(b) non-uniform cubic B-splines", cex.main = 1.5)
  xlim3 <- if (aligned) c(0, 1) else xt3[c(4, 9)]
  matplot(xg3, B3, ann = FALSE, axes = FALSE, type = "l", col = Bspl.col, bty = "n",
          lty = 1, lwd = 2, xlim = xlim3, xaxs = "i", yaxs = "i")
  abline(h = 0, col = 8)
  if (aligned) {
    segments(0, 0, xt3[4], 0, col = Bspl.col[1], lwd = 2)
    segments(xt3[9], 0, 1, 0, col = Bspl.col[8], lwd = 2)
    segments(xt3[4], 0, xt3[4], 1, col = Bspl.col[1], lty = 2, lwd = 2)
    segments(xt3[9], 0, xt3[9], 1, col = Bspl.col[8], lty = 2, lwd = 2)
  }
  ystack <- seq.int(0, 0.2, length.out = 4)
  points.default(xt3, c(ystack, numeric(4), ystack), pch = 19, cex = 1.5)
  polygon(xt3[c(4, 9, 9, 4)], c(0, 0, B3max, B3max), col = gray(0.3, 0.2), border = NA)
  title("(c) non-uniform & clamped cubic B-splines", cex.main = 1.5)
}

PlotNull1 <- function (n, k, uniform = FALSE, standard.difference = TRUE) {
  d <- 4
  m <- 2
  x <- qbeta(seq.int(0, 1, length.out = n), 3, 3)
  y <- rnorm(n, mean = x, sd = 0.2 * sd(x))
  xt <- PlaceKnots(x, d, k, uniform = uniform)
  xg <- seq.int(0, 1, length.out = 101)
  B <- splines::splineDesign(xt, x, d, sparse = TRUE)
  Bg <- splines::splineDesign(xt, xg, d, sparse = TRUE)
  if (standard.difference) {
    H <- NullD(seq.int(0, 1, length.out = length(xt)), d, m)
  } else {
    H <- NullD(xt, d, m)
  }
  X <- as_matrix(B %*% H)
  Xg <- as_matrix(Bg %*% H)
  fit <- .lm.fit(X, y)
  yh <- c(Xg %*% fit$coefficients)
  ylim <- range(y, yh)
  plot.default(x, y, col = 8, xlim = range(xt), ylim = ylim, ann = FALSE,
               xaxt = "n", yaxt = "n")
  if (uniform) {
    abline(v = xt, lty = 2, col = 8)
  } else {
    abline(v = xt[seq.int(d, length(xt) - d + 1)], lty = 2, col = 8)
  }
  lines.default(xg, yh, lwd = 2)
  b1 <- rep.int(yh[1], 101)
  b2 <- Xg[, 2]
  b2 <- ((0.5 - yh[1]) / (b2[101] - b2[1])) * (b2 - b2[1]) + yh[1]
  matlines(xg, cbind(b1, b2), type = "l", col = 1, lty = 2, lwd = 2)
}





PlotNull2 <- function (n, k, shape1 = 3, shape2 = 3) {
  d <- 4
  m <- 2
  x <- seq.int(0, 1, length.out = n)
  y <- rnorm(n, mean = x, sd = 0.2 * sd(x))
  if (shape1 < 1 && shape2 < 1) {
    xd <- seq.int(-1/8, 1/8, length.out = k + 2)
    xd <- sign(xd) * abs(xd) ^ (1 / 3) + 0.5
  } else {
    xd <- qbeta(seq.int(0, 1, length.out = k + 2), shape1, shape2)
  }
  xt <- c(numeric(d - 1), xd, rep.int(1, d - 1))
  xg <- seq.int(0, 1, length.out = 101)
  B <- splines::splineDesign(xt, x, d, sparse = TRUE)
  Bg <- splines::splineDesign(xt, xg, d, sparse = TRUE)
  H.ups <- NullD(seq.int(0, 1, length.out = length(xt)), d, m)
  H.gps <- NullD(xt, d, m)
  X.ups <- as_matrix(B %*% H.ups)
  X.gps <- as_matrix(B %*% H.gps)
  Xg.ups <- as_matrix(Bg %*% H.ups)
  Xg.gps <- as_matrix(Bg %*% H.gps)
  b.ups <- .lm.fit(X.ups, y)$coefficients
  b.gps <- .lm.fit(X.gps, y)$coefficients
  ups <- c(Xg.ups %*% b.ups)
  gps <- c(Xg.gps %*% b.gps)
  ylim <- range(y, ups, gps)
  plot.default(x, y, col = 8, ylim = ylim, ann = FALSE, xaxt = "n", yaxt = "n")
  abline(v = xd, lty = 2, col = 8)
  lines.default(xg, ups, lwd = 2)
  lines.default(xg, gps, lwd = 2, lty = 2)
  if (shape1 < 1 && shape2 < 1) {
    title("U-quadratic")
  } else {
    title(sprintf("Beta(%d, %d)", shape1, shape2))
  }
}

DemoNull <- function (n, k, type = 1) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  if (type == 1) {
    par(mfrow = c(1, 3), mar = c(0.25, 0.25, 2.5, 0.25), oma = c(0, 0, 1.5, 0))
    PlotNull1(n, k, uniform = TRUE, standard.difference = TRUE)
    title(c("(A) uniform B-splines", "+ standard difference penalty"))
    PlotNull1(n, k, uniform = FALSE, standard.difference = TRUE)
    title(c("(B) non-uniform B-splines", "+ standard difference penalty"))
    PlotNull1(n, k, uniform = FALSE, standard.difference = FALSE)
    title(c("(C) non-uniform B-splines", "+ general difference penalty [new]"))
    title(expression(paste("fitted P-spline at ", lambda == infinity)), outer = TRUE)
  } else if (type == 2) {
    par(mfrow = c(2, 3), mar = c(0.25, 0.25, 1.5, 0.25), oma = c(0, 0, 1.5, 0))
    PlotNull2(n, k, 3, 5)
    PlotNull2(n, k, 5, 5)
    PlotNull2(n, k, 5, 3)
    PlotNull2(n, k, 1, 3)
    PlotNull2(n, k, 0.5, 0.5)
    PlotNull2(n, k, 3, 1)
    title(expression(paste("fitted naive/general P-spline at ", lambda == infinity)),
          outer = TRUE)
  } else {
    stop("unknown 'type'!")
  }
}

DemoPBS <- function (uniform = TRUE) {
  if (uniform) xd <- seq.int(0, 0.6, by = 0.1)
  else xd <- c(0, 0.75, 1.25, 2, 3.5, 4.5, 5)
  x <- MakeGrid(xd, 21)
  PBsparse <- pbsDesign(x, xd, 4, 0, sparse = TRUE)
  PB <- Zero2NA(PBsparse)
  a <- xd[1L]
  b <- xd[7L]
  period <- b - a
  raux <- xd[2:4] + period
  xt <- c(xd, raux)
  x.wrapped <- x[1:63] + period
  B <- PB
  B[1:63, 4:6] <- NA_real_
  B.wrapped <- PB[1:63, 4:6]
  ip <- apply(PB, 2, which.max)
  xp.PBSpl <- x[ip]
  yp <- PB[cbind(ip, 1:6)]
  xp.PBSpl.modified <- xp.PBSpl[4:6]
  sub <- xp.PBSpl.modified < xd[4]
  xp.PBSpl.modified[sub] <- xp.PBSpl.modified[sub] + period
  xp.BSpl <- c(xp.PBSpl[1:3], xp.PBSpl.modified)
  xlim <- xt[c(1L, 10L)]
  ymax <- max(yp) * 1.15
  pch.knots <- rep.int(c(19, 15, 19, 15), c(1, 3, 3, 3))
  col.knots <- rep.int(c(1, 8), c(7, 3))
  lab1.knots <- vector("expression", 10)
  lab1.knots[1] <- expression(a)
  for (j in 1:5) {
    lab1.knots[j + 1] <- eval(substitute(expression(s[j]), list(j = j)))
  }
  lab1.knots[7] <- expression(b)
  for (j in 1:3) {
    lab1.knots[j + 7] <- eval(substitute(expression(s[j]^','), list(j = j)))
  }
  lab2.knots <- vector("expression", 10)
  for (j in 1:10) {
    lab2.knots[j] <- eval(substitute(expression(xi[j]), list(j = j)))
  }
  lab.BSpl <- vector("expression", 6)
  for (j in 1:6) {
    lab.BSpl[j] <- eval(substitute(expression(B[j]), list(j = j)))
  }
  op <- par(mfrow = c(2, 1), xpd = TRUE, mar = c(3, 1, 1.25, 1))
  on.exit(par(op))
  matplot(x, PB, type = "n", ann = FALSE, xlim = xlim, ylim = c(0, ymax), 
          xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
  segments(xt[1], 0, xt[10], 0, col = 8)
  polygon(c(a, b, b, a), c(0, 0, ymax, ymax), col = gray(0.3, 0.2), border = NA)
  matlines(x, B, lty = rep.int(c(1, 2), c(3, 3)), lwd = 2, col = 1:6)
  matlines(x.wrapped, B.wrapped, lty = 2, lwd = 2, col = 4:6)
  text.default(xp.BSpl, yp, labels = lab.BSpl, col = 1:6, adj = c(0.2, 0), cex = 1.5)
  points.default(xt, numeric(10), pch = pch.knots, col = col.knots)
  text.default(xt, numeric(10), labels = lab1.knots, col = col.knots, pos = 1, cex = 1.5)
  mtext(lab2.knots, side = 1, line = 2, at = xt, col = col.knots, cex = 1.5)
  title("ordinary cubic B-splines")
  matplot(x, PB, type = "n", ann = FALSE, xlim = xlim, ylim = c(0, ymax),
          xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
  segments(a, 0, b, 0, col = 8)
  polygon(c(a, b, b, a), c(0, 0, ymax, ymax), col = gray(0.3, 0.2), border = NA)
  matlines(x, B, col = 1:6, lty = rep.int(c(1, 2), c(3, 3)), lwd = 2)
  matlines(x.wrapped - period, B.wrapped, lty = 2, lwd = 2, col = 4:6)
  text.default(xp.PBSpl, yp, labels = lab.BSpl, col = 1:6, adj = c(0.2, 0), cex = 1.5)
  points.default(xd, numeric(7), pch = 19)
  text.default(xd, numeric(7), labels = lab1.knots[1:7], col = col.knots[1:7],
               pos = 1, cex = 1.5)
  mtext(lab2.knots[1:7], side = 1, line = 2, at = xd, col = col.knots[1:7], cex = 1.5)
  title("periodic cubic B-splines")
}


Diff <- function (n, x, xi = 1L, k = 1L) {
  if (!is.double(x)) stop("'x' is not in double-precision mode!")
  .Call("C_Diff", n, k, x, xi, PACKAGE = "gps")
}

SparseWtDelta <- function (h) {
  r <- length(h)
  x <- rep.int(c(-1, 1), r) * rep(1 / h, each = 2)
  i <- rep(seq.int(0L, r - 1L), each = 2)
  p <- c(0L, seq.int(1L, length.out = r, by = 2L), 2L * r)
  methods::new("dgCMatrix", i = i, p = p, Dim = c(r, r + 1L), x = x)
}

SparseD <- function (xt, d) {
  xt <- as.double(xt)
  K <- length(xt)
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  if (K < 2 * d) stop("length(xt) >= 2 * d required!", call. = FALSE)
  D <- vector("list", d - 1)
  h <- Diff(n = K - 2, x = xt, xi = 2, k = d - 1)
  D[[1]] <- SparseWtDelta(h)
  i <- 2
  while (i < d) {
    h <- Diff(n = K - 2 * i, x = xt, xi = i + 1, k = d - i)
    D[[i]] <- SparseWtDelta(h) %*% D[[i - 1]]
    i <- i + 1
  }
  D
}

ComputeLD <- function (xt, d) {
  xt <- as.double(xt)
  K <- length(xt)
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  if (K < 2 * d) stop("length(xt) >= 2 * d required!", call. = FALSE)
  .Call("C_ComputeLD", xt, d, PACKAGE = "gps")
}

OrthNullD <- function (ld, m = 1, orthogonal = TRUE) {
  if (m < 1 || m > ncol(ld)) stop("m = 1:(d - 1) required!", call. = FALSE)
  basis <- .Call("C_NullD", ld, m, PACKAGE = "gps")
  if (orthogonal && (m > 1)) {
    Q <- qr.Q(qr.default(basis[, m:1]))[, m:1]
    i <- sequence.default(1:(m - 1))
    j <- rep.int(2:m, 1:(m - 1))
    Q[cbind(i, j)] <- 0
    basis <- Q
  }
  basis
}

NullD <- function (xt, d, m) {
  ld <- ComputeLD(xt, d)
  OrthNullD(ld, m, orthogonal = FALSE)
}

DiffCoef <- function (b, xt, d, m) {
  K <- length(xt)
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  if (K < 2 * d) stop("length(xt) >= 2 * d required!", call. = FALSE)
  if (m < 0 || m >= d) stop("m = 0:(d - 1) required!", call. = FALSE)
  if (length(b) != K - d) stop("length(b) == length(xt) - d required!", call. = FALSE)
  if (m == 0) return(b)
  for (i in 1:m) {
    h <- Diff(n = K - 2 * i, x = xt, xi = i + 1, k = d - i)
    b <- Diff(length(b), b) / h
  }
  b
}

MakeGrid <- function (xd, n, rm.dup = FALSE) {
  if (is.unsorted(xd, strictly = TRUE)) {
    stop("'xd' is not strictly ascending!")
  }
  xd <- as.double(xd)
  if (n == 1) {
    lp <- xd[-length(xd)]
    rp <- xd[-1]
    mp <- 0.5 * (lp + rp)
    return(mp)
  }
  if (rm.dup && (n == 2)) return(xd)
  .Call("C_MakeGrid", xd, n, rm.dup, PACKAGE = "gps")
}


pbsDesign <- function (x, xd, d, nDeriv = 0, sparse = FALSE) {
  d <- as.integer(d)
  if (d < 2L) stop("d >= 2 required!")
  degree <- d - 1L
  k2 <- length(xd)
  if (k2 <= d) stop("length(xd) >= d + 1 required!")
  a <- xd[1L]
  b <- xd[k2]
  period <- b - a
  if (min(x) < a || max(x) > b) stop("domain does not contain all x-values!")
  p <- k2 - 1L
  raux <- xd[2:d] + period
  xt <- c(xd, raux)
  ind.x.wrapped <- which(x < xd[d])
  if (length(ind.x.wrapped) == 0L) {
    PB <- splines::splineDesign(xt, x, d, nDeriv, TRUE, sparse)
  } else {
    nx <- length(x)
    x.wrapped <- x[ind.x.wrapped] + period
    x.padded <- c(x, x.wrapped)
    Aknots <- c(rep.int(xt[1L], degree), xt, rep.int(xt[length(xt)], degree))
    B <- splines::splineDesign(Aknots, x.padded, d, nDeriv, sparse = TRUE)
    B.p <- B@p[seq.int(d, length(B@p) - degree)]
    PB.p <- B.p - B.p[1L]
    nnz <- PB.p[p + 1L]
    PB.i <- integer(nnz)
    PB.x <- numeric(nnz)
    num_unchanged_bsplines <- p - degree
    read_ind <- seq.int(B.p[1L] + 1L, B.p[num_unchanged_bsplines + 1L])
    num_unchanged_entries <- length(read_ind)
    write_ind <- seq_len(num_unchanged_entries)
    PB.i[write_ind] <- B@i[read_ind]
    PB.x[write_ind] <- B@x[read_ind]
    B.p.vital <- B.p[seq.int(num_unchanged_bsplines + 1L, p + 1L)]
    read_ind <- seq.int(B.p.vital[1L] + 1L, B.p.vital[d])
    unordered.i <- B@i[read_ind]
    unordered.x <- B@x[read_ind]
    sub <- (unordered.i >= nx)
    unordered.i[sub] <- ind.x.wrapped[unordered.i[sub] - (nx - 1L)] - 1L
    unordered.j <- rep.int(1:degree, diff.default(B.p.vital))
    reorder_ind <- order(unordered.j, unordered.i)
    write_ind <- seq.int(num_unchanged_entries + 1L, nnz)
    PB.i[write_ind] <- unordered.i[reorder_ind]
    PB.x[write_ind] <- unordered.x[reorder_ind]
    if (sparse) {
      PB <- methods::new("dgCMatrix", i = PB.i, p = PB.p, Dim = c(nx, p), x = PB.x)
    } else {
      PB.j <- rep.int(1:p, diff.default(PB.p))
      PB <- matrix(0, nx, p)
      PB[cbind(PB.i + 1L, PB.j)] <- PB.x
    }
  }
  PB
}

SparsePD <- function (xd, d) {
  d <- as.integer(d)
  if (d < 2L) stop("d >= 2 required!")
  degree <- d - 1L
  k2 <- length(xd)
  if (k2 <= d) stop("length(xd) >= d + 1 required!")
  a <- xd[1L]
  b <- xd[k2]
  period <- b - a
  laux <- xd[(k2 - degree):(k2 - 1L)] - period
  raux <- xd[2:d] + period
  xt <- c(laux, xd, raux)
  p <- k2 - 1L
  D <- SparseD(xt, d)
  PD <- vector("list", degree)
  for (m in 1:degree) {
    D.i <- D[[m]]@i
    D.p <- D[[m]]@p
    D.x <- D[[m]]@x
    m1 <- m + 1L
    nnz <- p * m1
    ind <- integer(nnz)
    done <- 0L
    for (j in 1:m) {
      ind1 <- seq.int(D.p[j] + 1L, D.p[j + 1L])
      ind2 <- seq.int(D.p[p + j] + 1L, length.out = m1 - j)
      ind[seq.int(done + 1L, length.out = m1)] <- c(ind1, ind2)
      done <- done + m1
    }
    ind[(done + 1L):nnz] <- seq.int(D.p[m1] + 1L, D.p[p + 1L])
    PD.i <- D.i[ind]
    PD.p <- seq.int(0L, by = m1, length.out = p + 1L)
    PD.x <- D.x[ind]
    PD[[m]] <- methods::new("dgCMatrix", i = PD.i, p = PD.p, Dim = c(p, p), x = PD.x)
  }
  PD
}

PlaceKnots <- function (x, d, k, domain = NULL, uniform = TRUE, periodic = FALSE) {
  xu <- sort(unique.default(x))
  nx <- length(xu)
  if (is.null(domain)) {
    a <- xu[1L]
    b <- xu[nx]
    domain <- c(a, b)
  } else {
    a <- domain[1L]
    b <- domain[2L]
    if (xu[1L] < a || xu[nx] > b) {
      stop("'domain' does not contain all x-values!")
    }
  }
  degree <- d - 1
  if (uniform) {
    xd <- seq.int(a, b, length.out = k + 2)
    h <- xd[2] - xd[1]
    laux <- seq.int(to = a - h, by = h, length.out = degree)
    raux <- seq.int(from = b + h, by = h, length.out = degree)
  } else {
    prob <- seq.int(0, 1, length.out = k + 2)
    xd <- quantile(xu, prob, names = FALSE)
    xd[c(1, k + 2)] <- domain
    laux <- rep.int(a, degree)
    raux <- rep.int(b, degree)
  }
  if (periodic) xd else c(laux, xd, raux)
}

QuadWts <- function (ord) {
  if (ord == 1) {
    return(2)
  }
  p <- seq.int(0, ord - 1)
  P <- outer(seq.int(-1, 1, length.out = ord), p, "^")
  Pinv <- solve.default(P)
  pow <- outer(1:ord, p, "+")
  H <- (1 - (-1) ^ pow) / pow
  base::crossprod(Pinv, H %*% Pinv)
}

SbarBlocks <- function (xt, d, m) {
  xt <- as.double(xt)
  K <- length(xt)
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  if (K < 2 * d) stop("length(xt) >= 2 * d required!", call. = FALSE)
  if (m < 0 || m >= d) stop("m = 0:(d - 1) required!", call. = FALSE)
  ord <- d - m
  if (ord == 1) {
    h <- Diff(K - 2 * (d - 1), xt, xi = d)
    return(h)
  }
  W <- QuadWts(ord)
  xd <- xt[seq.int(d, K - d + 1)]
  xg <- MakeGrid(xd, ord)
  xt.local <- xt[seq.int(1 + m, K - m)]
  B <- splines::splineDesign(xt.local, xg, ord, sparse = TRUE)
  .Call("C_SbarBlocks", xd, W, B@x, PACKAGE = "gps")
}

SparseSbar <- function (xt, d, m, do.chol = FALSE) {
  blocks <- SbarBlocks(xt, d, m)
  if (d - m == 1) {
    if (do.chol) blocks <- sqrt(blocks)
    return(blocks)
  }
  LTB <- .Call("C_SbarLTB", blocks, do.chol, PACKAGE = "gps")
  Matrix::bandSparse(n = ncol(LTB), k = seq.int(0, nrow(LTB) - 1),
                     diagonals = t.default(LTB), symmetric = !do.chol)
}

SparseS <- function (xt, d, m, root = FALSE) {
  xt <- as.double(xt)
  if (m == 0) return(SparseSbar(xt, d, 0, do.chol = root))
  U <- SparseSbar(xt, d, m, do.chol = TRUE)
  D <- SparseD(xt, d)[[m]]
  if (d - m == 1) {
    K <- D
    K@x <- U[D@i + 1L] * D@x
  } else {
    K <- U %*% D
  }
  if (root) K else Matrix::crossprod(K)
}

btSb <- function (b, xt, d, m) {
  db <- DiffCoef(b, xt, d, m)
  blocks <- SbarBlocks(xt, d, m)
  .Call("C_btSb", blocks, db, PACKAGE = "gps")
}

Zero2NA <- function (Bsparse) {
  if (!inherits(Bsparse, "dgCMatrix")) {
    stop("'Bsparse' is not a \"dgCMatrix\"!")
  }
  B <- matrix(NA_real_, Bsparse@Dim[1], Bsparse@Dim[2])
  i <- Bsparse@i + 1L
  j <- rep.int(seq_len(Bsparse@Dim[2]), diff.default(Bsparse@p))
  B[cbind(i, j)] <- Bsparse@x
  B
}

as_matrix <- function (A) {
  if (is.matrix(A)) return(A)
  sparse <- inherits(A, "dsparseMatrix")
  dense <- inherits(A, "ddenseMatrix")
  if (!sparse && !dense) {
    stop("'A' is not a \"dsparseMatrix\" or \"ddenseMatrix\"!")
  }
  nnz <- length(A@x)
  nr <- A@Dim[1]
  nc <- A@Dim[2]
  if (nnz == nr * nc) {
    denseA <- matrix(A@x, nr, nc)
  } else if (inherits(A, "dCsparseMatrix")) {
    i <- A@i
    j <- rep.int(seq.int(0L, nc - 1L), diff.default(A@p))
    denseA <- matrix(0, nr, nc)
    denseA[j * nr + (i + 1L)] <- A@x
    if (inherits(A, "dsCMatrix")) denseA[i * nc + (j + 1L)] <- A@x
  } else {
    stop("Not implemented yet!")
  }
  denseA
}

