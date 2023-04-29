Csparse2LTB <- function (A) {
  if (!inherits(A, "dCsparseMatrix")) {
    stop("'A' is not a \"dCsparseMatrix\"!")
  }
  n <- A@Dim[1L]
  k <- A@Dim[2L]
  i <- A@i
  j <- rep.int(1:k, diff.default(A@p))
  b1 <- sum(i == 0L)
  b <- b1 - 1L
  if (k < n || k > n + b) {
    stop(sprintf("nrow(A) <= ncol(A) <= nrow(A) + %d expected!", b))
  }
  l <- j - i
  if (any(l > b1 | l < 1L)) {
    stop("Nonzero entries detected outside the indended band region!")
  }
  N1 <- sum(seq.int(n, by = -1L, length.out = b1))
  N2 <- sum(seq.int(b, by = -1L, length.out = k - n))
  nnz.max <- N1 + N2
  x <- A@x
  nnz <- length(x)
  if (nnz == nnz.max) {
    .Call("C_Csparse2LTB", b1, n, k, x, PACKAGE = "gps")
  } else {
    At <- matrix(0, b1, k)
    At[cbind(l, i + 1L)] <- A@x
    At
  }
}

LTB2Csparse <- function (L, k = n, symmetric = FALSE) {
  b1 <- nrow(L); n <- ncol(L); b <- b1 - 1L
  k <- as.integer(k)
  if (symmetric) {
    k <- n
  } else if (k < n || k > n + b) {
    stop(sprintf("%d <= k <= %d expected!", n, n + b))
  }
  r <- n + b - k
  p <- rep.int(seq.int(b1, by = -1L, length.out = r + 1L), c(n - r, rep.int(1, r)))
  i <- sequence.default(nvec = p, from = seq.int(0L, n - 1L))
  p <- c(0L, cumsum(p))
  if (r) {
    x <- numeric(p[length(p)])
    x[seq_len(b1 * (k - b))] <- L[, seq_len(k - b)]
    columns <- seq.int(to = n, length.out = r)
    for (j in columns) {
      ind <- seq.int(p[j] + 1L, p[j + 1L])
      x[ind] <- L[seq_len(length(ind)), j]
    }
  } else {
    x <- c(L)
  }
  if (symmetric) {
    methods::new("dsCMatrix", i = i, p = p, x = x, Dim = c(k, n), uplo = "L")
  } else {
    methods::new("dgCMatrix", i = i, p = p, x = x, Dim = c(k, n))
  }
}

LTB2Dense <- function (L, k) {
  .Call("C_LTB2Dense", L, k, PACKAGE = "gps")
}

SolveLTB <- function (transA = FALSE, A, b, overwrite = FALSE) {
  .Call("C_SolveLTB", transA, A, b, overwrite, PACKAGE = "gps")
}

LPBTRF <- function (A, overwrite = FALSE) {
  .Call("C_LPBTRF", A, overwrite, PACKAGE = "gps")
}

DDt <- function (Dt) {
  .Call("C_DDt", Dt, PACKAGE = "gps")
}

DtD <- function (Dt) {
  .Call("C_DtD", Dt, PACKAGE = "gps")
}

IsMonoInc <- function (x, n = length(x), xi = 1L) {
  .Call("C_IsMonoInc", x, n, xi, PACKAGE = "gps") > 0L
}

MonoKnots <- function (xt, d) {
  xt <- as.double(xt)
  K <- length(xt)
  if (K < 2 * d) stop("length(xt) >= 2 * d required!", call. = FALSE)
  flag <- IsMonoInc(xt, n = K - 2 * (d - 1), xi = d)
  if (!flag) stop("Domain knots are not strictly ascending!", call. = FALSE)
  lbnd <- xt[seq.int(from = 1L, length.out = d)]
  rbnd <- xt[seq.int(to = K, length.out = d)]
  if (is.unsorted(lbnd) || is.unsorted(rbnd)) {
    stop("Boundary knots are not non-decreasing!", call. = FALSE)
  }
  xt
}

PlaceKnots <- function (x, d, k, domain = NULL, uniform = FALSE, periodic = FALSE) {
  xu <- sort.int(unique.default(x))
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
    h <- xd[2L] - xd[1L]
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

MakeGrid <- function (xd, n, rm.dup = FALSE) {
  xd <- as.double(xd)
  if (!IsMonoInc(xd)) stop("'xd' is not strictly ascending!")
  if (n == 1) {
    lp <- xd[-length(xd)]
    rp <- xd[-1]
    mp <- 0.5 * (lp + rp)
    return(mp)
  }
  if (rm.dup && (n == 2)) return(xd)
  .Call("C_MakeGrid", xd, n, rm.dup, PACKAGE = "gps")
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

ExBS <- function (d = 4, uniform = TRUE, clamped = FALSE) {
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

PlotExBS <- function (object) {
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
  spl3 <- ExBS(d = 4, uniform, clamped)
  spl2 <- ExBS(d = 3, uniform, clamped)
  spl1 <- ExBS(d = 2, uniform, clamped)
  if (uniform) {
    caption <- "uniform B-splines"
  } else if (clamped) {
    caption <- "non-uniform & clamped B-splines"
  } else {
    caption <- "non-uniform B-splines"
  }
  op <- par(mfrow = c(3, 1), mar = c(2, 0.75, 0, 0.75), oma = c(0, 0, 2, 0))
  on.exit(par(op))
  PlotExBS(spl3)
  PlotExBS(spl2)
  PlotExBS(spl1)
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
  ystack <- seq.int(0, by = 0.1, length.out = 4)
  points.default(xt3, c(ystack, numeric(4), ystack), pch = 19, cex = 1.5)
  polygon(xt3[c(4, 9, 9, 4)], c(0, 0, B3max, B3max), col = gray(0.3, 0.2), border = NA)
  title("(c) non-uniform & clamped cubic B-splines", cex.main = 1.5)
}





PlotNull <- function (x, y, k, shape1 = 3, shape2 = 3, gps = FALSE) {
  d <- 4
  m <- 2
  if (shape1 < 1 && shape2 < 1) {
    xd <- seq.int(-1/8, 1/8, length.out = k + 2)
    xd <- sign(xd) * abs(xd) ^ (1 / 3) + 0.5
    distr <- "U-quadratic(0, 1)"
  } else {
    xd <- qbeta(seq.int(0, 1, length.out = k + 2), shape1, shape2)
    distr <- sprintf("Beta(%d, %d)", shape1, shape2)
  }
  xt <- c(numeric(d - 1), xd, rep.int(1, d - 1))
  K <- length(xt)
  xg <- seq.int(0, 1, length.out = 101)
  B <- splines::splineDesign(xt, x, d, sparse = TRUE)
  Bg <- splines::splineDesign(xt, xg, d, sparse = TRUE)
  if (gps) {
    ld <- ComputeLD(xt, d)
  } else {
    ld <- matrix(1, k + d, d - 1)
  }
  H <- NullGD(ld, m, orthonormal = FALSE)
  X <- as_matrix(B %*% H)
  Xg <- as_matrix(Bg %*% H)
  XtX <- base::crossprod(X)
  Xty <- base::crossprod(X, y)
  U <- chol.default(XtX)
  f <- forwardsolve(U, Xty, upper.tri = TRUE, transpose = TRUE)
  b <- backsolve(U, f)
  yg <- c(Xg %*% b)
  ylim <- range(y, yg)
  plot.default(x, y, col = 8, ylim = ylim, ann = FALSE, xaxt = "n", yaxt = "n")
  abline(v = xd, lty = 2, col = 8)
  lines.default(xg, yg, lwd = 2)
  title(distr)
}

DemoNull <- function (n, k, gps = FALSE) {
  x <- seq.int(0, 1, length.out = n)
  y <- rnorm(n, mean = x, sd = 0.2 * sd(x))
  op <- par(mfrow = c(2, 3), mar = c(0.25, 0.25, 1.5, 0.25), oma = c(0, 0, 1.5, 0))
  on.exit(par(op))
  PlotNull(x, y, k, 3, 5, gps)
  PlotNull(x, y, k, 5, 5, gps)
  PlotNull(x, y, k, 5, 3, gps)
  PlotNull(x, y, k, 1, 3, gps)
  PlotNull(x, y, k, 0.5, 0.5, gps)
  PlotNull(x, y, k, 3, 1, gps)
  label <- sprintf("limiting %s P-spline fit at ", c("standard", "general")[gps + 1L])
  expr <- substitute(expression(paste(label, lambda == +infinity)), list(label = label))
  title(eval(expr), outer = TRUE)
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

GetE_periodic <- function (xd, d) {
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  degree <- d - 1
  k2 <- length(xd)
  if (k2 <= d) stop("length(xd) >= d + 1 required!", call. = FALSE)
  a <- rep.int(xd[1], degree)
  b <- rep.int(xd[k2], degree)
  local.knots.a <- c(a, xd[1:(d + 1)])
  local.knots.b <- c(xd[(k2 - d):k2], b)
  qs <- 0:(degree - 1)
  Ca <- splines::splineDesign(local.knots.a, a, d, derivs = qs)
  Cb <- splines::splineDesign(local.knots.b, b, d, derivs = qs)
  Ca <- Ca[, 1:degree]
  Cb <- Cb[, 2:d]
  det.a <- prod(Ca[seq.int(1, length.out = degree, by = d)])
  det.b <- abs(prod(Cb[seq.int(degree, length.out = degree, by = degree - 1)]))
  if (det.a >= det.b) {
    E <- forwardsolve(Ca, Cb)
    attr(E, "eliminate") <- "head"
  } else {
    perm <- degree:1
    E <- forwardsolve(Cb[, perm], Ca)[perm, ]
    attr(E, "eliminate") <- "tail"
  }
  E
}

AbsorbE_periodic <- function (X, E) {
  p <- ncol(X)
  degree <- nrow(E)
  X.head <- X[, 1:degree]
  X.tail <- X[, (p - degree + 1):p]
  X.mid <- X[, (degree + 1):(p - degree)]
  if (attr(E, "eliminate") == "head") {
    PX <- cbind(X.mid, X.tail + X.head %*% E)
  } else {
    PX <- cbind(X.head + X.tail %*% E, X.mid)
  }
  PX
}

DemoPBS2 <- function (uniform = TRUE, scaling = TRUE) {
  if (uniform) xd <- seq.int(0, 0.6, by = 0.1)
  else xd <- c(0, 0.75, 1.25, 2, 3.5, 4.5, 5)
  a <- xd[1]
  b <- xd[7]
  Aknots <- c(rep.int(a, 3), xd, rep.int(b, 3))
  x <- gps::MakeGrid(xd, n = 21)
  B <- splines::splineDesign(Aknots, x)
  E <- GetE_periodic(xd, 4)
  PB <- AbsorbE_periodic(B, E)
  lab.knots <- vector("expression", 7)
  for (j in 1:7) {
    lab.knots[j] <- eval(substitute(expression(xi[j]), list(j = j)))
  }
  lab.BSpl <- vector("expression", 9)
  for (j in 1:9) {
    lab.BSpl[j] <- eval(substitute(expression(B[j]), list(j = j)))
  }
  lab.PBSpl <- vector("expression", 9)
  for (j in 1:6) {
    lab.PBSpl[j] <- eval(substitute(expression(tilde(B)[j]), list(j = j)))
  }
  ip.BSpl <- apply(B, 2, which.max)
  xp.BSpl <- x[ip.BSpl]
  yp.BSpl <- B[cbind(ip.BSpl, 1:9)]
  ip.PBSpl <- apply(abs(PB), 2, which.max)
  xp.PBSpl <- x[ip.PBSpl]
  yp.PBSpl <- PB[cbind(ip.PBSpl, 1:6)]
  positive <- yp.PBSpl > 0
  negative <- yp.PBSpl < 0
  if (scaling) {
    big <- abs(yp.PBSpl) > 1
    big.yp.PBspl <- yp.PBSpl[big]
    fctr <- sign(big.yp.PBspl) / max(abs(big.yp.PBspl))
    PB[, big] <- PB[, big] * rep(fctr, each = length(x))
    yp.PBSpl[big] <- yp.PBSpl[big] * fctr
  }
  B[B == 0] <- NA
  PB[PB == 0] <- NA
  op <- par(mfrow = c(2, 1), xpd = TRUE, mar = c(2, 0.75, 1, 0.75))
  on.exit(par(op))
  matplot(x, B, type = "l", lty = rep(c(2, 1, 2), each = 3), lwd = 2, ann = FALSE,
          xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n", ylim = c(0, 1.1))
  segments(a, 0, b, 0, col = 8)
  text.default(xp.BSpl, yp.BSpl, labels = lab.BSpl, col = 1:6, adj = c(0.5, 0), cex = 1.5)
  points.default(xd, numeric(7), pch = 19)
  ystack <- seq.int(0.1, by = 0.1, length.out = 3)
  points.default(rep.int(c(a, b), c(3, 3)), c(ystack, ystack), pch = 19, col = 8)
  text.default(xd, numeric(7), labels = lab.knots, pos = 1, cex = 1.5)
  title("ordinary cubic B-splines")
  if (attr(E, "eliminate") == "head") {
    lty.PBSpl <- rep(1:2, each = 3)
  } else {
    lty.PBSpl <- rep(2:1, each = 3)
  }
  ylim.PBSpl <- range(-0.03, yp.PBSpl) * c(1, 1.1)
  matplot(x, PB, type = "l", lty = lty.PBSpl, lwd = 2, ann = FALSE, ylim = ylim.PBSpl,
          xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
  segments(a, 0, b, 0, col = 8)
  text.default(xp.PBSpl[positive], yp.PBSpl[positive], labels = lab.PBSpl[positive],
       col = (1:6)[positive], adj = c(0.5, 0), cex = 1.5)
  text.default(xp.PBSpl[negative], yp.PBSpl[negative], labels = lab.PBSpl[negative],
       col = (1:6)[negative], pos = 1, cex = 1.5)
  points.default(xd, numeric(7), pch = 19)
  text.default(xd, numeric(7), labels = lab.knots, pos = 1, cex = 1.5)
  title("periodic cubic B-splines")
}

DemoSpl <- function (uniform = TRUE) {
  if (uniform) {
    xd <- 1:6
    xt <- seq.int(-2, 9)
    xg <- MakeGrid(xt, 21)
    b <- c(0.44, 1.11, 1.66, 0.25, 1.6, 1.43, 1.49, 2.52)
  } else {
    xd <- c(1, 2.06, 2.65, 3.22, 3.91, 6)
    xt <- xd[c(1, 1, 1, 1, 2:5, 6, 6, 6, 6)]
    xg <- MakeGrid(xd, 21)
    b <- c(1.9, 1.84, 1.35, 1.67, 0.48, 1.85, 1.59, 1.35)
  }
  Bg <- splines::splineDesign(xt, xg, outer.ok = uniform, sparse = TRUE)
  Bgdense <- Zero2NA(Bg)
  yg <- (Bg %*% b)@x
  if (uniform) {
    i1 <- 21 * 3 + 1
    i2 <- 21 * 8
    x <- xg[i1:i2]
    y <- yg[i1:i2]
  } else {
    x <- xg
    y <- yg
  }
  ymax <- max(y)
  yd <- y[c(1, 1:5 * 21)]
  op <- par(xpd = TRUE, mar = c(2, 2, 1.5, 0.5))
  on.exit(par(op))
  plot.default(x, y, type = "n", xlim = c(-2, 9), ylim = c(0, ymax),
       xaxs = "i", yaxs = "i", ann = FALSE, bty = "n")
  polygon(c(1, 6, 6, 1), c(0, 0, ymax, ymax), col = gray(0.3, 0.2), border = NA)
  lines.default(x, y, lwd = 2)
  points.default(xd, yd, pch = c(17, 19, 19, 19, 19, 17))
  matlines(xg, Bgdense, type = "l", lty = 1, lwd = 2, col = c(1:6, 8, 7))
  points.default(xd[2:5], numeric(4), pch = 19)
  if (uniform) {
    points.default(c(-2:0, 7:9), numeric(6), pch = 15)
  } else {
    segments(-2, 0, 1, 0, lwd = 2, col = 1)
    segments(6, 0, 9, 0, lwd = 2, col = 7)
    segments(1, 1, 1, 0, lwd = 2, col = 1, lty = 2)
    segments(6, 1, 6, 0, lwd = 2, col = 7, lty = 2)
    ystack <- seq.int(0.1, by = 0.1, length.out = 3)
    points.default(c(1, 1, 1, 6, 6, 6), c(ystack, ystack), pch = 15)
  }
  points.default(c(1, 6), numeric(2), pch = 17)
  type <- c("non-uniform", "uniform")[uniform + 1L]
  title(sprintf("cubic spline with %s knots", type))
  dB <- splines::splineDesign(xt, rep.int(xd[1:5], 3), derivs = rep(1:3, each = 5))
  pc <- c(dB %*% b) * rep(1 / factorial(1:3), each = 5)
  lab1 <- paste0("f", 1:5)
  lab2 <- c("constant", "linear", "quadratic", "cubic")
  pc <- matrix(c(yd[1:5], pc), ncol = 4, dimnames = list(lab1, lab2))
  list(xd = xd, bspl.coef = b, piecepoly.coef = pc)
}

CheckSupp <- function (x, xt, d) {
  degree <- d - 1
  K <- length(xt)
  xd <- xt[seq.int(d, K - degree)]
  nI <- length(xd) - 1L
  id <- findInterval(x, xd, rightmost.closed = TRUE)
  if (min(id) < 1L || max(id) > nI) {
    stop("Domain does not contain all x-values!", call. = FALSE)
  }
  ni <- tabulate(id, nbins = nI)
  .Call("C_CheckSupp", ni, d, PACKAGE = "gps") > 0L
}

gps2bspl <- function (x, xt, d, m, gps = TRUE, periodic = FALSE,
                      unique.x.provided = FALSE) {
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  if (m < 1 || m > d - 1) stop("1 <= m <= d - 1 required!", call. = FALSE)
  xt <- MonoKnots(xt, d)
  p <- length(xt) - d
  xu <- if (unique.x.provided) x else unique.default(x)
  if (p > length(xu)) {
    stop("There are more B-splines than unique x-values!", call. = FALSE)
  }
  if (CheckSupp(xu, xt, d)) {
    stop("There are no x-values on some B-spline's support!", call. = FALSE)
  }
  B <- splines::splineDesign(xt, x, d, sparse = TRUE)
  D <- SparseD(xt, d, m, gps)[[1L]]
  Dt <- Csparse2LTB(D)
  S <- DtD(Dt)
  if (ncol(S) > p) S <- S[, 1:p]
  if (periodic) {
    C <- GetPBC(xt, d, compact = FALSE, sparse = FALSE)
  } else {
    C <- matrix(0, nrow = 0, ncol = p)
  }
  list(B = B, Dt = Dt, S = S, xt = xt, d = d, m = m, C = C)
}

gps2red <- function (bspl, y, w = NULL, scalePen = TRUE) {
  y <- as.double(y)
  X <- bspl$B
  n <- X@Dim[1L]
  if (length(y) != n) stop("Incorrect number of y-values!")
  if (!is.null(w)) {
    if (length(w) != n) stop("Incorrect number of weights!")
    w <- (1 / sum(w)) * w
    w <- sqrt(w)
    X@x <- w[X@i + 1L] * X@x
    y <- w * y
  }
  if (scalePen) {
    scalePen.fctr <- VecDot(X@x) / VecDot(bspl$Dt)
    Dt <- VecScal(sqrt(scalePen.fctr), bspl$Dt, overwrite = FALSE)
    S <- VecScal(scalePen.fctr, bspl$S, overwrite = FALSE)
  } else {
    scalePen.fctr <- 1
    Dt <- bspl$Dt
    S <- bspl$S
  }
  Z <- Matrix::crossprod(X)
  z <- Matrix::crossprod(X, y)@x
  Z <- Csparse2LTB(Z)
  L <- try(LPBTRF(Z, overwrite = FALSE), silent = TRUE)
  if (inherits(L, "try-error")) {
    stop("B-spline design matrix does not have full column rank!")
  }
  f <- SolveLTB(transA = FALSE, L, z, overwrite = FALSE)
  minRSS <- VecDot(y) - VecDot(f)
  E <- .Call("C_FormE", L, Dt, PACKAGE = "gps")
  list(n = n, Z = Z, z = z, L = L, f = f, minRSS = minRSS,
       scalePen.fctr = scalePen.fctr, Dt = Dt, S = S, E = E)
}

gps2DR <- function (red, tol = 1e-6) {
  t0 <- Sys.time()
  q <- ncol(red$E)
  dbar <- .Call("C_MeanDR", red$E, PACKAGE = "gps")
  dq <- .Call("C_MinDR", red$E, tol, PACKAGE = "gps")
  d1 <- .Call("C_MaxDR", red$L, red$Dt, tol, PACKAGE = "gps")
  dq.cutoff <- d1 * (.Machine$double.eps / 2)
  if (dq < dq.cutoff) {
    dq <- dq.cutoff
    warning("The smallest eigenvalue is erratic as E'E is numerically singular!")
  }
  approx.DR <- .Call("C_ApproxDR", q, dq, d1, dbar, tol, PACKAGE = "gps")
  t1 <- Sys.time()
  secs <- difftime(t1, t0, units = "secs")[[1L]]
  list(mean = dbar, min = dq, max = d1, approx.values = approx.DR, secs = secs)
}

ExactDR <- function (E) {
  t0 <- Sys.time()
  A <- .Call("C_LAUUM", E, PACKAGE = "gps")
  exact.DR <- base::eigen(A, symmetric = TRUE, only.values = TRUE)$values
  q <- length(exact.DR)
  d1 <- exact.DR[1L]
  dq.cutoff <- d1 * (.Machine$double.eps / 2)
  ind <- which(exact.DR < dq.cutoff)
  if (length(ind)) {
    exact.DR[ind] <- dq.cutoff
    warning("Small eigenvalues are erratic as E'E is numerically singular!")
  }
  t1 <- Sys.time()
  secs <- difftime(t1, t0, units = "secs")[[1L]]
  list(values = exact.DR, secs = secs)
}

Rho2REDF <- function (DR, rho) {
  .Call("C_Rho2REDF", DR, rho, PACKAGE = "gps")
}

REDF2Rho <- function (DR, redf, rho0, MaxNewton, tol = 1e-6) {
  .Call("C_REDF2Rho", DR, redf, rho0, MaxNewton, tol, PACKAGE = "gps")
}

gps2RhoLim <- function (DR, kappa = 0.01) {
  lam.min <- kappa / (1 - kappa) / DR$mean
  rho.min <- log(lam.min)
  lam.max <- (1 - kappa) / kappa / DR$min
  rho.max <- log(lam.max)
  approx.DR <- DR$approx.values
  if (is.na(approx.DR[1L])) {
    rho.max.heuristic <- NA_real_
  } else {
    redf <- length(approx.DR) * kappa
    rho.mid <- 0.5 * (rho.min + rho.max)
    MaxNewton <- 0.25 * (rho.max - rho.min)
    rho.max.heuristic <- REDF2Rho(approx.DR, redf, rho.mid, MaxNewton)
  }
  list(min = rho.min, max = rho.max, max.heuristic = rho.max.heuristic)
}

Q1fun <- function (z, a = 0, b = 1) {
  if (a >= b) stop("a < b required!")
  h <- z * z - z
  theta <- exp(a + (b - a) * z)
  bounds <- c(0, b - a)
  list(h = h, theta = theta, bounds = bounds)
}

Q2fun <- function (z, a = 0, b = 1) {
  if (a >= b) stop("a < b required!")
  y <- 1 - z
  B0 <- y * y * y
  B1 <- 3 * z * y * y
  B2 <- 3 * z * z * y
  B3 <- z * z * z
  h <- B1 - B2
  theta <- exp((B0 + B2) * a + (B2 + B3) * b)
  bounds <- c(a, (2 * a + b) / 3)
  list(h = h, theta = theta, bounds = bounds)
}

map01 <- function (x) (1 / (x[1L] - x[length(x)])) * (x - x[length(x)])

ApproxDR <- function (DR, Qfun, gamma.grid = seq.int(0, 1, 0.05), verbose = TRUE) {
  if (!is.function(Qfun)) stop("'Qfun' is not a function!")
  Qz <- Qfun(z = 0.5, a = 0, b = 1)
  if (anyNA(match(c("h", "theta", "bounds"), names(Qz)))) {
    stop("'Qfun' needs to return elements 'h', 'theta' and 'bounds'!")
  }
  q <- length(DR)
  a <- log(DR[q])
  b <- log(DR[1L])
  dbar <- sum(DR) / q
  logDR <- map01(log(DR))
  step <- 1 / (q + 1)
  p <- seq.int(from = step, by = step, length.out = q)
  aggregated.DR <- numeric(q)
  success <- 0
  op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5))
  on.exit(par(op))
  for (gamma in gamma.grid) {
    z <- map01(log(1 - p) - gamma * log(p))
    zlab <- sprintf("z(%.2f)", gamma)
    Qz <- Qfun(z, a, b)
    plot.default(p, logDR, type = "l", ann = FALSE)
    lines.default(p, z, lty = 2)
    title(sprintf("log(d) and %s", zlab))
    plot.default(z, logDR, type = "l", ann = FALSE)
    title(sprintf("log(d) ~ %s", zlab))
    theta <- Qz$theta
    h <- Qz$h
    eta <- theta * h
    bounds <- Qz$bounds
    tmpDR <- .Call("C_RootApproxDR", theta, h, eta, dbar, bounds, PACKAGE = "gps")
    if (tmpDR[1L]) {
      success <- success + 1
      aggregated.DR <- aggregated.DR + tmpDR
      lines.default(z, map01(log(tmpDR)), lty = 2)
    }
    if (verbose) readline("Hit <Enter> to continue: ")
    if (success && tmpDR[1L] == 0) break
  }
  list(success = success, aggregated.DR = aggregated.DR, gamma.exit = gamma)
}

DemoApproxDR <- function (DR, Qfun = NULL, verbose = TRUE) {
  gamma.grid <- seq.int(0, 1, 0.05)
  if (is.null(Qfun)) {
    out1 <- ApproxDR(DR, Q1fun, gamma.grid, verbose)
    gamma.exit <- out1$gamma.exit
    if (gamma.exit < 1) {
      gamma.grid <- seq.int(gamma.exit, 1, 0.05)
    } else {
      gamma.grid <- numeric(0)
    }
    out2 <- ApproxDR(DR, Q2fun, gamma.grid, verbose)
    success <- out1$success + out2$success
    approx.DR <- (out1$aggregated.DR + out2$aggregated.DR) / success
  } else {
    out <- ApproxDR(DR, Qfun, gamma.grid, verbose)
    success <- out$success
    approx.DR <- out$aggregated.DR / success
  }
  if (success) {
    op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5))
    on.exit(par(op))
    q <- length(DR)
    step <- 1 / (q + 1)
    p <- seq.int(from = step, by = step, length.out = q)
    logDR <- map01(log(DR))
    result <- list(min = DR[q], mean = sum(DR) / q, approx.values = approx.DR)
    rho.lim <- gps2RhoLim(result)
    rho.grid <- seq.int(rho.lim$min, rho.lim$max, length.out = 100)
    redf <- Rho2REDF(DR, rho.grid)
    approx.redf <- Rho2REDF(approx.DR, rho.grid)
    redf.heuristic <- Rho2REDF(DR, rho.lim$max.heuristic)
    mapping <- sprintf("loose: %.2f -> %.2f\nheuristic: %.2f -> %.2f",
                       rho.lim$max, redf[100L], rho.lim$max.heuristic, redf.heuristic)
    plot.default(p, logDR, type = "l", ann = FALSE)
    lines.default(p, map01(log(approx.DR)), lty = 2)
    title("log(d)")
    plot.default(rho.grid, redf, type = "l", ann = FALSE)
    lines.default(rho.grid, approx.redf, lty = 2)
    points.default(rho.lim$max.heuristic, redf.heuristic, pch = 19)
    text.default(rho.lim$max, redf[1L], labels = mapping, adj = c(1, 1))
    title("redf ~ rho")
  } else {
    warning("Unable to approximate Demmler-Reinsch eigenvalues!")
  }
}

GridPWLS <- function (bspl, red, rho) {
  .Call("C_GridPWLS", red$Z, red$L, red$S, red$z, t.default(bspl$C), rho,
        PACKAGE = "gps")
}

GridGCV <- function (red, pwls) {
  .Call("C_GridGCV", red$n, red$L, red$f, red$minRSS, pwls$beta, pwls$edf,
        PACKAGE = "gps")
}

gps2GS <- function (x, y, w = NULL, xt, d = 4, m = 2, gps = TRUE,
                    periodic = FALSE, ng = 20, scalePen = TRUE) {
  bspl <- gps2bspl(x, xt, d, m, gps, periodic)
  red <- gps2red(bspl, y, w, scalePen)
  DR <- gps2DR(red)
  rho.lim <- gps2RhoLim(DR)
  if (ng > 0) {
    if (is.na(rho.lim$max.heuristic)) {
      rho <- seq.int(rho.lim$min, rho.lim$max, length.out = ng)
    } else {
      rho <- seq.int(rho.lim$min, rho.lim$max.heuristic, length.out = ng)
    }
    t0 <- Sys.time()
    pwls <- GridPWLS(bspl, red, rho)
    pwls$gcv <- GridGCV(red, pwls)
    t1 <- Sys.time()
    pwls$secs <- difftime(t1, t0, units = "secs")[[1L]]
  } else {
    pwls <- NULL
  }
  S <- LTB2Csparse(bspl$S, symmetric = TRUE)
  BWB <- LTB2Csparse(red$Z, symmetric = TRUE)
  eqn <- list(B = bspl$B, BWB = BWB, BWy = red$z,
              omega = red$scalePen.fctr, S = S, C = bspl$C)
  list(eqn = eqn, eigen = DR, rho.lim = rho.lim, E = red$E, pwls = pwls)
}

DemoRhoLim <- function (fit, plot = TRUE) {
  DR <- ExactDR(fit$E)
  approx.DR <- fit$eigen[c("approx.values", "secs")]
  names(approx.DR) <- c("values", "secs")
  q <- length(DR$values)
  rho.min.star <- fit$rho.lim$min
  rho.max.star <- fit$rho.lim$max
  redf.min <- 0.01 * q
  redf.max <- 0.99 * q
  rho.mid <- (rho.min.star + rho.max.star) / 2
  MaxNewton <- 0.25 * (rho.max.star - rho.min.star)
  rho.min <- REDF2Rho(DR$values, redf.max, rho.mid, MaxNewton)
  rho.max <- REDF2Rho(DR$values, redf.min, rho.mid, MaxNewton)
  ng <- 100L
  rho.grid <- seq.int(rho.min.star, rho.max.star, length.out = ng)
  redf.exact <- Rho2REDF(DR$values, rho.grid)
  redf.min.star <- redf.exact[ng]
  redf.max.star <- redf.exact[1L]
  rho.max.heuristic <- fit$rho.lim$max.heuristic
  if (is.na(rho.max.heuristic)) {
    redf.approx <- rep.int(NA_real_, ng)
    redf.min.heuristic <- NA_real_
    warning("No approximated eigenvalues and heuristic upper bound for 'rho'!")
  } else {
    redf.approx <- Rho2REDF(approx.DR$values, rho.grid)
    redf.min.heuristic <- Rho2REDF(DR$values, rho.max.heuristic)
  }
  map.exact <- sprintf("[%.1f, %.1f] -> [%.1f, %.1f]",
                       rho.min, rho.max, redf.min, redf.max)
  map.wider <- sprintf("[%.1f, %.1f] -> [%.1f, %.1f]",
                       rho.min.star, rho.max.star, redf.min.star, redf.max.star)
  map.heuristic <- sprintf("%.1f -> %.1f", rho.max.heuristic, redf.min.heuristic)
  if (plot) {
    op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5))
    on.exit(par(op))
    p <- (1:q) / (q + 1)
    plot.default(p, log(DR$values), type = "l", ann = FALSE)
    lines.default(p, log(approx.DR$values), lty = 2)
    text.default(p[q], log(DR$values[1L]), adj = c(1, 1),
                 labels = sprintf("exact: solid\napprox: dashed"))
    title("log eigenvalues")
    plot.default(rho.grid, redf.exact, type = "l", ann = FALSE)
    lines.default(rho.grid, redf.approx, lty = 2)
    segments(rho.min, q, rho.min, 0.75 * q, col = 8, lty = 2)
    segments(rho.max, 0, rho.max, 0.25 * q, col = 8, lty = 2)
    points.default(rho.max.heuristic, redf.min.heuristic, pch = 19)
    text.default(rho.max.star, redf.max.star, adj = c(1, 1),
                 labels = sprintf("exact: %s\nwider: %s\nheuristic: %s",
                                  map.exact, map.wider, map.heuristic))
    title("redf ~ rho")
  }
  eigen <- list(approx = approx.DR, exact = DR)
  limit <- list(rho = list(exact = c(rho.min, rho.max),
                           wider = c(rho.min.star, rho.max.star),
                           heuristic = rho.max.heuristic),
                redf = list(exact = c(redf.min, redf.max),
                            wider = c(redf.min.star, redf.max.star),
                            heuristic = redf.min.heuristic))
  mapping <- list(exact = map.exact, wider = map.wider, heuristic = map.heuristic)
  redf <- list(rho = rho.grid, exact = redf.exact, approx = redf.approx)
  list(eigen = eigen, limit = limit, mapping = mapping, redf = redf)
}


Diff <- function (x, k = 1L, n = length(x), xi = 1L) {
  .Call("C_Diff", x, k, n, xi, PACKAGE = "gps")
}

SparseDelta <- function (r) {
  r <- as.integer(r)
  x <- rep.int(c(-1, 1), r)
  i <- rep(seq.int(0L, r - 1L), each = 2L)
  p <- c(0L, seq.int(1L, length.out = r, by = 2L), 2L * r)
  methods::new("dgCMatrix", i = i, p = p, Dim = c(r, r + 1L), x = x)
}

SparseWtDelta <- function (h) {
  r <- length(h)
  x <- rep.int(c(-1, 1), r) * rep(1 / h, each = 2)
  i <- rep(seq.int(0L, r - 1L), each = 2L)
  p <- c(0L, seq.int(1L, length.out = r, by = 2L), 2L * r)
  methods::new("dgCMatrix", i = i, p = p, Dim = c(r, r + 1L), x = x)
}

SparseGD <- function (xt, d, m) {
  K <- length(xt)
  D <- vector("list", m)
  h <- Diff(x = xt, k = d - 1L, n = K - 2L, xi = 2L)
  D[[1L]] <- SparseWtDelta(h)
  i <- 2L
  while (i <= m) {
    h <- Diff(x = xt, k = d - i, n = K - 2L * i, xi = i + 1L)
    D[[i]] <- SparseWtDelta(h) %*% D[[i - 1L]]
    i <- i + 1L
  }
  D
}

DiffCoef <- function (b, xt, d, m) {
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  if (m < 1 || m >= d) stop("1 <= m <= d - 1 required!", call. = FALSE)
  xt <- MonoKnots(xt, d)
  K <- length(xt)
  if (length(b) != K - d) {
    stop("length(b) == length(xt) - d required!", call. = FALSE)
  }
  b <- as.double(b)
  for (i in 1:m) {
    h <- Diff(x = xt, k = d - i, n = K - 2L * i, xi = i + 1L)
    b <- Diff(b) / h
  }
  b
}

ComputeLD <- function (xt, d) {
  .Call("C_ComputeLD", xt, d, PACKAGE = "gps")
}

NullGD <- function (ld, m = 1, orthonormal = TRUE) {
  basis <- .Call("C_NullGD", ld, m, PACKAGE = "gps")
  if (orthonormal && (m > 1)) {
    Q <- qr.Q(qr.default(basis[, m:1]))[, m:1]
    i <- sequence.default(1:(m - 1))
    j <- rep.int(2:m, 1:(m - 1))
    Q[cbind(i, j)] <- 0
    basis <- Q
  }
  basis
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
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  if (m < 0 || m >= d) stop("0 <= m <= d - 1 required!", call. = FALSE)
  xt <- MonoKnots(xt, d)
  K <- length(xt)
  ord <- d - m
  if (ord == 1) {
    h <- Diff(x = xt, n = K - 2 * (d - 1), xi = d)
    return(h)
  }
  W <- QuadWts(ord)
  xd <- xt[seq.int(d, K - d + 1)]
  xg <- MakeGrid(xd, ord)
  xt.local <- xt[seq.int(1 + m, K - m)]
  B <- splines::splineDesign(xt.local, xg, ord, sparse = TRUE)
  .Call("C_SbarBlocks", xd, W, B@x, PACKAGE = "gps")
}

GramBS <- function (xt, d) {
  blocks <- SbarBlocks(xt, d, m = 0)
  LTB <- .Call("C_SbarLTB", blocks, PACKAGE = "gps")
  LTB2Csparse(LTB, symmetric = TRUE)
}

btSb <- function (b, xt, d, m) {
  db <- DiffCoef(b, xt, d, m)
  blocks <- SbarBlocks(xt, d, m)
  .Call("C_btSb", blocks, db, PACKAGE = "gps")
}


SparseD <- function (xt, d, m = NULL, gps = TRUE) {
  d <- as.integer(d)
  if (d < 2L) stop("d >= 2 required!", call. = FALSE)
  if (is.null(m)) {
    m <- seq_len(d - 1L)
  } else {
    m <- sort.int(as.integer(m), method = "radix")
  }
  m.min <- m[1L]
  m.max <- m[length(m)]
  if (m.min < 1L || m.max >= d) {
    stop("1 <= m <= d - 1 required!", call. = FALSE)
  }
  NAMES <- sprintf("order.%d", m)
  xt <- MonoKnots(xt, d)
  D <- SparseGD(xt, d, m.max)
  D <- D[m]
  names(D) <- NAMES
  if (gps) return(D)
  nm <- length(m)
  Sbar <- vector("list", nm); names(Sbar) <- NAMES
  K <- vector("list", nm); names(K) <- NAMES
  for (i in 1:nm) {
    m_i <- m[i]
    D_i <- D[[i]]
    blocks <- SbarBlocks(xt, d, m_i)
    if (is.array(blocks)) {
      Sbar.LTB <- .Call("C_SbarLTB", blocks, PACKAGE = "gps")
      Sbar.Matrix <- LTB2Csparse(Sbar.LTB, symmetric = TRUE)
      L.LTB <- LPBTRF(Sbar.LTB, overwrite = TRUE)
      L.Matrix <- LTB2Csparse(L.LTB)
      K_i <- Matrix::crossprod(L.Matrix, D_i)
    } else {
      k1 <- length(blocks)
      Sbar.Matrix <- methods::new("ddiMatrix", diag = "N", x = blocks, Dim = c(k1, k1))
      L.diag <- sqrt(blocks)
      newx <- L.diag[D_i@i + 1L] * D_i@x
      K_i <- methods::new("dgCMatrix", i = D_i@i, p = D_i@p, Dim = D_i@Dim, x = newx)
    }
    Sbar[[i]] <- Sbar.Matrix
    K[[i]] <- K_i
  }
  sandwich <- list(D = D, Sbar = Sbar)
  structure(K, sandwich = sandwich)
}



PriorCoef1 <- function (n, D) {
  q <- D@Dim[1L]
  S <- as_matrix(Matrix::crossprod(D))
  ei <- base::eigen(S, symmetric = TRUE)
  qs <- seq_len(q)
  d <- ei$values[qs]
  U <- ei$vectors[, qs]
  di <- 1 / d
  e <- rnorm(n * q)
  if (n > 1) e <- ChangeDim(e, c(q, n))
  b <- U %*% (sqrt(di) * e)
  if (n == 1) b <- c(b)
  b
}

PriorCoef2 <- function (n, D) {
  Dim <- D@Dim
  q <- Dim[1L]
  p <- Dim[2L]
  m <- p - q
  e <- rnorm(n * q)
  if (n > 1) e <- ChangeDim(e, c(q, n))
  Dt <- Csparse2LTB(D)
  y <- SolveLTB(transA = TRUE, Dt, e, overwrite = TRUE)
  R <- as_matrix(D[, seq.int(q + 1, p), drop = FALSE])
  X <- SolveLTB(transA = TRUE, Dt, R, overwrite = TRUE)
  XtX <- base::crossprod(X)
  Xty <- base::crossprod(X, y)
  U <- chol.default(XtX + diag(m))
  f <- forwardsolve(U, Xty, upper.tri = TRUE, transpose = TRUE)
  b <- backsolve(U, f)
  r <- y - X %*% b
  if (n > 1) rbind(r, b) else c(r, b)
}

PriorCoef <- function (n, D) {
  PriorCoef2(n, D)
}

Dt2ThinQR <- function (Dt) {
  A <- DDt(Dt)
  Rt <- LPBTRF(A, overwrite = TRUE)
  denseDt <- LTB2Dense(Dt, sum(dim(Dt)) - 1L)  ## does not work for O-spline
  Qt <- SolveLTB(transA = FALSE, Rt, t.default(denseDt), overwrite = TRUE)
  list(Qt = Qt, Rt = Rt)
}

MPinvUUt <- function (D, method = "qr") {
  if (method == "qr") {
    Dt <- Csparse2LTB(D)
    QR <- Dt2ThinQR(Dt)
    Qt <- QR$Qt
    Rt <- QR$Rt
    SolveLTB(transA = TRUE, Rt, Qt, overwrite = TRUE)
  } else if (method == "eigen") {
    q <- D@Dim[1L]
    S <- as_matrix(Matrix::crossprod(D))
    ei <- base::eigen(S, symmetric = TRUE)
    qs <- seq_len(q)
    d <- ei$values[qs]
    Q <- ei$vectors[, qs]
    di <- 1 / d
    sqrt(di) * t.default(Q)
  } else {
    stop("'method' must be either \"qr\" or \"eigen\"!")
  }
}

MPinv <- function (D, only.diag = FALSE) {
  Ut <- MPinvUUt(D, method = "qr")
  if (only.diag) {
    base::colSums(Ut * Ut)
  } else {
    base::crossprod(Ut)
  }
}



GetPBC <- function (xt, d, compact = TRUE, transpose = FALSE, sparse = TRUE) {
  d <- as.integer(d)
  if (d < 2L) stop("d >= 2 required!")
  degree <- d - 1L
  K <- length(xt)
  if (length(xt) < degree + 2L * d) {
    stop("length(xt) >= 3 * d - 1 required!")
  }
  p <- K - d
  a <- xt[d]
  b <- xt[K - degree]
  nDeriv <- seq.int(0L, length.out = degree)
  knots.a <- xt[seq.int(from = 1L, length.out = 2L * d)]
  knots.b <- xt[seq.int(to = K, length.out = 2L * d)]
  Ca <- splines::splineDesign(knots.a, rep.int(a, degree), d, nDeriv)
  Cb <- splines::splineDesign(knots.b, rep.int(b, degree), d, nDeriv)
  Ca <- Ca[, seq.int(1L, degree)]
  Cb <- Cb[, seq.int(2L, d)]
  C.compact <- cbind(Ca, -Cb, deparse.level = 0L)
  Ct.compact <- t.default(C.compact)
  if (compact) {
    if (transpose) Ct.compact else C.compact
  } else {
    if (sparse) {
      if (transpose) {
        Ct.i <- rep.int(c(nDeriv, seq.int(to = p - 1L, length.out = degree)), degree)
        Ct.p <- seq.int(0L, by = 2L * degree, length.out = d)
        ChangeDim(Ct.compact, NULL)
        Ct.x <- Ct.compact
        methods::new("dgCMatrix", i = Ct.i, p = Ct.p, Dim = c(p, degree), x = Ct.x)
      } else {
        blocksize <- degree * degree
        C.i <- rep.int(nDeriv, 2L * degree)
        C.p <- c(seq.int(0, by = degree, length.out = d),
                 rep.int(blocksize, p - 2L * degree),
                 seq.int(blocksize + degree, by = degree, length.out = degree))
        ChangeDim(C.compact, NULL)
        C.x <- C.compact
        methods::new("dgCMatrix", i = C.i, p = C.p, Dim = c(degree, p), x = C.x)
      }
    } else {
      ind <- c(seq.int(1L, degree), seq.int(to = p, length.out = degree))
      if (transpose) {
        Ct <- matrix(0, p, degree)
        Ct[ind, ] <- Ct.compact
        Ct
      } else {
        C <- matrix(0, degree, p)
        C[, ind] <- C.compact
        C
      }
    }
  }
}

NullPBC <- function (xt, d, compact = TRUE) {
  d <- as.integer(d)
  degree <- d - 1L
  Ct <- GetPBC(xt, d, compact = TRUE, transpose = TRUE)
  QR <- qr.default(Ct)
  Q <- qr.Q(QR, complete = TRUE)
  Q.compact <- Q[, seq.int(d, 2L * degree), drop = FALSE]
  if (compact) {
    Q.compact
  } else {
    p <- length(xt) - d
    pp <- p - degree
    N <- pp - degree
    Q.compact.i <- c(seq.int(0L, length.out = degree),
                     seq.int(to = p - 1L, length.out = degree))
    Q.i <- c(seq.int(degree, length.out = N), rep.int(Q.compact.i, degree))
    Q.p <- c(seq.int(0, length.out = N),
             seq.int(N, by = 2L * degree, length.out = d))
    Q.x <- c(rep.int(1, N), Q.compact)
    methods::new("dgCMatrix", i = Q.i, p = Q.p, Dim = c(p, pp), x = Q.x)
  }
}

pbsDesign <- function (x, xd, d, nDeriv = 0, sparse = FALSE, wrap = TRUE) {
  d <- as.integer(d)
  if (d < 2L) stop("d >= 2 required!")
  degree <- d - 1L
  if (nDeriv < 0 || nDeriv > degree) stop("0 <= nDeriv <= d - 1 required!")
  xd <- as.double(xd)
  if (!IsMonoInc(xd)) stop("'xd' is not strictly ascending!")
  k2 <- length(xd)
  if (k2 <= d) stop("length(xd) >= d + 1 required!")
  a <- xd[1L]
  b <- xd[k2]
  if (min(x) < a || max(x) > b) stop("domain does not contain all x-values!")
  if (wrap) {
    pbsDesign1(x, xd, d, nDeriv, sparse)
  } else {
    pbsDesign2(x, xd, d, nDeriv, sparse)
  }
}

pbsDesign1 <- function (x, xd, d, nDeriv = 0, sparse = FALSE) {
  degree <- d - 1L
  k2 <- length(xd)
  p <- k2 - 1L
  a <- xd[1L]
  b <- xd[k2]
  period <- b - a
  raux <- xd[2:d] + period
  xt <- c(xd, raux)
  b.aux <- raux[degree]
  ind.x.wrapped <- which(x < xd[d])
  if (length(ind.x.wrapped) == 0L) {
    PB <- splines::splineDesign(xt, x, d, nDeriv, sparse = sparse)
  } else {
    nx <- length(x)
    x.wrapped <- x[ind.x.wrapped] + period
    x.wrapped[x.wrapped < b] <- b
    x.wrapped[x.wrapped > b.aux] <- b.aux
    x.padded <- c(x, x.wrapped)
    xt0 <- c(rep.int(a, degree), xt, rep.int(b.aux, degree))
    B <- splines::splineDesign(xt0, x.padded, d, nDeriv, sparse = TRUE)
    B.p <- B@p[seq.int(d, length(B@p) - degree)]
    ind <- seq.int(B.p[1L] + 1L, B.p[p - degree + 1L])
    unchanged.i <- B@i[ind]
    unchanged.x <- B@x[ind]
    B.p.vital <- B.p[seq.int(to = p + 1L, length.out = d)]
    ind <- seq.int(B.p.vital[1L] + 1L, B.p.vital[d])
    unordered.i <- B@i[ind]
    unordered.x <- B@x[ind]
    ind <- (unordered.i >= nx)
    unordered.i[ind] <- ind.x.wrapped[unordered.i[ind] + 1L - nx] - 1L
    unordered.j <- rep.int(1:degree, diff.default(B.p.vital))
    ind <- order(unordered.j, unordered.i, method = "radix")
    reordered.i <- unordered.i[ind]
    reordered.x <- unordered.x[ind]
    PB.p <- B.p - B.p[1L]
    PB.i <- c(unchanged.i, reordered.i)
    PB.x <- c(unchanged.x, reordered.x)
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

pbsDesign2 <- function (x, xd, d, nDeriv = 0, sparse = FALSE) {
  degree <- d - 1L
  k2 <- length(xd)
  xt <- c(rep.int(xd[1L], degree), xd, rep.int(xd[k2], degree))
  B <- splines::splineDesign(xt, x, d, nDeriv, sparse = TRUE)
  p <- ncol(B)
  ind <- seq.int(B@p[d] + 1L, B@p[p + 1L - degree])
  PB.i <- B@i[ind]
  PB.x <- B@x[ind]
  PB.p <- B@p[seq.int(d, p + 1L)]
  PB.p <- PB.p - PB.p[1L]
  IND <- seq.int(to = k2, length.out = degree)
  PB.p[IND] <- PB.p[k2 - degree]
  bnd_cols <- c(seq.int(1L, degree), seq.int(k2, p))
  nnz_per_col <- B@p[bnd_cols + 1L] - B@p[bnd_cols]
  active_col <- nnz_per_col > 0L
  bnd_cols <- bnd_cols[active_col]
  num_active_cols <- length(bnd_cols)
  if (num_active_cols > 0L) {
    Q <- NullPBC(xt, d)
    Q <- Q[active_col, , drop = FALSE]
    i_lst <- x_lst <- vector("list", num_active_cols)
    for (i in 1:num_active_cols) {
      j <- bnd_cols[i]
      ind <- seq.int(B@p[j] + 1L, B@p[j + 1L])
      i_lst[[i]] <- B@i[ind]
      x_lst[[i]] <- B@x[ind]
    }
    active_row <- sort.int(unique.default(unlist(i_lst)), method = "radix")
    num_active_rows <- length(active_row)
    block.i <- unlist(lapply(i_lst, match, active_row))
    block.j <- rep.int(1:num_active_cols, nnz_per_col[active_col])
    block.x <- unlist(x_lst)
    B.boundary <- matrix(0, num_active_rows, num_active_cols)
    B.boundary[cbind(block.i, block.j)] <- block.x
    B.transformed <- B.boundary %*% Q
    PB.i <- c(PB.i, rep.int(active_row, degree))
    PB.x <- c(PB.x, B.transformed)
    PB.p[IND] <- PB.p[IND] + num_active_rows * (1:degree)
  }
  nr <- length(x)
  nc <- k2 - 1L
  if (sparse) {
    PB <- methods::new("dgCMatrix", i = PB.i, p = PB.p, Dim = c(nr, nc), x = PB.x)
  } else {
    PB.j <- rep.int(1:nc, diff.default(PB.p))
    PB <- matrix(0, nr, nc)
    PB[cbind(PB.i + 1L, PB.j)] <- PB.x
  }
  PB
}

SparsePD <- function (xd, d, wrap = TRUE) {
  d <- as.integer(d)
  if (d < 2L) stop("d >= 2 required!")
  degree <- d - 1L
  xd <- as.double(xd)
  if (!IsMonoInc(xd)) stop("'xd' is not strictly ascending!")
  k2 <- length(xd)
  if (k2 <= d) stop("length(xd) >= d + 1 required!")
  if (wrap) {
    SparsePD1(xd, d)
  } else {
    SparsePD2(xd, d)
  }
}

SparsePD1 <- function (xd, d) {
  degree <- d - 1L
  k2 <- length(xd)
  p <- k2 - 1L
  a <- xd[1L]
  b <- xd[k2]
  period <- b - a
  laux <- xd[(k2 - degree):p] - period
  raux <- xd[2:d] + period
  xt <- c(laux, xd, raux)
  D <- SparseD(xt, d)
  PD <- vector("list", degree)
  
  for (m in 1:degree) {
    PD[[m]] <- D[[m]][, 1:p]
    PD[[m]][, 1:degree] <- PD[[m]][, 1:degree] + D[[m]][, p + 1:degree]
  }
  PD
}

SparsePD2 <- function (xd, d) {
  degree <- d - 1L
  k2 <- length(xd)
  xt <- c(rep.int(xd[1L], degree), xd, rep.int(xd[k2], degree))
  D <- SparseD(xt, d)
  Q <- NullPBC(xt, d, compact = FALSE)
  PD <- vector("list", degree)
  for (m in 1:degree) PD[[m]] <- D[[m]] %*% Q
  m <- 1L
  while (m < degree) {
    Qm <- NullPBC(xt[seq.int(1L + m, length(xt) - m)], d - m, compact = FALSE)
    PD[[m]] <- Matrix::crossprod(Qm, PD[[m]])
    m <- m + 1L
  }
  PD
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

VecDot <- function (x, y = NULL) {
  if (is.null(y)) y <- x
  .Call("C_VecDot", x, y, PACKAGE = "gps")
}

VecScal <- function (a, x, overwrite = FALSE) {
  .Call("C_VecScal", a, x, overwrite, PACKAGE = "gps")
}


ChangeDim <- function (x, Value) {
  .Call("C_SetDim", x, Value, PACKAGE = "gps")
}
