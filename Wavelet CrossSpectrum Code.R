require(WaveletComp)

################
# load data
################
cholera_df <- read.csv("cholera.csv")

#####################################
# set some basic plotting parameters
#####################################
color_key = "interval"    # equispaced breakpoints in color scale (= "quantile" for default)
which_image = "wp"        # plot cross-power (= "wc" to plot coherence)
sel_lower = 8             # cut-off scale/period for reconstruction (= 2 for default)
dj = 1/100                # frequency resolution (1 / no. of voices per octave)
n_sim = 10                # no. of simulations
ticks  <- seq(2*12+6, nrow(cholera_df), by = 12*10)   # x-tick marks for plotting
labels <- seq(1895, 1940, by = 10)                    # x-tick labels for plotting

##########################################################
# compute cross-wavelet power of original pair-wise series
##########################################################
# (Dacca, Rangpur) =>
w1 <- analyze.coherency(my.data = cholera_df, my.pair = c("Dacca", "Rangpur"), loess.span = 0,
  dt = 1, dj = dj, lowerPeriod = 2, upperPeriod = 128,
  make.pval = TRUE, n.sim = n_sim, verbose = FALSE)

# (Dacca, Mymensing) =>
w2 <- analyze.coherency(my.data = cholera_df, my.pair = c("Dacca", "Mymensing"), loess.span = 0,
  dt = 1, dj = dj, lowerPeriod = 2, upperPeriod = 128,
  make.pval = TRUE, n.sim = n_sim, verbose = FALSE)

# (Dacca, Parganas) =>
w3 <- analyze.coherency(my.data = cholera_df, my.pair = c("Dacca", "Parganas"), loess.span = 0,
  dt = 1, dj = dj, lowerPeriod = 2, upperPeriod = 128,
  make.pval = TRUE, n.sim = n_sim, verbose = FALSE)

# (Dacca, Calcutta) =>
w4 <- analyze.coherency(my.data = cholera_df, my.pair = c("Dacca", "Calcutta"), loess.span = 0,
  dt = 1, dj = dj, lowerPeriod = 2, upperPeriod = 128,
  make.pval = TRUE, n.sim = n_sim, verbose = FALSE)

#######################################################
# plot cross-wavelet power of original pair-wise series
#######################################################
# (Dacca, Rangpur) =>
wc.image(w1, which.image = which_image, color.key = color_key,
  main = "(Dacca, Rangpur) (original)",
  n.levels = 250, legend.params = list(lab = "power levels"),
  periodlab = "period (month)", timelab = "",
  plot.contour = TRUE, col.contour = "black", plot.ridge = FALSE, plot.arrow = FALSE,
  spec.time.axis = list(at = ticks, labels = labels))

# (Dacca, Mymensing) =>
dev.new()
wc.image(w2, which.image = which_image, color.key = color_key,
  main = "(Dacca, Mymensing) (original)",
  n.levels = 250, legend.params = list(lab = "power levels"),
  periodlab = "period (month)", timelab = "",
  plot.contour = TRUE, col.contour = "black", plot.ridge = FALSE, plot.arrow = FALSE,
  spec.time.axis = list(at = ticks, labels = labels))

# (Dacca, Parganas) =>
dev.new()
wc.image(w3, which.image = which_image, color.key = color_key,
  main = "(Dacca, Parganas) (original)",
  n.levels = 250, legend.params = list(lab = "power levels"),
  periodlab = "period (month)", timelab = "",
  plot.contour = TRUE, col.contour = "black", plot.ridge = FALSE, plot.arrow = FALSE,
  spec.time.axis = list(at = ticks, labels = labels))

# (Dacca, Calcutta) =>
dev.new()
wc.image(w4, which.image = which_image, color.key = color_key,
  main = "(Dacca, Calcutta) (original)",
  n.levels = 250, legend.params = list(lab = "power levels"),
  periodlab = "period (month)", timelab = "",
  plot.contour = TRUE, col.contour = "black", plot.ridge = FALSE, plot.arrow = FALSE,
  spec.time.axis = list(at = ticks, labels = labels))

#######################################################
# reconstruct each series with scales < 8-month removed
#######################################################
# first step: compute standard wavelet transform of each series
iw1 <- analyze.wavelet(cholera_df, "Dacca", loess.span = 0, dt = 1, dj = 1/250,
  lowerPeriod = 2, upperPeriod = 128, make.pval = TRUE, n.sim = 20, verbose = FALSE)
iw2 <- analyze.wavelet(cholera_df, "Rangpur", loess.span = 0, dt = 1, dj = 1/250,
  lowerPeriod = 2, upperPeriod = 128, make.pval = TRUE, n.sim = 20, verbose = FALSE)
iw3 <- analyze.wavelet(cholera_df, "Mymensing", loess.span = 0, dt = 1, dj = 1/250,
  lowerPeriod = 2, upperPeriod = 128, make.pval = TRUE, n.sim = 20, verbose = FALSE)
iw4 <- analyze.wavelet(cholera_df, "Parganas", loess.span = 0, dt = 1, dj = 1/250,
  lowerPeriod = 2, upperPeriod = 128, make.pval = TRUE, n.sim = 20, verbose = FALSE)
iw5 <- analyze.wavelet(cholera_df, "Calcutta", loess.span = 0, dt = 1, dj = 1/250,
  lowerPeriod = 2, upperPeriod = 128, make.pval = TRUE, n.sim = 20, verbose = FALSE)
# next step: reconstruct each series
r1 <- reconstruct(iw1, only.sig = FALSE, sel.lower = sel_lower,
  plot.waves = FALSE, plot.rec = FALSE, verbose = FALSE)
r2 <- reconstruct(iw2, only.sig = FALSE, sel.lower = sel_lower,
  plot.waves = FALSE, plot.rec = FALSE, verbose = FALSE)
r3 <- reconstruct(iw3, only.sig = FALSE, sel.lower = sel_lower,
  plot.waves = FALSE, plot.rec = FALSE, verbose = FALSE)
r4 <- reconstruct(iw4, only.sig = FALSE, sel.lower = sel_lower,
  plot.waves = FALSE, plot.rec = FALSE, verbose = FALSE)
r5 <- reconstruct(iw5, only.sig = FALSE, sel.lower = sel_lower,
  plot.waves = FALSE, plot.rec = FALSE, verbose = FALSE)

####################################################################
# now compute cross-wavelet power of reconstructed pair-wise series
####################################################################
# (Dacca, Rangpur) =>
df1 = data.frame(x = r1$series$Dacca.r, y = r2$series$Rangpur.r)
w5 <- analyze.coherency(my.data = df1, my.pair = c("x", "y"), loess.span = 0,
  dt = 1, dj = dj, lowerPeriod = 2, upperPeriod = 128,
  make.pval = TRUE, n.sim = n_sim, verbose = FALSE)

# (Dacca, Mymensing) =>
df2 = data.frame(x = r1$series$Dacca.r, y = r3$series$Mymensing.r)
w6 <- analyze.coherency(my.data = df2, my.pair = c("x", "y"), loess.span = 0,
  dt = 1, dj = dj, lowerPeriod = 2, upperPeriod = 128,
  make.pval = TRUE, n.sim = n_sim, verbose = FALSE)

# (Dacca, Parganas) =>
df3 = data.frame(x = r1$series$Dacca.r, y = r4$series$Parganas.r)
w7 <- analyze.coherency(my.data = df3, my.pair = c("x", "y"), loess.span = 0,
  dt = 1, dj = dj, lowerPeriod = 2, upperPeriod = 128,
  make.pval = TRUE, n.sim = n_sim, verbose = FALSE)

# (Dacca, Calcutta) =>
df4 = data.frame(x = r1$series$Dacca.r, y = r5$series$Calcutta.r)
w8 <- analyze.coherency(my.data = df4, my.pair = c("x", "y"), loess.span = 0,
  dt = 1, dj = dj, lowerPeriod = 2, upperPeriod = 128,
  make.pval = TRUE, n.sim = n_sim, verbose = FALSE)

####################################################################
# finally plot cross-wavelet power of reconstructed pair-wise series
####################################################################
# (Dacca, Rangpur) =>
dev.new()
wc.image(w5, which.image = which_image, color.key = color_key,
  main = "(Dacca, Rangpur) (reconstructed)",
  n.levels = 250, legend.params = list(lab = "power levels"),
  periodlab = "period (month)", timelab = "",
  plot.contour = TRUE, col.contour = "black", plot.ridge = FALSE, plot.arrow = FALSE,
  spec.time.axis = list(at = ticks, labels = labels))

# (Dacca, Mymensing) =>
dev.new()
wc.image(w6, which.image = which_image, color.key = color_key,
  main = "(Dacca, Mymensing) (reconstructed)",
  n.levels = 250, legend.params = list(lab = "power levels"),
  periodlab = "period (month)", timelab = "",
  plot.contour = TRUE, col.contour = "black", plot.ridge = FALSE, plot.arrow = FALSE,
  spec.time.axis = list(at = ticks, labels = labels))

# (Dacca, Parganas) =>
dev.new()
wc.image(w7, which.image = which_image, color.key = color_key,
  main = "(Dacca, Parganas) (reconstructed)",
  n.levels = 250, legend.params = list(lab = "power levels"),
  periodlab = "period (month)", timelab = "",
  plot.contour = TRUE, col.contour = "black", plot.ridge = FALSE, plot.arrow = FALSE,
  spec.time.axis = list(at = ticks, labels = labels))

# (Dacca, Calcutta) =>
dev.new()
wc.image(w8, which.image = which_image, color.key = color_key,
  main = "(Dacca, Calcutta) (reconstructed)",
  n.levels = 250, legend.params = list(lab = "power levels"),
  periodlab = "period (month)", timelab = "",
  plot.contour = TRUE, col.contour = "black", plot.ridge = FALSE, plot.arrow = FALSE,
  spec.time.axis = list(at = ticks, labels = labels))
