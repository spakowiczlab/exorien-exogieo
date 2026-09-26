# Hex sticker and graphical abstract for the exogieo manuscript repository.
# Run from the repository root:
#   Rscript scripts/figures_hex-and-abstract.R

library(ggplot2)
library(hexSticker)
library(grid)
library(showtext)

font_add(
  "Arial",
  regular = "/System/Library/Fonts/Supplemental/Arial.ttf",
  bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
)
showtext_auto()

ink <- "#14343C"
teal <- "#0E6B6B"
gold <- "#E0A106"
coral <- "#D4654A"
cream <- "#F7F4EF"
paper <- "#FFFcf7"

# --- hex sticker ----------------------------------------------------------

clock <- data.frame(
  x = cos(seq(0, 2 * pi, length.out = 120)),
  y = sin(seq(0, 2 * pi, length.out = 120))
)
tick_ang <- seq(0, 2 * pi, length.out = 13)[-13]
ticks <- data.frame(
  x = 0.78 * cos(tick_ang),
  y = 0.78 * sin(tick_ang),
  xend = 0.98 * cos(tick_ang),
  yend = 0.98 * sin(tick_ang)
)
# Chronologic hand near 10 o'clock; methylation-age hand pushed ahead.
hands <- data.frame(
  x = 0, y = 0,
  xend = c(0.42 * cos(2.3), 0.72 * cos(0.35)),
  yend = c(0.42 * sin(2.3), 0.72 * sin(0.35))
)
capsule <- function(x, y, len, wid, ang, n = 30) {
  theta <- seq(-pi / 2, pi / 2, length.out = n)
  right <- cbind(len / 2 + (wid / 2) * cos(theta), (wid / 2) * sin(theta))
  theta2 <- seq(pi / 2, 3 * pi / 2, length.out = n)
  left <- cbind(-len / 2 + (wid / 2) * cos(theta2), (wid / 2) * sin(theta2))
  pts <- rbind(right, left)
  rot <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), 2)
  pts <- t(rot %*% t(pts))
  data.frame(x = pts[, 1] + x, y = pts[, 2] + y)
}

bugs <- rbind(
  cbind(capsule(-0.58, -1.22, 0.46, 0.16, -0.55), id = 1),
  cbind(capsule(0.02, -1.50, 0.50, 0.16, 0.05), id = 2),
  cbind(capsule(0.62, -1.22, 0.46, 0.16, 0.55), id = 3)
)

icon <- ggplot() +
  geom_polygon(data = bugs, aes(x, y, group = id), fill = coral, color = NA) +
  geom_path(data = clock, aes(x, y), color = gold, linewidth = 1.15) +
  geom_segment(
    data = ticks,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "white",
    linewidth = 0.45
  ) +
  geom_segment(
    data = hands,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "white",
    linewidth = 1.15,
    lineend = "round"
  ) +
  annotate("point", x = 0, y = 0, color = gold, size = 2.4) +
  coord_fixed(xlim = c(-1.45, 1.45), ylim = c(-1.75, 1.35), clip = "off") +
  theme_void() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  )

sticker(
  icon,
  package = "exogieo",
  p_size = 22,
  p_y = 1.48,
  p_color = "white",
  p_family = "Aller_Rg",
  p_fontface = "plain",
  s_x = 1,
  s_y = 0.82,
  s_width = 1.05,
  s_height = 0.95,
  h_fill = ink,
  h_color = gold,
  h_size = 1.4,
  url = "EOCRC",
  u_size = 4.6,
  u_color = gold,
  u_family = "Aller_Rg",
  white_around_sticker = FALSE,
  filename = "figures/hex_exogieo.png",
  dpi = 320
)

# --- graphical abstract ---------------------------------------------------

png(
  "figures/graphical_abstract.png",
  width = 13.2,
  height = 7.15,
  units = "in",
  res = 180,
  bg = cream
)
grid.newpage()

round_box <- function(x, y, w, h, fill, border, lwd = 1.6) {
  grid.roundrect(
    x = x, y = y, width = w, height = h,
    r = unit(0.14, "inches"),
    gp = gpar(fill = fill, col = border, lwd = lwd)
  )
}

label <- function(txt, x, y, just = "left", size = 11, face = "plain", col = ink) {
  grid.text(
    txt, x = x, y = y, just = just,
    gp = gpar(fontfamily = "Arial", fontface = face, fontsize = size, col = col)
  )
}

# Header
grid.roundrect(
  x = 0.5, y = 0.915, width = 0.96, height = 0.13,
  r = unit(0.12, "inches"),
  gp = gpar(fill = ink, col = NA)
)
label(
  "EPIGENETIC MODULATION, INTRATUMORAL MICROBIOME, AND IMMUNITY",
  0.04, 0.945, size = 15.5, face = "bold", col = "white"
)
label(
  "Early-onset colorectal cancer  ·  Jin et al., Cancer Research Communications (2025)",
  0.04, 0.885, size = 10.5, col = "#F4C95D"
)

# Cohort strip
round_box(0.5, 0.785, 0.96, 0.085, paper, "#E4DDD2", 1)
label("Cohorts", 0.04, 0.785, size = 11, face = "bold", col = teal)
label("EOCRC  < 50 years", 0.13, 0.805, size = 9.5, face = "bold")
label("AOCRC  ≥ 50 years", 0.13, 0.765, size = 9.5)
label("TCGA COAD/READ    358 tumors    54 early  ·  304 average onset", 0.32, 0.805, size = 9.5)
label("HM450 methylation  +  RNA-seq", 0.32, 0.765, size = 9, col = "#5C6B72")
label("ORIEN Avatar    453 tumors    120 early  ·  333 average onset", 0.68, 0.805, size = 9.5)
label("RNA-seq, MSI-filtered  ·  16S on a paired subset", 0.68, 0.765, size = 9, col = "#5C6B72")

panels <- data.frame(
  x = c(0.185, 0.50, 0.815),
  fill = c("#FFF8EC", "#F3FAF8", "#FFF4F1"),
  border = c(gold, teal, coral)
)
for (i in seq_len(nrow(panels))) {
  round_box(panels$x[i], 0.45, 0.30, 0.50, panels$fill[i], panels$border[i], 2)
}

# Panel A — clocks
label("A", 0.055, 0.66, size = 16, face = "bold", col = gold)
label("Epigenetic age", 0.085, 0.66, size = 14, face = "bold")
label("DNA methylation age is older\nin early-onset tumors.", 0.055, 0.585, size = 10.5)

# Three small clocks
clock_x <- c(0.09, 0.185, 0.28)
clock_lab <- c("Horvath", "Hannum", "PhenoAge")
for (i in seq_along(clock_x)) {
  grid.circle(
    x = clock_x[i], y = 0.455, r = 0.038,
    gp = gpar(fill = "white", col = gold, lwd = 2)
  )
  ang <- c(2.4, 0.9, 0.2)[i]
  grid.segments(
    x0 = clock_x[i], y0 = 0.455,
    x1 = clock_x[i] + 0.028 * cos(ang),
    y1 = 0.455 + 0.028 * sin(ang),
    gp = gpar(col = ink, lwd = 1.8, lineend = "round")
  )
  label(clock_lab[i], clock_x[i], 0.395, just = "center", size = 8, col = "#5C6B72")
}
label("+12 years", 0.185, 0.335, just = "center", size = 16, face = "bold", col = gold)
label("mean DNAm-age acceleration\nversus average-onset CRC", 0.185, 0.275, just = "center", size = 8.5, col = "#5C6B72")
label("CREB  ·  GPCR  ·  phagosome  ·  S100", 0.185, 0.225, just = "center", size = 8.5, face = "bold", col = ink)

# Panel B — microbes
label("B", 0.37, 0.66, size = 16, face = "bold", col = teal)
label("Intratumoral microbes", 0.40, 0.66, size = 14, face = "bold")
label("{exotic} counts non-human reads\nin bulk tumor RNA-seq.", 0.37, 0.585, size = 10.5)

# Three dataset glyphs, drawn as microbes rather than bars.
glyph_x <- c(0.43, 0.50, 0.57)
glyph_col <- c(teal, "#3D8B6E", coral)
for (i in seq_along(glyph_x)) {
  grid.roundrect(
    x = glyph_x[i], y = 0.485, width = 0.055, height = 0.028,
    r = unit(0.08, "inches"),
    gp = gpar(fill = glyph_col[i], col = NA)
  )
  grid.circle(
    x = glyph_x[i] - 0.012, y = 0.485, r = 0.006,
    gp = gpar(fill = "white", col = NA)
  )
}
label("TCGA", 0.43, 0.43, just = "center", size = 8, col = "#5C6B72")
label("ORIEN", 0.50, 0.43, just = "center", size = 8, col = "#5C6B72")
label("16S", 0.57, 0.43, just = "center", size = 8, col = "#5C6B72")
label("No shared enrichment", 0.50, 0.325, just = "center", size = 13, face = "bold", col = teal)
label("No taxon was higher in early-onset\ntumors across TCGA, ORIEN, and 16S.\nSite did not structure the community.", 0.50, 0.25, just = "center", size = 8.7, col = "#5C6B72")

# Panel C — immunity
label("C", 0.685, 0.66, size = 16, face = "bold", col = coral)
label("Microbe–immune axis", 0.715, 0.66, size = 14, face = "bold")
label("Cell fractions are similar.\nCorrelations are not.", 0.685, 0.585, size = 10.5)

# Mini correlation tiles: AO muted, EO stronger
tile_cols_ao <- colorRampPalette(c("#F6E7E2", "#E7B2A4"))(8)
tile_cols_eo <- colorRampPalette(c("#F6D2C8", "#D4654A"))(8)
set.seed(25)
ao_vals <- sample(tile_cols_ao, 12, replace = TRUE)
eo_vals <- sample(tile_cols_eo, 12, replace = TRUE)
for (k in seq_len(12)) {
  row <- (k - 1) %% 4
  col <- (k - 1) %/% 4
  grid.rect(
    x = 0.73 + col * 0.028, y = 0.50 - row * 0.032,
    width = 0.024, height = 0.026,
    gp = gpar(fill = ao_vals[k], col = "white", lwd = 0.4)
  )
  grid.rect(
    x = 0.84 + col * 0.028, y = 0.50 - row * 0.032,
    width = 0.024, height = 0.026,
    gp = gpar(fill = eo_vals[k], col = "white", lwd = 0.4)
  )
}
label("AOCRC", 0.758, 0.555, just = "center", size = 8, col = "#5C6B72")
label("EOCRC", 0.868, 0.555, just = "center", size = 8, face = "bold", col = coral)
label("Larger positive correlations", 0.815, 0.33, just = "center", size = 12.5, face = "bold", col = coral)
label("Activated mast cells across microbes.\nFusobacterium with neutrophils,\nlarger in early-onset tumors.", 0.815, 0.255, just = "center", size = 8.7, col = "#5C6B72")

# Conclusion
grid.roundrect(
  x = 0.5, y = 0.085, width = 0.96, height = 0.10,
  r = unit(0.12, "inches"),
  gp = gpar(fill = ink, col = NA)
)
label(
  "Early-onset tumors look epigenetically older, and their microbes engage immune cells more strongly,",
  0.5, 0.105, just = "center", size = 11, face = "bold", col = "white"
)
label(
  "without a microbe that consistently marks early-onset disease.",
  0.5, 0.065, just = "center", size = 11, col = "#F4C95D"
)

dev.off()
message("Wrote figures/hex_exogieo.png and figures/graphical_abstract.png")
