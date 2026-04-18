library(ggplot2)
library(patchwork)
library(scales)

# ── Data ──────────────────────────────────────────────────────────────────────

ne <- data.frame(
  Population = factor(c("DICF","KGCF","LICF","SICF"),
                      levels = c("DICF","KGCF","LICF","SICF")),
  Ne   = c(169.8, 798.7,  86.0, 103.8),
  low  = c(151.9, 340.1,  78.1,  56.0),
  high = c(187.8, 1257.4, 93.9, 151.6),
  n    = c(22, 19, 16, 13)
)

sens <- data.frame(
  Population = factor(
    rep(c("DICF","KGCF","LICF","SICF"), each = 3),
    levels = c("DICF","KGCF","LICF","SICF")
  ),
  Pcrit    = rep(c(0.05, 0.02, 0.01), 4),
  Ne       = c(169.8, 212.3, 212.3,
               798.7, 5545.1, 5545.1,
               86.0,  106.8, 106.8,
               103.8, 103.8, 103.8),
  reliable = c(TRUE,  TRUE,  TRUE,
               TRUE,  FALSE, FALSE,
               TRUE,  TRUE,  TRUE,
               TRUE,  TRUE,  TRUE)
)

cols <- c(DICF = "#1D9E75", KGCF = "#378ADD",
          LICF = "#7F77DD", SICF = "#D85A30")

# ── Panel A ───────────────────────────────────────────────────────────────────
# Fix 1: Split into two y-axis ranges using facet or annotation
# Best solution: use a broken axis via two separate plots stacked,
# OR simply show KGCF separately with a note and cap y at 300 for clarity.
# Here we use a clean approach: show all populations but use a sqrt scale
# so small Ne values are visible, and annotate KGCF value explicitly.

pA <- ggplot(ne, aes(x = Population, y = Ne, color = Population)) +
  # Shaded band for LICF/SICF range to aid reading
  annotate("rect", xmin = 0.5, xmax = 4.5,
           ymin = 0, ymax = 200,
           fill = "grey95", alpha = 0.5) +
  geom_hline(yintercept = c(100, 200, 500, 800),
             color = "grey88", linewidth = 0.4) +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    size      = 0.7,
    linewidth = 1.2,
    fatten    = 3
  ) +
  # Annotate KGCF upper CI explicitly since it goes off scale
  geom_text(
    data = subset(ne, Population == "KGCF"),
    aes(y = high + 35, label = paste0("CI: ", low, "–", high)),
    size  = 2.8,
    color = "#378ADD",
    fontface = "plain"
  ) +
  # Annotate Ne values next to each point
  geom_text(
    aes(y = Ne, label = Ne),
    nudge_x = 0.28,
    size     = 2.8,
    fontface = "plain",
    show.legend = FALSE
  ) +
  # Sample size below x axis
  geom_text(
    aes(y = -55, label = paste0("n=", n)),
    size  = 2.6,
    color = "grey50"
  ) +
  scale_color_manual(values = cols) +
  scale_y_continuous(
    limits = c(-80, 1350),
    breaks = c(0, 100, 200, 500, 800),
    labels = comma,
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x     = NULL,
    y     = expression(italic(N)[e]),
    title = "(a)  Contemporary Ne — Pcrit = 0.05"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.title         = element_text(size = 10, face = "plain", hjust = 0),
    axis.text.x        = element_text(face = "bold", size = 10),
    axis.text.y        = element_text(size = 9),
    plot.margin        = margin(t = 8, r = 12, b = 20, l = 8)
  )

# ── Panel B ───────────────────────────────────────────────────────────────────
# Fix 2: Move legend outside plot (bottom), fix line direction,
#         add clear point shapes for reliable vs unreliable

pB <- ggplot(sens, aes(x = Pcrit, y = Ne,
                        color    = Population,
                        group    = Population,
                        linetype = Population)) +
  geom_line(linewidth = 1.0) +
  geom_point(aes(shape = reliable, size = reliable)) +
  # Annotate KGCF inflation with arrow label
  annotate("text",
           x = 0.025, y = 4000,
           label = "KGCF inflates\nat lower Pcrit",
           size = 2.8, color = "#378ADD",
           hjust = 0, lineheight = 0.9) +
  annotate("segment",
           x = 0.025, xend = 0.021,
           y = 3500, yend = 5545,
           color = "#378ADD",
           arrow = arrow(length = unit(0.15, "cm")),
           linewidth = 0.5) +
  scale_color_manual(values = cols, name = "Population") +
  scale_linetype_manual(
    values = c(DICF = "solid", KGCF = "dashed",
               LICF = "solid", SICF = "solid"),
    name = "Population"
  ) +
  scale_shape_manual(
    values = c("TRUE" = 19, "FALSE" = 1),
    guide  = "none"
  ) +
  scale_size_manual(
    values = c("TRUE" = 2.5, "FALSE" = 2.5),
    guide  = "none"
  ) +
  scale_x_reverse(
    breaks = c(0.05, 0.02, 0.01),
    labels = c("0.05", "0.02", "0.01")
  ) +
  scale_y_log10(
    breaks = c(100, 500, 1000, 5000),
    labels = comma,
    limits = c(70, 8000)
  ) +
  annotation_logticks(sides = "l", size = 0.3,
                      short = unit(0.1, "cm"),
                      mid   = unit(0.15, "cm"),
                      long  = unit(0.2, "cm")) +
  labs(
    x     = "Pcrit threshold",
    y     = expression(italic(N)[e]~"(log"[10]~"scale)"),
    title = "(b)  Pcrit sensitivity"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 9, face = "plain"),
    legend.text        = element_text(size = 9),
    legend.key.width   = unit(1.4, "cm"),
    legend.key.height  = unit(0.4, "cm"),
    legend.margin      = margin(t = 0),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.title         = element_text(size = 10, face = "plain", hjust = 0),
    axis.text          = element_text(size = 9),
    plot.margin        = margin(t = 8, r = 8, b = 8, l = 12)
  ) +
  guides(
    color    = guide_legend(nrow = 1, override.aes = list(linewidth = 1)),
    linetype = guide_legend(nrow = 1)
  )

# ── Combine and save ──────────────────────────────────────────────────────────

combined <- pA + pB +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    caption = paste0(
      "Open circles = unreliable estimates (wide CI). ",
      "Solid circles = reliable estimates. ",
      "Dashed line (KGCF) indicates sensitivity to Pcrit threshold.\n",
      "Ne estimated using LD method (Waples 2006). ",
      "95% CI by block jackknife (200 blocks)."
    ),
    theme = theme(
      plot.caption = element_text(size = 8, color = "grey50",
                                  hjust = 0, lineheight = 1.3)
    )
  )

ggsave("/nesi/nobackup/uoo04082/Amphi/ne_results/Figure_Ne_multipanel.pdf",
       combined, width = 8, height = 4.2, device = cairo_pdf)

ggsave("/nesi/nobackup/uoo04082/Amphi/ne_results/Figure_Ne_multipanel.png",
       combined, width = 8, height = 4.2, dpi = 300)

cat("Saved:\n")
cat("  ne_results/Figure_Ne_multipanel.pdf\n")
cat("  ne_results/Figure_Ne_multipanel.png\n")
