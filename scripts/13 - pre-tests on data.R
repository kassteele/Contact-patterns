#
# Pre-tests on data
#

#
# Init
#

# Load packages
library(RColorBrewer)

# Load data
load(file = "results/polymod.data.bin")

#
# Test if number of participants depends on day of the week
#

tab <- with(subset(polymod.data, !duplicated(local.id)), table(part.day))
chisq.test(tab) # Conclusion: no evidence for difference (NL)

#
# Crude number of contacts per age group and sex
#

tmp <- within(polymod.tab, m <- y/t)
tmp.agg <- with(tmp, aggregate(m ~ part.age + part.sex, FUN = sum))
tmp.agg <- within(tmp.agg, part.age.cat <- cut(as.numeric(as.character(part.age)), breaks = seq(0, 81, by = 10), include.lowest = TRUE))

# Settings
cols.line <- brewer.pal(n = 8, "Accent")[c(5, 6)] # Color for males and females: lines
cols.bord <- adjustcolor(cols.line, alpha = 0.5)

# Plot
par(mar = c(4, 4, 0.3, 0.3), mgp = c(2.5, 1, 0), xaxs = "r", yaxs = "r", cex.axis = 0.9, bty = "n")
with(tmp.agg, boxplot(m ~ part.sex + part.age.cat,
  xlab = "Age category", ylab = "Reported mean number of contacts per day",
  ylim = c(0, 31),
  at = (1:16)+c(0.1, -0.1),
  col = cols.bord, border = cols.line,
  xaxt = "n", yaxt = "n"))
axis(1, at = seq(1.5, 15.5, by = 2), labels = levels(tmp.agg$part.age.cat), lwd = 0.5)
axis(2, lwd = 0.5)
par(bty = "o"); box(lwd = 0.5)
legend("topright", legend = c("Males", "Females"), fill = cols.bord, border = cols.line, bty = "n", cex = 0.9)

#
# Distribution of total number of contacts
#

# Check Negative Binomial assumption
tmp <- within(polymod.tab, {
  part.age.cat <- cut(as.numeric(as.character(part.age)), breaks = seq(0, 81, by = 10), include.lowest = TRUE)
  cont.age.cat <- cut(as.numeric(as.character(cont.age)), breaks = seq(0, 81, by = 10), include.lowest = TRUE)
})
Ey <- with(tmp, aggregate(y ~ part.age.cat + cont.age.cat + part.sex + cont.sex, FUN = mean, na.rm = TRUE)$y)
Vy <- with(tmp, aggregate(y ~ part.age.cat + cont.age.cat + part.sex + cont.sex, FUN = var, na.rm = TRUE)$y)
plot(Ey, sqrt(Vy))


