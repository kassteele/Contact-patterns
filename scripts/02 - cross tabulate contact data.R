#
# Read data
#

# Read
load(file = "results/polymod.data.bin")
load(file = "results/pop.data.bin")

#
# Cross tabulation
#

# Cross tabulate number of participants for each combination of age and sex
# t = number of participants
polymod.part.tab <- with(subset(polymod.data, subset = !duplicated(local.id)),
  as.data.frame(table(part.age = part.age, part.sex = part.sex),
    responseName = "t"))

# Cross tabulate number of contacts for each combination of age and sex
# y = number of participants x number of contacts
polymod.cont.tab <- with(polymod.data,
  as.data.frame(table(part.age = part.age, cont.age = cont.age, part.sex = part.sex, cont.sex = cont.sex),
    responseName = "y"))

# Make age an integer again (was factor in previous section)
polymod.part.tab <- within(polymod.part.tab, {
  part.age <- as.numeric(as.character(part.age))
})
polymod.cont.tab <- within(polymod.cont.tab, {
  part.age <- as.numeric(as.character(part.age))
  cont.age <- as.numeric(as.character(cont.age))
})

# Merge polymod.cont.tab, polymod.part.tab and pop.data
# polymod.tab is the tabulated version of polymod.data
polymod.tab <- merge(polymod.part.tab, merge(polymod.cont.tab, pop.data, all = TRUE), all = TRUE)

# Reorder polymod.tab
polymod.tab <- polymod.tab[c("part.age", "cont.age", "part.sex", "cont.sex", "y", "t", "w")]
polymod.tab <- with(polymod.tab, polymod.tab[order(cont.sex, part.sex, cont.age, part.age), ])

# Calculate denominator U = number of participants t x population w
# Divide by 1e6 for better numerical properties (inflates contact rate c)
# If t or w = 0, then set U = 1 and y = NA (record will not contribute to likelihood)
polymod.tab <- within(polymod.tab, {
  U <- t*w/1e6
  y <- ifelse(U == 0, yes = NA, no = y)
  U <- ifelse(U == 0, yes = 1, no = U)
})

#
# Save result
#

# Save
save(polymod.tab, file = "results/polymod.tab.bin")
