#
# Read POLYMOD and population data
#

# Read data
part.data <- read.delim(file = "data/contact data/polymod/Participants_20071108_NL.txt")
cont.data <- read.delim(file = "data/contact data/polymod/Contacts_20071108_NL.txt")
pop.data <- read.delim(file = "data/population data/CBS age sex pop 2007.txt")

# Select relevant columns
part.data  <- subset(part.data, select = c("local_id", "participant_age", "participant_gender", "dayofweek"))
cont.data  <- subset(cont.data, select = c("local_id", "cnt_age_l", "cnt_age_r", "cnt_sex", "cnt_touch"))

# Optional: subset cnt_touch. 1 = Yes, 2 = No
# cont.data <- subset(cont.data, subset = cnt_touch == 1 & !is.na(cnt_touch))

#
# Data pre-processing
#

# Clean pop.data
pop.data <- within(pop.data, {
  # Rename pop to w
  w <- pop
  # Rename age to cont.age (needed for merge with polymod data later)
  cont.age <- age
  # Rename sex to cont.sex and reorder levels (also needed later)
  cont.sex <- factor(sex, levels = c("man", "vrouw"), labels = c("male", "female"))
  # Remove old variabels
  rm(age, sex, pop)
})

# Clean part.data
names(part.data) <- gsub(names(part.data), pattern = "_", replacement = "." )
names(part.data) <- gsub(names(part.data), pattern = "participant", replacement = "part")
part.data <- within(part.data, {
  # Rename part.gender to part.sex and relabel levels
  part.sex <- factor(part.gender, levels = c("M", "F"), labels = c("male", "female"))
  # Give dayofweek informative labels and call it part.day
  part.day <- factor(dayofweek, levels = 0:6, labels = c("sunday", "monday", "tuesday", "wednesday", "thursday", "friday", "saturday"))
  # Remove old variabels
  rm(part.gender, dayofweek)
})

# Modify cont.data
# Replace _ by . and cnt by cont
names(cont.data) <- gsub(names(cont.data), pattern = "_", replacement = ".")
names(cont.data) <- gsub(names(cont.data), pattern = "cnt", replacement = "cont")
cont.data <- within(cont.data, {
  # Relabel levels in cont.sex
  cont.sex <- factor(cont.sex, levels = c("M", "F"), labels = c("male", "female"))
  # Make cont.age.l numeric and make NA if onbekend
  cont.age.l <- as.numeric(ifelse(cont.age.l == "onbekend", yes = NA, no = as.character(cont.age.l)))
})

#
# Construct polymod.data for analysis
#

# Construct contact age in cont.data
# Two problems: 1) Contact age is sometimes given as a range
#               2) Digit preferencing in contact age

# Set seed for reproducibility
set.seed(1)

# 1) Take care of the age range
cont.data <- within(cont.data, {
  # Start with empty cont.age
  cont.age <- NA
  # All ok? Then cont.age is cont.age.l
  ix <- (!is.na(cont.age.l) & is.na(cont.age.r)) | (!is.na(cont.age.l) & !is.na(cont.age.r) & cont.age.r<=cont.age.l)
  cont.age[ix] <- cont.age.l[ix]
  # Missing cont.age.l and cont.age.r reported? Then cont.age is cont.age.r
  ix <- is.na(cont.age.l) & !is.na(cont.age.r)
  cont.age[ix] <- cont.age.r[ix]
  # Age given in a range? Then sample uniformly from range
  ix <- !is.na(cont.age.l) & !is.na(cont.age.r) & cont.age.r > cont.age.l
  cont.age[ix] <- mapply(FUN = function(l, r) sample(l:r, size = 1), l = cont.age.l[ix], r = cont.age.r[ix])
  # Remove old variabels
  rm(cont.age.l, cont.age.r, ix)
})

# 2) Correction for digit preferencing in cont.age
# Extract age from cont.data
age <- cont.data$cont.age
# age.tab is a table with number of contacts for the non-pref ages
age.tab <- table(age[age>=20 & !is.element(age, seq(20, 80, 5))])
# Fit smooth spline through number of contacts for the non-pref ages
age.smt <- smooth.spline(x = as.numeric(names(age.tab)), y = age.tab, cv = TRUE)
# For each pref age, calculate the difference between the peak and the spline
# Then redistribute the difference over the neighbouring ages -2:2
for (age.i in seq(from = 20, to = 80, by = 5)) {
  ix1 <- age == age.i & !is.na(age)
  n.obs <- sum(ix1)                             # Observed number at age i
  n.exp <- round(predict(age.smt, x = age.i)$y) # Expected number at age i
  ix2 <- sample(1:n.obs, size = n.obs - n.exp)
  correction <- try(sample(-2:2, size = n.obs - n.exp, replace = TRUE))
  age[ix1][ix2] <- age[ix1][ix2] + correction
}
# Put it pack in cont.data
cont.data$cont.age <- age

# Remove unused variables
rm(age, age.smt, age.tab, age.i, correction, ix1, ix2, n.exp, n.obs)

#
# Construct polymod.data from part.data and cont.data
#

# Set age range
age <- 0:80
n <- length(age)

# Merge part.data and cont.data into polymod.data
polymod.data <- merge(part.data, cont.data, all = TRUE)

# Count number of participants and contacts
cat("Before\n")
print(length(unique(polymod.data$local.id))) # Participants
print(length(polymod.data$local.id))         # Contacts

# Omit missing age, sex and age > 80
polymod.data <- droplevels(subset(polymod.data,
  !(is.na(part.age) | is.na(part.sex) | part.age > max(age) | is.na(part.day) |
      is.na(cont.age) | is.na(cont.sex) | cont.age > max(age))))

# Make age a categorical variable (needed for cross tabulation in next section)
polymod.data <- within(polymod.data, {
  part.age <- factor(part.age, levels = age)
  cont.age <- factor(cont.age, levels = age)
})

# Count number of participants and contacts again
cat("After\n")
print(length(unique(polymod.data$local.id))) # Participants
print(length(polymod.data$local.id))           # Contacts

# Also in pop.data, keep age <= max(age)
pop.data <- subset(pop.data, cont.age <= max(age))

# Reorder pop.data
pop.data <- with(pop.data, pop.data[order(cont.sex, cont.age), ])

#
# Save results as R binaries
#

# Save
save(polymod.data, file = "results/polymod.data.bin")
save(pop.data, file = "results/pop.data.bin")
