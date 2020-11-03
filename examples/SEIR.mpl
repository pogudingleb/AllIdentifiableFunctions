# Taken from
# N. Tuncer, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.2) with cumulative incidence observations
read "../ComputeIdentifiableFunctions.mpl";

# Prevalence observations
model := [
  diff(S(t), t) - b * S(t) * In(t) / N,
  diff(E(t), t) - b * S(t) * In(t) / N + nu * E(t),
  diff(In(t), t) - nu * E(t) + a * In(t),
  y(t) - In(t),
  y2(t) - N
]:

me := MultiExperimentIdentifiableFunctions(model, simplified_generators=true):
printf("The bound for the number of experiments is %a, this means that the fields of single-experiment and multi-experiment identifiable functions coincide, and they are equal to: %a\n", me[1], me[3]):

# Cummulative incidence observations
# using change of variables Ninv := 1 / N to simplify the computation
model := [
  diff(S(t), t) - b * S(t) * In(t) * Ninv,
  diff(E(t), t) - b * S(t) * In(t) * Ninv + nu * E(t),
  diff(In(t), t) - nu * E(t) + a * In(t),
  diff(Cu(t), t) - nu * E(t),
  y(t) - Cu(t),
  y2(t) - Ninv
]:

me := MultiExperimentIdentifiableFunctions(model, simplified_generators=true, infolevel=1):
printf("The bound for the number of experiments is %a, this means that the field of multi-experiment identifiable functions is the same, %a\n", me[1], me[3]):
