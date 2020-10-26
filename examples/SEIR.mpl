# Taken from
# N. Tuncer, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.2) with cumulative incidence observations
read "../ComputeIdentifiableFunctions.mpl";

model := [
  diff(S(t), t) - b * S(t) * In(t) / N,
  diff(E(t), t) - b * S(t) * In(t) / N + nu * E(t),
  diff(In(t), t) - nu * E(t) + a * In(t),
  y(t) - In(t),
  y2(t) - N
]:

# Computing function identifiable from a single experiment
printf("Single-experiment identifiable functions: %a\n", SingleExperimentIdentifiableFunctions(model)):

me := MultiExperimentIdentifiableFunctions(model, simplified_generators=true):
printf("The bound for the numer of experiments is %a, this means that the field of multi-experiment identifiable functions is the same, %a\n", me[1], me[3]):
