# Example from Section 5.1 of the paper

read "../ComputeIdentifiableFunctions.mpl";

model := [
    diff(x1(t), t) - (a * x1(t) - b * x1(t) * x2(t)),
    diff(x2(t), t) - (-c * x2(t) + d * x1(t) * x2(t) + e * u(t)),
    y(t) - x1(t)
]:

# Computing function identifiable from a single experiment
printf("Single-experiment identifiable functions: %a\n", SingleExperimentIdentifiableFunctions(model)):

me := MultiExperimentIdentifiableFunctions(model, simplified_generators=true):
printf("The bound for the number of experiments is %a, this means that the field of multi-experiment identifiable functions is the same, %a\n", me[1], me[3]):
