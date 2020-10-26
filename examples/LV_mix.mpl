# Example from Section 5.3
read "../ComputeIdentifiableFunctions.mpl";

# We are interested in the following model with extra assumption that f is known parameter
# different for different experiments
model := [
    diff(x1(t), t) - (a * x1(t) - b * x1(t) * x2(t)),
    diff(x2(t), t) - (-c * x2(t) + d * x1(t) * x2(t)),
    y(t) - (e * x1(t) + f * x2(t))
]:

# First we compute the bound and IO-equations considering f just as a parameter
result := MultiExperimentIdentifiableFunctions(model, simplified_generators=true):
printf("If f is considered as a parameter, the bound from the theorem will be %a\n", result[1]):
printf("The coefficients of the only input-output equation have the following degrees in f: %a\n", map(p -> degree(p, f), result[2][1])):
# In the notation of Remark 22
ns := sort(map(p -> degree(p, f), result[2][1])):
n := ns[1]:
printf("Remark 22 implies that the number of experiments for the case when f is a known parameter different in different experiments is %a\n", max(n + min(ns[2...]), max(ns[2..]))):
