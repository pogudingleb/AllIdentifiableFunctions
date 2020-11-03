# Example from Section 5.2 of the paper

read "../ComputeIdentifiableFunctions.mpl";

model := [
     diff(xa(t), t) + k1 * xa(t),
     diff(xb(t), t) - k1 * xa(t) + k2 * xb(t),
     diff(xc(t), t) - k2 * xb(t),
     diff(eA(t), t),
     diff(eC(t), t),
     y1(t) - eA(t) * xa(t) - eB * xb(t) - eC(t) * xc(t),
     y2(t) - xc(t),
     y3(t) - eA(t), 
     y4(t) - eC(t)
]:

# Computing function identifiable from a single experiment
printf("Single-experiment identifiable functions: %a\n", SingleExperimentIdentifiableFunctions(model)):

me := MultiExperimentIdentifiableFunctions(model, simplified_generators=true):
printf("The multi-experiment identifiable functions are %a; This is not the same as for a single experiment!\n", me[3]):
printf("The bound from the theorem is %a showing that %a experiments are sufficient to identify these functions\n", me[1], me[1]):
