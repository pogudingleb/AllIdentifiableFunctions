read "../ComputeIdentifiableFunctions.mpl";

cases := [
    [
        [
            diff(x1(t), t) = (a * x1(t) - b * x1(t) * x2(t) + u(t)),
            diff(x2(t), t) = (-c * x2(t) + d * x1(t) * x2(t)),
            y(t) = x1(t)
        ], {a, c, d}
    ],
    [
        [
            diff(xa(t), t) = -k1 * xa(t),
            diff(xb(t), t) = k1 * xa(t) - k2 * xb(t),
            diff(xc(t), t) = k2 * xb(t),
            diff(eA(t), t) = 0,
            diff(eC(t), t) = 0,
            y1(t) = eA(t) * xa(t) + eB * xb(t) + eC(t) * xc(t),
            y2(t) = xc(t),
            y3(t) = eA(t), 
            y4(t) = eC(t)
        ], {k1 + k2, k1 * k2}
    ] 
]:

num_passed := 0:
num_failed := 0:

for case in cases do
    input := case[1]:
    correct := case[2]:
    result := SingleExperimentIdentifiableFunctions(input, infolevel=0):
    print(result):
    if CompareFields(result, correct) then
        printf("PASSED\n");
        num_passed := num_passed + 1:
    else
        num_failed := num_failed + 1:
        printf("FAILED\n");
        printf("Expected %a, got %a", correct, result);
    end if:
end do:

printf("Passed: %a, failed %a \n", num_passed, num_failed);
