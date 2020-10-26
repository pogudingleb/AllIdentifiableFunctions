read "../ComputeIdentifiableFunctions.mpl":

cases := [
    [
        [[a * b, b], [a]],
        [1 / a]
    ],
    [
        [[a + b, a * b], [b]],
        []
    ],
    [
        [[a * b + c, b], [a, c]],
        []
    ],
    [
        [[a * b, c * b, d * b], [a, c, d]],
        [a / c, a / d]
    ],
    [
        [[a + b + c, a^2 + b^2 + c^2, a^3 + b^3 + c^3, c], [a, b]],
        [a + b, a * b]
    ],
    [
        [[a + b, a * b, x, a * x + b * c], [a, b, c]],
        [a + b, a * b]
    ]
]:

num_passed := 0:
num_failed := 0:

for case in cases do
    input := case[1]:
    correct := case[2]:
    if IdealsEq(FieldIntersection(op(input)), FieldToIdeal(correct)) then
        printf("PASSED\n");
        num_passed := num_passed + 1:
    else
        printf("FAILED\n");
        num_failed := num_failed + 1:
        print("Expected: ", correct);
        print("Got: ", FieldIntersection(op(input)));
    end if:
end do:

printf("Passed: %a, failed %a \n", num_passed, num_failed);
