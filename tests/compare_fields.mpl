read "../ComputeIdentifiableFunctions.mpl":

cases := [
    [
        [[a / b, b], [a, b]], true
    ],
    [
        [[a * b, a + b], [a, b]], false
    ],
    [
        [[a + b^2 - c^2, b * c, b / c], [a, b^2, c^2, b / c]], true
    ]
];

num_passed := 0:
num_failed := 0:

for case in cases do
    input := case[1]:
    correct := case[2]:
    if CompareFields(op(input)) = correct then
        printf("PASSED");
        num_passed := num_passed + 1:
    else
        printf("FAILED");
        num_failed := num_failed + 1:
    end if:
end do:

printf("Passed: %a, failed %a \n", num_passed, num_failed);
