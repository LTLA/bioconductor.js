import * as bioc from "../src/index.js";

test("presplitFactor works as expected", () => {
    var factor = ["A", "B", "C", "A", "C", "B", "D"];
    let split = bioc.presplitFactor(factor);
    expect(split.A).toEqual([0, 3]);
    expect(split.B).toEqual([1, 5]);
    expect(split.C).toEqual([2, 4]);
    expect(split.D).toEqual([6]);
})

test("which works as expected", () => {
    let truth = [ 0, 1, 1, 0, 0, 1, -1 ];
    expect(bioc.which(truth)).toEqual([1,2,5,6]);
    expect(bioc.which(truth, { not: true })).toEqual([0, 3, 4]);
})
