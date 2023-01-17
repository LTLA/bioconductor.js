import * as bioc from "../src/index.js";
import * as overlap from "../src/overlap-utils.js";

test("position to rank conversion works correctly", () => {
    let start = [5,1,20,30,40];
    let end = [6,4,25,40,44];
    let conversion = overlap.convertPositionToRank(start, end);

    let sorted = [ 1, 4, 5, 6, 20, 25, 30, 40, 44 ];
    expect(conversion.rank2position).toEqual(sorted);
    expect(Array.from(conversion.startRanks).map(i => sorted[i])).toEqual(start);
    expect(Array.from(conversion.endRanks).map(i => sorted[i])).toEqual(end);
})

test("position to rank conversion works correctly with slicing", () => {
    let start = [5,1,5,30,398,188,40,50];
    let end = [99,4,25,40, 400, 250, 44, 99];
    let slice = [0, 2, 4, 6];
    let conversion = overlap.convertPositionToRank(start, end, { slice: slice });

    let sorted = [ 5, 25, 40, 44, 99, 398, 400 ]
    expect(conversion.rank2position).toEqual(sorted);
    expect(Array.from(conversion.startRanks).map(i => sorted[i])).toEqual(bioc.SLICE(start, slice));
    expect(Array.from(conversion.endRanks).map(i => sorted[i])).toEqual(bioc.SLICE(end, slice));
})

function slowReference(ref_start, ref_end, query_start, query_end) {
    let overlaps = [];
    for (var r = 0; r < ref_start.length; r++) {
        if ((ref_start[r] < query_end && query_start < ref_end[r]) || ref_start[r] == query_start) { // handle zero-length intervals correctly.
            overlaps.push(r);
        }
    }
    return overlaps;
}

test("tree searching works correctly for simple cases", () => {
    let starts = [];
    let ends = [];
    for (var i = 0; i < 100; i++) {
        let s = Math.floor(Math.random() * 1000);
        starts.push(s);
        ends.push(s + Math.floor(Math.random() * 20));
    }
    let tree = overlap.buildIntervalTree(starts, ends);

    for (var i = 0; i < 20; i++) {
        let qs = Math.floor(Math.random() * 1000);
        let qe = qs + Math.floor(Math.random() * 20);

        let results = overlap.queryIntervalTree(qs, qe, tree);
        results.sort((a, b) => a - b);
        let expected = slowReference(starts, ends, qs, qe);
        expect(results).toEqual(expected);
    }
})

