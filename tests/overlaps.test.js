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

test("tree searching works correctly for the general case", () => {
    let starts = [];
    let ends = [];
    for (var i = 0; i < 1000; i++) {
        let s = Math.floor(Math.random() * 1000);
        starts.push(s);
        ends.push(s + Math.floor(Math.random() * 20));
    }
    let tree = overlap.buildIntervalTree(starts, ends);

    for (var i = 0; i < 200; i++) {
        let qs = Math.floor(Math.random() * 1000);
        let qe = qs + Math.floor(Math.random() * 20);

        let results = overlap.queryIntervalTree(qs, qe, tree);
        results.sort((a, b) => a - b);
        let expected = slowReference(starts, ends, qs, qe);

//        console.log([
//            qs, 
//            qe, 
//            [expected.map(i => starts[i]), expected.map(i => ends[i])], 
//            [results.map(i => starts[i]), results.map(i => ends[i])], 
//        ]);

        expect(results).toEqual(expected);
    }
})

test("tree searching works correctly for zero-length intervals", () => {
    let starts = [];
    let ends = [];
    for (var i = 0; i < 100; i++) {
        let s = Math.floor(Math.random() * 50);
        starts.push(s);
        ends.push(s);
    }
    let tree = overlap.buildIntervalTree(starts, ends);

    for (var i = 0; i < 20; i++) {
        let q = Math.floor(Math.random() * 50);

        let results = overlap.queryIntervalTree(q, q, tree);
        results.sort((a, b) => a - b);
        let expected = slowReference(starts, ends, q, q);
        expect(results).toEqual(expected);
    }
})

test("tree searching works correctly with heavy overlaps", () => {
    let starts = [];
    let ends = [];
    for (var i = 0; i < 100; i++) {
        let s = Math.floor(Math.random() * 1000);
        starts.push(s);
        ends.push(s + Math.floor(50 + Math.random() * 50));
    }
    let tree = overlap.buildIntervalTree(starts, ends);

    for (var i = 0; i < 20; i++) {
        let qs = Math.floor(Math.random() * 1000);
        let qe = qs + Math.floor(50 + Math.random() * 50);

        let results = overlap.queryIntervalTree(qs, qe, tree);
        results.sort((a, b) => a - b);
        let expected = slowReference(starts, ends, qs, qe);
        expect(results).toEqual(expected);
    }
})

test("tree searching works correctly after slicing", () => {
    let starts = [];
    let ends = [];
    for (var i = 0; i < 100; i++) {
        let s = Math.floor(Math.random() * 1000);
        starts.push(s);
        ends.push(s + Math.floor(Math.random() * 20));
    }

    let slice = [0, 5, 9, 11, 14, 17, 22, 33, 41, 42, 44, 49, 57];
    let tree = overlap.buildIntervalTree(starts, ends, { slice: slice });
    let ref = overlap.buildIntervalTree(bioc.SLICE(starts, slice), bioc.SLICE(ends, slice));

    for (var i = 0; i < 200; i++) {
        let qs = Math.floor(Math.random() * 1000);
        let qe = qs + Math.floor(Math.random() * 20);

        let results = overlap.queryIntervalTree(qs, qe, tree);
        let expected = overlap.queryIntervalTree(qs, qe, ref).map(j => slice[j]);
        expect(results).toEqual(expected);
    }
})

test("tree searching skips correctly for out-of-range queries", () => {
    let tree = overlap.buildIntervalTree([10, 20, 30], [30, 40, 50]);

    expect(overlap.queryIntervalTree(0, 5, tree)).toEqual([]);
    expect(overlap.queryIntervalTree(60, 65, tree)).toEqual([]);

    expect(overlap.queryIntervalTree(10, 10, tree)).toEqual([0]);
    expect(overlap.queryIntervalTree(9, 10, tree)).toEqual([]);
    expect(overlap.queryIntervalTree(9, 9, tree)).toEqual([]);

    expect(overlap.queryIntervalTree(50, 50, tree)).toEqual([]);
    expect(overlap.queryIntervalTree(49, 50, tree)).toEqual([2]);
    expect(overlap.queryIntervalTree(49, 49, tree)).toEqual([2]);
})

