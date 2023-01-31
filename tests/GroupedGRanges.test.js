import * as bioc from "../src/index.js";
import * as utils from "./utils.js";

export function spawn_GRanges_list(N) {
    let output = [];
    for (var i = 0; i < N; i++) {
        output.push(utils.spawn_random_GRanges(Math.ceil(Math.random() * 10)));
    }
    return output;
}

test("constructing a GroupedGRanges works", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    expect(ggr.numberOfGroups()).toEqual(10);
    expect(ggr.ranges().seqnames()).toEqual(bioc.COMBINE(grl_raw.map(x => x.seqnames())));
    expect(ggr.ranges().start()).toEqual(bioc.COMBINE(grl_raw.map(x => x.start())));
    expect(ggr.ranges().width()).toEqual(bioc.COMBINE(grl_raw.map(x => x.width())));

    expect(ggr.group(0).start()).toEqual(grl_raw[0].start());
    expect(ggr.group(9).start()).toEqual(grl_raw[9].start());

    // Passes along arguments to the Vector constructor.
    let ggr2 = new bioc.GroupedGRanges(grl_raw, {
        metadata: { foo: 1 },
        names: [ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J" ],
        elementMetadata: new bioc.DataFrame({ foo: utils.spawn_random_vector(10) })
    });

    expect(ggr2.metadata().foo).toBe(1);
    expect(ggr2.names().length).toBe(10);
    expect(ggr2.elementMetadata().column("foo").length).toBe(10);

    // Fails if the array contains non-GRanges.
    expect(() => new bioc.GroupedGRanges([1,2,3])).toThrow("must either be a 'GRanges'");
})

test("constructing a GroupedGRanges works with runLengths", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    let grl_len = grl_raw.map(bioc.LENGTH);
    let grl_com = bioc.COMBINE(grl_raw);
    let ggr2 = new bioc.GroupedGRanges(grl_com, { rangeLengths: grl_len });

    expect(ggr.rangeStarts()).toEqual(ggr2.rangeStarts());
    expect(ggr.rangeLengths()).toEqual(ggr2.rangeLengths());
    expect(ggr.ranges().seqnames()).toEqual(ggr2.ranges().seqnames());
    expect(ggr.ranges().start()).toEqual(ggr2.ranges().start());
    expect(ggr.ranges().width()).toEqual(ggr2.ranges().width());

    // Fails for inconsistent lengths.
    expect(() => new bioc.GroupedGRanges(grl_com)).toThrow("rangeLengths");
    expect(() => new bioc.GroupedGRanges(grl_com, { rangeLengths: [1,2,3] })).toThrow("rangeLengths");
    expect(() => new bioc.GroupedGRanges(grl_com, { rangeLengths: [1,1,-1,1,1,1,1,1,1,1] })).toThrow("negative");
    expect(() => new bioc.GroupedGRanges(grl_com, { rangeLengths: [1,1,1,1,1,1,1,1,1,1] })).toThrow("sum of 'rangeLengths'");
})

test("setting the inner ranges works as expected", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    let rr = bioc.CLONE(ggr.ranges());
    let str = new Int8Array(bioc.LENGTH(rr));
    str.fill(-1);
    rr.$setStrand(str);

    let ggr2 = ggr.setRanges(rr);
    expect(ggr.ranges().strand()[0]).toEqual(0); // doesn't mutate the original.
    expect(ggr2.ranges().strand()[0]).toEqual(-1);

    ggr.$setRanges(rr);
    expect(ggr.ranges().strand()[0]).toEqual(-1);

    // Errors
    expect(() => ggr.$setRanges(1)).toThrow("GRanges");
    expect(() => ggr.$setRanges(bioc.SLICE(rr, [0,1,2,3]))).toThrow("length");
})

test("setting an individual group works as expected", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    let current = ggr.group(5);
    let before = ggr.group(4);
    let after = ggr.group(6);

    let replacement = utils.spawn_random_GRanges(bioc.LENGTH(current) + 5);
    let ggr2 = ggr.setGroup(5, replacement);
    expect(ggr2.group(5).start()).toEqual(replacement.start());
    expect(ggr.group(5).start()).toEqual(current.start()); // doesn't mutate the original.

    ggr.$setGroup(5, replacement);
    expect(ggr.group(5).start()).toEqual(replacement.start());
    expect(ggr.group(5).start()).not.toEqual(current.start());

    // Groups before and after are preserved.
    expect(ggr.group(4).start()).toEqual(before.start());
    expect(ggr.group(6).start()).toEqual(after.start());
})

test("setting multiple different groups works as expected", () => {
    let grl_raw = spawn_GRanges_list(10);

    { 
        let ggr = new bioc.GroupedGRanges(grl_raw);
        let grl_raw2 = spawn_GRanges_list(10);
        for (var i = 9; i >= 0; i--) { // setting it in reverse, to check the sorting.
            ggr.$setGroup(i, grl_raw2[i]);
        }

        let ggr2 = new bioc.GroupedGRanges(grl_raw2);
        expect(ggr2.ranges().start()).toEqual(ggr.ranges().start());
        expect(ggr2.rangeStarts()).toEqual(ggr.rangeStarts());
        expect(ggr2.rangeLengths()).toEqual(ggr.rangeLengths());
    }

    // Trying again, only replacing even indices.
    {
        let ggr = new bioc.GroupedGRanges(grl_raw);
        let updated = [...grl_raw];
        let grl_raw2 = spawn_GRanges_list(5);
        for (var i = 0; i < 5; i++) {
            ggr.$setGroup(i * 2, grl_raw2[i]);
            updated[i * 2] = grl_raw2[i];
        }

        let ggr2 = new bioc.GroupedGRanges(updated);
        expect(ggr2.ranges().start()).toEqual(ggr.ranges().start());
        expect(ggr2.rangeStarts()).toEqual(ggr.rangeStarts());
        expect(ggr2.rangeLengths()).toEqual(ggr.rangeLengths());
    }

    // Checking for immutability.
    {
        let ggr = new bioc.GroupedGRanges(grl_raw);
        let copy = bioc.CLONE(ggr);
        let updated = [...grl_raw];
        let grl_raw2 = spawn_GRanges_list(5);

        for (var i = 0; i < 4; i++) { 
            let ggr2 = ggr.setGroup(i * 2 + 1, grl_raw2[i]);
            updated[i * 2 + 1] = grl_raw2[i];

            let ggr3 = new bioc.GroupedGRanges(updated);
            expect(ggr2.ranges().start()).toEqual(ggr3.ranges().start());
            expect(ggr2.rangeStarts()).toEqual(ggr3.rangeStarts());
            expect(ggr2.rangeLengths()).toEqual(ggr3.rangeLengths());

            expect(copy.ranges().start()).toEqual(ggr.ranges().start());
            expect(copy.rangeStarts()).toEqual(ggr.rangeStarts());
            expect(copy.rangeLengths()).toEqual(ggr.rangeLengths());
        }
    }
})

test("setting an individual group multiple times works as expected", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    let current = ggr.group(3);
    let before = ggr.group(0);
    let after = ggr.group(9);

    // Second $setGroup clobbers the first, as it should.
    let replacement1 = utils.spawn_random_GRanges(bioc.LENGTH(current) + 5);
    ggr.$setGroup(3, replacement1);

    let replacement2 = utils.spawn_random_GRanges(bioc.LENGTH(current) + 3);
    ggr.$setGroup(3, replacement2);

    expect(ggr.group(3).start()).toEqual(replacement2.start());
    expect(ggr.group(3).start()).not.toEqual(current.start());
    expect(ggr.group(3).start()).not.toEqual(replacement1.start());

    // Groups before and after are preserved.
    expect(ggr.group(0).start()).toEqual(before.start());
    expect(ggr.group(9).start()).toEqual(after.start());
})

test("LENGTH generic works as expected", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);
    expect(bioc.LENGTH(ggr)).toEqual(10);
})

test("SLICE generic works as expected with indices", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    let indices = [0, 2, 4, 6, 8];
    let sliced = bioc.SLICE(ggr, indices);
    expect(sliced.numberOfGroups()).toEqual(indices.length);
    for (var i = 0; i < indices.length; i++) {
        expect(sliced.group(i).start()).toEqual(grl_raw[indices[i]].start());
    }

    expect(sliced.rangeStarts()[0]).toBe(0);
    expect(Array.from(sliced.rangeLengths())).toEqual(indices.map(i => bioc.LENGTH(grl_raw[i])));
})

test("SLICE generic works as expected with a range", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    let start = 5;
    let end = 10;
    let sliced = bioc.SLICE(ggr, { start, end });
    expect(sliced.numberOfGroups()).toEqual(end - start);

    let sliced_lengths = [];
    for (var i = start; i < end; i++) {
        expect(sliced.group(i - start).start()).toEqual(grl_raw[i].start());
        sliced_lengths.push(bioc.LENGTH(grl_raw[i]));
    }

    expect(sliced.rangeStarts()[0]).toBe(0);
    expect(Array.from(sliced.rangeLengths())).toEqual(sliced_lengths);
})

test("COMBINE generic works as expected", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    let grl_raw2 = spawn_GRanges_list(15);
    let ggr2 = new bioc.GroupedGRanges(grl_raw2);

    let combined = bioc.COMBINE([ggr, ggr2]);
    expect(combined.ranges().seqnames()).toEqual(bioc.COMBINE([ggr.ranges().seqnames(), ggr2.ranges().seqnames()]));
    expect(combined.ranges().start()).toEqual(bioc.COMBINE([ggr.ranges().start(), ggr2.ranges().start()]));
    expect(combined.ranges().end()).toEqual(bioc.COMBINE([ggr.ranges().end(), ggr2.ranges().end()]));

    expect(combined.numberOfGroups()).toEqual(25);
    expect(combined.group(0).start()).toEqual(ggr.group(0).start());
    expect(combined.group(10).start()).toEqual(ggr2.group(0).start());
    expect(combined.group(24).start()).toEqual(ggr2.group(14).start());
})

test("CLONE generic works as expected", () => {
    let grl_raw = spawn_GRanges_list(10);
    let ggr = new bioc.GroupedGRanges(grl_raw);

    // Deep copy.
    let cloned = bioc.CLONE(ggr);
    cloned.ranges().start()[0] = 1000;
    expect(cloned.ranges().start()[0]).toBe(1000);
    expect(ggr.ranges().start()[0]).not.toBe(1000);

    // Shallow copy.
    let shallow = bioc.CLONE(ggr, { deepCopy: false });
    shallow.ranges().start()[0] = 1000;
    expect(shallow.ranges().start()[0]).toBe(1000);
    expect(ggr.ranges().start()[0]).toBe(1000);
})

test("overlaps work as expected", () => {
    let grl_raw = new bioc.GRanges(["A", "A", "B", "B"], new bioc.IRanges([1, 10, 20, 50], [5, 2, 7, 60]));
    let ggr = new bioc.GroupedGRanges(grl_raw, { rangeLengths: [2, 2] });
    let index = ggr.buildOverlapIndex();

    // Overlap with a GRanges.
    let query_gr = new bioc.GRanges(["A", "A", "A", "B", "B", "B"], new bioc.IRanges([1, 6, 11, 10, 30, 45], [3, 4, 5, 20, 10, 20]));
    {
        let olap = index.overlap(query_gr);
        expect(olap[0]).toEqual([0]);
        expect(olap[1]).toEqual([]);
        expect(olap[2]).toEqual([0]);

        expect(olap[3]).toEqual([1]);
        expect(olap[4]).toEqual([]);
        expect(olap[5]).toEqual([1]);

        // Redundant hits are uniquified.
        let full = new bioc.GRanges(["A", "B"], new bioc.IRanges([1, 1], [1000, 1000]));
        let olap2 = index.overlap(full);
        expect(olap2[0]).toEqual([0]);
        expect(olap2[1]).toEqual([1]);
    }

    // Overlap with a GRangesList.
    {
        let query_grl1 = new bioc.GroupedGRanges(query_gr, { rangeLengths: [ 3, 3 ] });
        let olap1 = index.overlap(query_grl1);
        expect(olap1[0]).toEqual([0]);
        expect(olap1[1]).toEqual([1]);

        let query_grl2 = new bioc.GroupedGRanges(query_gr, { rangeLengths: [ 4, 2 ] });
        let olap2 = index.overlap(query_grl2);
        expect(olap2[0]).toEqual([0, 1]);
        expect(olap2[1]).toEqual([1]);

        let query_grl3 = new bioc.GroupedGRanges(query_gr, { rangeLengths: [ 1, 1, 2, 1, 1 ] });
        let olap3 = index.overlap(query_grl3);
        expect(olap3[0]).toEqual([0]);
        expect(olap3[1]).toEqual([]);
        expect(olap3[2]).toEqual([0, 1]);
        expect(olap3[3]).toEqual([]);
        expect(olap3[4]).toEqual([1]);
    }
})
