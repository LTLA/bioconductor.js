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

    expect(ggr.ranges().strand()[0]).toEqual(0);
    ggr.$setRanges(rr);
    expect(ggr.ranges().strand()[0]).toEqual(-1);

    // Errors
    expect(() => ggr.$setRanges(1)).toThrow("GRanges");
    expect(() => ggr.$setRanges(bioc.SLICE(rr, [0,1,2,3]))).toThrow("length");
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

