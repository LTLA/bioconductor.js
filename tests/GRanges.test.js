import * as bioc from "../src/index.js";

function spawnObject() {
    return { 
        "ranges": new bioc.IRanges([2, 4, 6, 8], [ 10, 20, 30, 40 ]),
        "seqnames": [ "chrA", "chrA", "chrB", "chrB" ]
    };
}

test("constructing a GRanges works (simple)", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);

    expect(x.start()).toEqual(obj.ranges.start());
    expect(x.width()).toEqual(obj.ranges.width());

    let str = x.strand();
    expect(str instanceof Int8Array).toBe(true);
    expect(str.every(x => x == 0)).toBe(true);

    expect(x.names()).toBeNull();
    expect(x.elementMetadata().numberOfRows()).toBe(obj.seqnames.length);
    expect(x.elementMetadata().numberOfColumns()).toBe(0);
    expect(x.metadata()).toEqual({});

    // Fails if the start and ranges don't match up.
    expect(() => new bioc.GRanges(obj.seqnames.slice(0, 2), obj.ranges)).toThrow("length equal to");
})

test("constructing a GRanges works with names", () => {
    let obj = spawnObject();
    let names = ["A", "B", "C", "D"];
    let x = new bioc.GRanges(obj.seqnames, obj.ranges, { names: names });

    expect(x.names()).toEqual(names);

    // Fails if the names are of different length.
    expect(() => new bioc.GRanges(obj.seqnames, obj.ranges, { names: ["A"] })).toThrow("should have length equal to");

    // Fails if the names are not strings.
    expect(() => new bioc.GRanges(obj.seqnames, obj.ranges, { names: [0] })).toThrow("strings");
})

test("constructing a GRanges works with extra elementMetadata", () => {
    let obj = spawnObject();
    let df = new bioc.DataFrame({ thing: ["A", "B", "C", "D"] });
    let x = new bioc.GRanges(obj.seqnames, obj.ranges, { elementMetadata: df });

    expect(x.elementMetadata().columnNames()).toEqual([ "thing" ]);
    expect(x.elementMetadata().column("thing")).toEqual([ "A", "B", "C", "D" ]);

    // Fails if the elementMetadata is not a DataFrame.
    expect(() => new bioc.GRanges(obj.seqnames, obj.ranges, { elementMetadata: 1 })).toThrow("DataFrame");

    // Fails if the elementMetadata is not of the right length.
    expect(() => new bioc.GRanges(obj.seqnames, obj.ranges, { elementMetadata: new bioc.DataFrame({ thing: [ "A" ] }) })).toThrow("number of rows");
})

test("constructing a GRanges works with non-default strands", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges, { strand: [ 1, 0, 1, -1 ] });

    let str = x.strand();
    expect(str).toEqual(new Int8Array([1, 0, 1, -1]));

    // Fails if not of the right length.
    expect(() => new bioc.GRanges(obj.seqnames, obj.ranges, { strand: [ 1 ] })).toThrow("should have length");

    // Fails for non-allowed values.
    expect(() => new bioc.GRanges(obj.seqnames, obj.ranges, { strand: [ 2, 3, 4, 5 ] })).toThrow("'strand' must be");
})

test("GRanges setters work for start positions", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);

    x.$setStart([ 4, 5, 6, 7 ]);
    expect(x.start()).toEqual(new Int32Array([4, 5, 6, 7]));
    expect(x.end()).toEqual(new Int32Array([ 14, 25, 36, 47 ]));
})

test("GRanges setters work for widths", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);

    x.$setWidth([ 4, 5, 6, 7 ]);
    expect(x.width()).toEqual(new Int32Array([4, 5, 6, 7]));
    expect(x.end()).toEqual(new Int32Array([ 6, 9, 12, 15 ]));
})

test("GRanges setters work for seqnames", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);

    x.$setSeqnames([ "chrX", "chrX", "chrY", "chrY" ]);
    expect(x.seqnames()).toEqual([ "chrX", "chrX", "chrY", "chrY" ]);

    // Fails if not of the right length.
    expect(() => x.$setSeqnames([ "A", "B" ])).toThrow("should have length")

    // Fails if they're not strings.
    expect(() => x.$setSeqnames([ 0, 1, 2, 3 ])).toThrow("strings")
})

test("GRanges setters work for names", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);

    x.$setNames(["A", "B", "C", "D"]);
    expect(x.names()).toEqual(["A", "B", "C", "D"]);

    x.$setNames(null);
    expect(x.names()).toBeNull();

    // Fails if not of the right length.
    expect(() => x.$setNames(["A", "B"])).toThrow("replacement 'names'");

    // Fails if not strings.
    expect(() => x.$setNames([1,2,3,4])).toThrow("strings");
})

test("GRanges setters work for the ranges", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);

    let replacement = new bioc.IRanges([9,7,5,3], [100, 200, 300, 400]);
    x.$setRanges(replacement);
    expect(x.start()).toEqual(replacement.start());
    expect(x.width()).toEqual(replacement.width());

    // Fails if not an IRanges.
    expect(() => x.$setRanges(1)).toThrow("IRanges");

    // Fails if not the right length.
    expect(() => x.$setRanges(bioc.SLICE(replacement, [0, 1]))).toThrow("should have length");
})

test("GRanges setters work for the strands", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);

    x.$setStrand([1, 0, 1, -1]);
    expect(x.strand()).toEqual(new Int8Array([1, 0, 1, -1]));

    // Fails if not the right length.
    expect(() => x.$setStrand([1])).toThrow("should have length");

    // Fails with illegal values.
    expect(() => x.$setStrand([2,3,4,5])).toThrow("'strand' must be");
})

test("GRanges setters work for the elementMetadata", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);
    let df = new bioc.DataFrame({ thing: ["A", "B", "C", "D"], foo: new Float64Array([1,2,3,4]) });

    x.$setElementMetadata(df);
    expect(x.elementMetadata().numberOfRows()).toEqual(4);
    expect(x.elementMetadata().numberOfColumns()).toEqual(2);

    // Setting to null wipes the information.
    x.$setElementMetadata(null);
    expect(x.elementMetadata().numberOfRows()).toEqual(4);
    expect(x.elementMetadata().numberOfColumns()).toEqual(0);

    // Fails if it's not a DataFrame.
    expect(() => x.$setElementMetadata(1)).toThrow("DataFrame")

    // Fails if it's not of the right rows.
    expect(() => x.$setElementMetadata(new bioc.DataFrame({ thing: [1] }))).toThrow("number of rows")
})

test("GRanges setters work for the metadata", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges, { metadata: { thing: 2 } });
    expect(x.metadata()).toEqual({ thing: 2 });

    x.$setMetadata({ boo: 5 });
    expect(x.metadata()).toEqual({ boo: 5 });
})

test("LENGTH generic works correctly for GRanges", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges);
    expect(bioc.LENGTH(x)).toEqual(obj.seqnames.length);
})

test("SLICE generic works correctly for GRanges", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges, { elementMetadata: new bioc.DataFrame({ thing: [1,2,3,4] }) });

    let y = bioc.SLICE(x, [0, 3]);
    expect(y.start()).toEqual(new Int32Array([2, 8]));
    expect(y.width()).toEqual(new Int32Array([10, 40]));
    expect(y.elementMetadata().column("thing")).toEqual([1, 4]);

    // Works correctly for the names.
    x.$setNames(["A", "B", "C", "D"]);
    y = bioc.SLICE(x, [2, 1]);
    expect(y.names()).toEqual(["C", "B"]);
})

test("COMBINE generic works correctly for GRanges", () => {
    let obj = spawnObject();
    let x1 = new bioc.GRanges(obj.seqnames, obj.ranges, { elementMetadata: new bioc.DataFrame({ thing: [1,2,3,4] }) });
    let x2 = new bioc.GRanges(["chrX", "chrY"], new bioc.IRanges([100, 200], [50, 70]), { elementMetadata: new bioc.DataFrame({ thing: [100, 200] }) });

    let combined = bioc.COMBINE([x1, x2]);
    expect(combined.start()).toEqual(new Int32Array([2,4,6,8,100,200]));
    expect(combined.width()).toEqual(new Int32Array([10,20,30,40,50,70]));
    expect(combined.elementMetadata().column("thing")).toEqual([1,2,3,4,100,200]);

    // Works with partial names.
    x1.$setNames(["A", "B", "C", "D"]);
    let combined2 = bioc.COMBINE([x1, x2]);
    expect(combined2.names()).toEqual(["A", "B", "C", "D", "", ""]);
})

test("CLONE generic works correctly for GRanges", () => {
    let obj = spawnObject();
    let x = new bioc.GRanges(obj.seqnames, obj.ranges, { elementMetadata: new bioc.DataFrame({ thing: [1,2,3,4] }) });
    let y = bioc.CLONE(x);

    expect(x.start()).toEqual(y.start());
    expect(x.width()).toEqual(y.width());
    expect(x.elementMetadata().column("thing")).toEqual(y.elementMetadata().column("thing"));

    // Modifying one object doesn't affect the other.
    x.$setStart([9,8,7,6]);
    expect(x.start()).toEqual(new Int32Array([9,8,7,6]));
    expect(y.start()).toEqual(new Int32Array([2,4,6,8]));

    // Same for the metadata.
    x.metadata().foo = 5;
    expect("foo" in x.metadata()).toBe(true);
    expect("foo" in y.metadata()).toBe(false);
})

test("overlap method works correctly for GRanges", () => {
    // Ignoring the strand.
    let x = new bioc.GRanges(["A", "B", "A", "B"], new bioc.IRanges([1, 10, 20, 50], [5, 2, 7, 60]));
    let y = new bioc.GRanges(["A", "B", "A"], new bioc.IRanges([2, 0, 21], [8, 60, 30]));

    let idx = x.buildOverlapIndex();
    let results = idx.overlap(y);
    expect(results).toEqual([[0], [3, 1], [2]]);

    // Strand is effectively ignored when dealing with '*' references.
    y.$setStrand([ 1, -1, -1 ]);
    let results2 = idx.overlap(y);
    expect(results2).toEqual(results);

    // Now respecting the strand.
    x.$setStrand([ -1, -1, 1, 1 ]);
    let idx3 = x.buildOverlapIndex();
    let results3 = idx3.overlap(y, { ignoreStrand: false });
    expect(results3).toEqual([[], [1], []]);
})
