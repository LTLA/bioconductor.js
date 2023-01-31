import * as bioc from "../src/index.js";

function spawnObject() {
    return { 
        "start": [1,3,5,7],
        "width": [10, 12, 5, 6]
    };
}

test("constructing a IRanges works (simple)", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width);

    expect(x.start()).toEqual(new Int32Array(obj.start));
    expect(x.width()).toEqual(new Int32Array(obj.width));
    expect(x.end()).toEqual(new Int32Array([ 11, 15, 10, 13 ]));

    expect(x.names()).toBeNull();
    expect(x.elementMetadata().numberOfRows()).toBe(obj.start.length);
    expect(x.elementMetadata().numberOfColumns()).toBe(0);

    expect(x.metadata()).toEqual(new Map);

    // Fails if the start and widths don't match up.
    expect(() => new bioc.IRanges(obj.start.slice(0, 2), obj.width)).toThrow("same length");

    // Fails if the start or widths are negative.
    expect(() => new bioc.IRanges([-1], [0])).toThrow("detected a negative entry");
    expect(() => new bioc.IRanges([0], [-1])).toThrow("detected a negative entry");
})

test("constructing a IRanges works with names", () => {
    let obj = spawnObject();
    let names = ["A", "B", "C", "D"];
    let x = new bioc.IRanges(obj.start, obj.width, { names: names });

    expect(x.names()).toEqual(names);
    expect(x.elementMetadata().rowNames()).toBeNull();

    // Fails if the names are of different length.
    expect(() => new bioc.IRanges(obj.start, obj.width, { names: ["A"] })).toThrow("should have length equal to");

    // Fails if the names are not strings.
    expect(() => new bioc.IRanges(obj.start, obj.width, { names: [0] })).toThrow("strings");
})

test("constructing a IRanges works with extra elementMetadata", () => {
    let obj = spawnObject();
    let df = new bioc.DataFrame({ thing: ["A", "B", "C", "D"] });
    let x = new bioc.IRanges(obj.start, obj.width, { elementMetadata: df });

    expect(x.elementMetadata().columnNames()).toEqual([ "thing" ]);
    expect(x.elementMetadata().column("thing")).toEqual([ "A", "B", "C", "D" ]);

    // Fails if the elementMetadata is not a DataFrame.
    expect(() => new bioc.IRanges(obj.start, obj.width, { elementMetadata: 1 })).toThrow("DataFrame");

    // Fails if the elementMetadata is not of the right length.
    expect(() => new bioc.IRanges(obj.start, obj.width, { elementMetadata: new bioc.DataFrame({ thing: [ "A" ] }) })).toThrow("number of rows");
})

test("IRanges setters work for start positions", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width);

    let y = x.setStart([9,8,7,5]);
    expect(y.start()[0]).toBe(9);
    expect(x.start()[0]).toBe(1); // same as before.

    x.$setStart([ 4, 5, 6, 7 ]);
    expect(x.start()).toEqual(new Int32Array([4, 5, 6, 7]));
    expect(x.end()).toEqual(new Int32Array([ 14, 17, 11, 13 ]));

    // Fails if not of the right length.
    expect(() => x.$setStart([ 1, 2, 3 ])).toThrow("same length")

    // Fails on negative values.
    expect(() => x.$setStart([ -1, 1, 1, 1 ])).toThrow("negative");
})

test("IRanges setters work for widths", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width);

    let y = x.setWidth([100, 200, 300, 400]);
    expect(y.width()[0]).toBe(100);
    expect(x.width()[0]).toBe(10); // same as before.

    x.$setWidth([ 4, 5, 6, 7 ]);
    expect(x.width()).toEqual(new Int32Array([4, 5, 6, 7]));
    expect(x.end()).toEqual(new Int32Array([ 5, 8, 11, 14 ]));

    // Fails if not of the right length.
    expect(() => x.$setWidth([ 1, 2, 3 ])).toThrow("same length")

    // Fails on negative values.
    expect(() => x.$setWidth([ -1, 1, 1, 1 ])).toThrow("negative");
})

test("IRanges setters work for names", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width);

    x.$setNames(["A", "B", "C", "D"]);
    expect(x.names()).toEqual(["A", "B", "C", "D"]);
    expect(x.elementMetadata().rowNames()).toBeNull();

    x.$setNames(null);
    expect(x.names()).toBeNull();

    // Fails if not of the right length.
    expect(() => x.$setNames(["A", "B"])).toThrow("replacement 'names'");

    // Fails if not strings.
    expect(() => x.$setNames([1,2,3,4])).toThrow("strings");
})

test("IRanges setters work for the elementMetadata", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width);
    let df = new bioc.DataFrame({ thing: ["A", "B", "C", "D"], foo: new Float64Array([1,2,3,4]) });

    let y = x.setElementMetadata(df);
    expect(x.elementMetadata().numberOfColumns()).toEqual(0);
    expect(y.elementMetadata().numberOfColumns()).toEqual(2);

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

test("IRanges setters work for the metadata", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width, { metadata: { thing: 2 } });
    expect(x.metadata().get("thing")).toEqual(2);

    let y = x.setMetadata({ boo: 5 });
    expect(y.metadata().get("boo")).toEqual(5);
    expect(x.metadata().has("boo")).toBe(false); // doesn't mutate the original.

    x.$setMetadata({ boo: 5 });
    expect(x.metadata().get("boo")).toEqual(5);
})

test("LENGTH generic works correctly for IRanges", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width);
    expect(bioc.LENGTH(x)).toEqual(obj.start.length);
})

test("SLICE generic works correctly for IRanges", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width, { elementMetadata: new bioc.DataFrame({ thing: [1,2,3,4] }) });

    let y = bioc.SLICE(x, [0, 3]);
    expect(y.start()).toEqual(new Int32Array([1, 7]));
    expect(y.width()).toEqual(new Int32Array([10, 6]));
    expect(y.elementMetadata().column("thing")).toEqual([1, 4]);

    // Works correctly for the names.
    x.$setNames(["A", "B", "C", "D"]);
    y = bioc.SLICE(x, [2, 1]);
    expect(y.names()).toEqual(["C", "B"]);
})

test("COMBINE generic works correctly for IRanges", () => {
    let obj = spawnObject();
    let x1 = new bioc.IRanges(obj.start, obj.width, { elementMetadata: new bioc.DataFrame({ thing: [1,2,3,4] }) });
    let x2 = new bioc.IRanges([100, 200], [50, 70], { elementMetadata: new bioc.DataFrame({ thing: [100, 200] }) });

    let combined = bioc.COMBINE([x1, x2]);
    expect(combined.start()).toEqual(new Int32Array([1,3,5,7,100,200]));
    expect(combined.width()).toEqual(new Int32Array([10,12,5,6,50,70]));
    expect(combined.elementMetadata().column("thing")).toEqual([1,2,3,4,100,200]);

    // Works with partial names.
    x1.$setNames(["A", "B", "C", "D"]);
    let combined2 = bioc.COMBINE([x1, x2]);
    expect(combined2.names()).toEqual(["A", "B", "C", "D", "", ""]);
})

test("CLONE generic works correctly for IRanges", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width, { elementMetadata: new bioc.DataFrame({ thing: [1,2,3,4] }) });
    let y = bioc.CLONE(x);

    expect(x.start()).toEqual(y.start());
    expect(x.width()).toEqual(y.width());
    expect(x.elementMetadata().column("thing")).toEqual(y.elementMetadata().column("thing"));

    // Modifying one object doesn't affect the other.
    x.$setStart([9,8,7,6]);
    expect(x.start()).toEqual(new Int32Array([9,8,7,6]));
    expect(y.start()).toEqual(new Int32Array([1,3,5,7]));

    // Same for the metadata.
    x.metadata().set("foo", 5);
    expect(x.metadata().has("foo")).toBe(true);
    expect(y.metadata().has("foo")).toBe(false);
})

test("overlap method works correctly for IRanges", () => {
    let x = new bioc.IRanges([1, 10, 20, 50], [5, 2, 7, 60]);
    let y = new bioc.IRanges([2, 30, 21], [8, 10, 30]);

    let idx = x.buildOverlapIndex();
    let results = idx.overlap(y);
    expect(results).toEqual([[0], [], [2, 3]]);
})
