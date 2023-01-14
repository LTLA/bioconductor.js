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

    expect(x.hasNames()).toBe(false);
    expect(x.rangeMetadata().numberOfRows()).toBe(obj.start.length);
    expect(x.rangeMetadata().numberOfColumns()).toBe(0);

    expect(x.metadata()).toEqual({});

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

    expect(x.hasNames()).toBe(true);
    expect(x.names()).toEqual(names);

    // Fails if the names are of different length.
    expect(() => new bioc.IRanges(obj.start, obj.width, { names: ["A"] })).toThrow("same length");

    // Fails if the names are not strings.
    expect(() => new bioc.IRanges(obj.start, obj.width, { names: [0] })).toThrow("strings");
})

test("constructing a IRanges works with extra rangeMetadata", () => {
    let obj = spawnObject();
    let df = new bioc.DataFrame({ thing: ["A", "B", "C", "D"] });
    let x = new bioc.IRanges(obj.start, obj.width, { rangeMetadata: df });

    expect(x.rangeMetadata().columnNames()).toEqual([ "thing" ]);
    expect(x.rangeMetadata().column("thing")).toEqual([ "A", "B", "C", "D" ]);

    // Fails if the rangeMetadata is not a DataFrame.
    expect(() => new bioc.IRanges(obj.start, obj.width, { rangeMetadata: 1 })).toThrow("DataFrame");

    // Fails if the rangeMetadata is not of the right length.
    expect(() => new bioc.IRanges(obj.start, obj.width, { rangeMetadata: new bioc.DataFrame({ thing: [ "A" ] }) })).toThrow("number of rows");
})

test("IRanges setters work for start positions", () => {
    let obj = spawnObject();
    let x = new bioc.IRanges(obj.start, obj.width);

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

    x.$setWidth([ 4, 5, 6, 7 ]);
    expect(x.width()).toEqual(new Int32Array([4, 5, 6, 7]));
    expect(x.end()).toEqual(new Int32Array([ 5, 8, 11, 14 ]));

    // Fails if not of the right length.
    expect(() => x.$setWidth([ 1, 2, 3 ])).toThrow("same length")

    // Fails on negative values.
    expect(() => x.$setWidth([ -1, 1, 1, 1 ])).toThrow("negative");
})

