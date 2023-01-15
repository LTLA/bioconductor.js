import * as bioc from "../src/index.js";

function spawnObject() {
    return  { 
        "A": new Int32Array([ 1, 2, 3, 4 ]),
        "B": [ 'x', 'y', 'z', 'aa' ],
    };
}

test("constructing a DataFrame works (simple)", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);

    expect(x.columnNames().slice().sort()).toEqual(["A", "B"]);
    expect(x.numberOfRows()).toEqual(4);
    expect(x.numberOfColumns()).toEqual(2);
    expect(x.rowNames()).toBeNull();

    expect(x.column("A").constructor.name).toEqual("Int32Array");
    expect(Array.from(x.column("A"))).toEqual([1,2,3,4]);;
    expect(x.column(1)).toEqual(["x", "y", "z", "aa"]);

    // Still works if the number of rows is set.
    x = new bioc.DataFrame(obj, { numberOfRows: 4 });
    expect(x.numberOfRows()).toEqual(4);

    expect(() => new bioc.DataFrame(obj, { numberOfRows: 5 })).toThrow("expected all arrays in 'columns' to have equal length");
    obj.B.push(5);
    expect(() => new bioc.DataFrame(obj)).toThrow("expected all arrays in 'columns' to have equal length");
    obj.B.pop();
})

test("constructing a DataFrame works (custom column order)", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj, { columnOrder: [ "B", "A" ] });
    expect(x.columnNames()).toEqual(["B", "A"]);

    expect(() => new bioc.DataFrame(obj, { columnOrder: [ "B" ] })).toThrow("'columnOrder' should have the same length");
    expect(() => new bioc.DataFrame(obj, { columnOrder: [ "B", "C" ] })).toThrow("values of 'columnOrder' should be the same");
})

test("constructing a DataFrame works (with rownames)", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj, { rowNames: [ "alpha", "bravo", "charlie", "delta"] }); 

    expect(x.numberOfRows()).toEqual(4);
    expect(x.rowNames()).toEqual(["alpha", "bravo", "charlie", "delta"]);

    expect(() => new bioc.DataFrame(obj, { rowNames: [ "B" ] })).toThrow("'rowNames' array");
})

test("constructing a DataFrame works (with empty objects)", () => {
    let x = new bioc.DataFrame({}, { rowNames: [ "a", "b", "c", "d" ] });
    expect(x.numberOfRows()).toEqual(4);
    expect(x.rowNames().length).toEqual(4);

    x = new bioc.DataFrame({}, { numberOfRows: 4 });
    expect(x.numberOfRows()).toEqual(4);
    expect(x.rowNames()).toBeNull();

    x = new bioc.DataFrame({});
    expect(x.numberOfRows()).toEqual(0);
    expect(x.rowNames()).toBeNull();
})

test("setting a column works correctly", () => {
    let x = new bioc.DataFrame({}, { numberOfRows: 4 });
    x.$setColumn("A", [1,2,3,4])
    expect(x.column("A")).toEqual([1,2,3,4]);
    expect(x.columnNames()).toEqual(["A"]);

    x.$setColumn("B", new Uint8Array([1,0,1,1]));
    expect(Array.from(x.column("B"))).toEqual([1, 0, 1, 1]);
    expect(x.columnNames()).toEqual(["A", "B"]);

    x.$setColumn(0, ["aa", "b", "cc", "d"]);
    expect(x.column("A")).toEqual(["aa", "b", "cc", "d"]);

    x.$setColumn("B", [5,6,7,8]);
    expect(x.column("B")).toEqual([5,6,7,8]);

    expect(() => x.$setColumn("C", [1])).toThrow("same length");
    expect(() => x.$setColumn(2, [1,2,3,4])).toThrow("out of range");
});

test("removing a column works correctly", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);
    x.$removeColumn("A");
    expect(x.columnNames()).toEqual(["B"]);
    expect(() => x.column("A")).toThrow("A");

    obj = spawnObject();
    x = new bioc.DataFrame(obj);
    x.$removeColumn(1);
    expect(x.columnNames()).toEqual(["A"]);
    expect(() => x.column("B")).toThrow("B");
})

test("setting column names works correctly", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);
    x.$setColumnNames(["alpha", "bravo"]);

    expect(x.column("alpha")).toEqual(obj.A);
    expect(x.column("bravo")).toEqual(obj.B);
    expect(x.columnNames()).toEqual(["alpha", "bravo"]);
})

test("slicing columns works correctly", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);
    x.$sliceColumns(["B", "A"]);

    expect(x.column("A")).toEqual(obj.A);
    expect(x.column("B")).toEqual(obj.B);
    expect(x.columnNames()).toEqual(["B", "A"]);

    x.$sliceColumns([1]);
    expect(x.column("A")).toEqual(obj.A);
    expect(() => x.column("B")).toThrow("B");
    expect(x.columnNames()).toEqual(["A"]);

    expect(() => x.$sliceColumns([0, 0])).toThrow("duplicate");
})

test("cloning an array collection works", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);
    let y = bioc.CLONE(x);
    expect(y.numberOfRows()).toEqual(4);

    y.column("A")[0] = 2;
    y.$setColumn("C", [5,4,3,2]);

    expect(x.column("A")[0]).toBe(1);
    expect(x.hasColumn("C")).toBe(false);

    // Works with row names.
    x = new bioc.DataFrame(obj, { rowNames: [ "a", "b", "c", "d" ] });
    expect(x.rowNames().length).toEqual(4);

    y = bioc.CLONE(x);
    y.$setRowNames(null);

    expect(y.rowNames()).toBeNull();
    expect(x.rowNames().length).toEqual(4);
})

test("cloning an array collection works (shallow)", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj, { rowNames: [ "a", "b", "c", "d" ] });
    let y = bioc.CLONE(x, { deepCopy: false });

    y.column("A")[0] = 2;
    y.$setColumn("C", [5,4,3,2]);
    y.$setRowNames(null);

    expect(x.column("A")[0]).toBe(2); // inner values are still referenced.
    expect(x.hasColumn("C")).toBe(false);
    expect(x.rowNames().length).toEqual(4);
})

test("slicing a DataFrame works (indices)", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);
    let out = bioc.SLICE(x, [3, 1, 2]);

    expect(out.numberOfRows()).toBe(3);
    expect(x.column("A").constructor.name).toEqual("Int32Array");
    expect(Array.from(out.column("A"))).toEqual([4,2,3]);
    expect(out.column("B")).toEqual(['aa','y','z']);

    // Works with row names in the house.
    x = new bioc.DataFrame(obj, { rowNames: [ "a", "b", "c", "d" ] });
    out = bioc.SLICE(x, [3, 1, 2]);
    expect(out.rowNames()).toEqual(["d", "b", "c"]);
})

test("slicing a DataFrame works (range)", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);
    let out = bioc.SLICE(x, {start: 1, end: 3});

    expect(out.numberOfRows()).toBe(2);
    expect(Array.from(out.column("A"))).toEqual([2,3]);
})

test("slicing a DataFrame works (view)", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);
    let out = bioc.SLICE(x, {start: 1, end: 3}, { allowView: true });

    out.column("A")[0] = 100;
    expect(x.column("A")[1]).toBe(100);
})

test("combining multiple array collections works", () => {
    let obj = spawnObject();
    let stub = { "A": new Int32Array([ 5, 6 ]), "B":[ 'bb', 'cc' ]};

    // Simple case.
    {
        let x = new bioc.DataFrame(obj);
        let y = new bioc.DataFrame(stub);
        let out = bioc.COMBINE([x, y]);

        expect(out.numberOfRows()).toBe(6);
        expect(x.column("A").constructor.name).toEqual("Int32Array");
        expect(Array.from(out.column("A"))).toEqual([1,2,3,4,5,6]);
        expect(out.column("B")).toEqual(["x", "y", "z", "aa", "bb", "cc"]);
    }

    // Multiple objects.
    {
        let x = new bioc.DataFrame(obj);
        let y = new bioc.DataFrame(stub);
        let out = bioc.COMBINE([y, x, y]);

        expect(out.numberOfRows()).toBe(8);
        expect(x.column("A").constructor.name).toEqual("Int32Array");
        expect(Array.from(out.column("A"))).toEqual([5,6,1,2,3,4,5,6]);
        expect(out.column("B")).toEqual(["bb", "cc", "x", "y", "z", "aa", "bb", "cc"]);
    }

    // Works for empty objects.
    {
        let x = new bioc.DataFrame({}, { numberOfRows: 5 });
        let y = new bioc.DataFrame({}, { numberOfRows: 4 });
        let out = bioc.COMBINE([x, y]);
        expect(out.numberOfRows()).toBe(9);
        expect(out.numberOfColumns()).toBe(0);
    }

    // Works for row names.
    {
        let x = new bioc.DataFrame(obj, { rowNames: [ "a", "b", "c", "d" ]});
        let y = new bioc.DataFrame(stub, { rowNames: [ "AA", "BB" ] });
        let out = bioc.COMBINE([x, y]);
        expect(out.rowNames()).toEqual(["a", "b", "c", "d", "AA", "BB"]);

        y.$setRowNames(null);
        out = bioc.COMBINE([x, y]);
        expect(out.rowNames()).toEqual(["a", "b", "c", "d", "", ""]);
    }
})

test("combining multiple array collections preserves TypedArray types", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj);
    let stub = { "A": new Float32Array([ 5, 6 ]), "B":[ 'bb', 'cc' ]};
    let y = new bioc.DataFrame(stub);

    let out = bioc.COMBINE([x, y]);
    expect(out.column("A").constructor.name).toBe("Float64Array"); // promoted.
    expect(Array.from(out.column("A"))).toEqual([1,2,3,4,5,6]);

    stub.A = [null, null];
    out = bioc.COMBINE([x, y]);
    expect(out.column("A").constructor.name).toBe("Array"); // promoted.
    expect(Array.from(out.column("A"))).toEqual([1,2,3,4,null,null]);

    // Handles BigInts with some grace.
    stub.A = new BigInt64Array([5n, 6n]);
    out = bioc.COMBINE([x, y]);
    expect(out.column("A").constructor.name).toBe("Array"); // promoted.
    expect(Array.from(out.column("A"))).toEqual([1,2,3,4,5n,6n]);

    obj.A = new BigInt64Array([1n,2n,3n,4n]);
    out = bioc.COMBINE([x, y]);
    expect(out.column("A").constructor.name).toBe("BigInt64Array"); // promoted.
    expect(Array.from(out.column("A"))).toEqual([1n,2n,3n,4n,5n,6n]);
})

test("fliexbly combining DataFrames by row", () => {
    let x = new bioc.DataFrame({ "A": [ 1, 2, 3, 4 ], "B":[ 'x', 'y', 'z', 'aa' ]});
    let y = new bioc.DataFrame({ "A": [ 5, 6 ], "B":[ 'bb', 'cc' ]});
    let z = new bioc.DataFrame({ "A": [5, 6] });

    let out = bioc.flexibleCombineRows([x, z, y]);
    expect(out.column("A")).toEqual([1,2,3,4,5,6,5,6]);
    expect(out.column("B")).toEqual(["x", "y", "z", "aa", null, null, "bb", "cc"]);

    z.$removeColumn("A");
    out = bioc.flexibleCombineRows([x, z, y]);
    expect(out.column("A")).toEqual([1,2,3,4,null,null,5,6]);
    expect(out.column("B")).toEqual(["x", "y", "z", "aa", null, null, "bb", "cc"]);
})

test("metadata operations are respected for DataFrames", () => {
    let obj = spawnObject();
    let x = new bioc.DataFrame(obj, { metadata: { "foo": 2, "bar": 3 } });
    expect(x.metadata().foo).toBe(2);
    expect(x.$setMetadata({ whee: 1 }).metadata().whee).toBe(1);

    // Picked up by the clone.
    let y = bioc.CLONE(x);
    y.metadata().ark = 2;
    expect(y.metadata().ark).toBe(2);
    expect(x.metadata().ark).toBeUndefined();

    // Preserved in the COMBINE operation.
    let z = bioc.COMBINE([x, y]);
    expect(z.numberOfRows()).toBe(x.numberOfRows() + y.numberOfRows());
    expect(z.metadata().whee).toBe(1);
})

test("splitting works correctly for DataFrames", () => {
    let df = new bioc.DataFrame({ 
        foo: [ 9, 8, 7, 6, 5, 4 ], 
        bar: [ "chino", "cocoa", "rize", "syaro", "chiya", "tippy" ]
    });

    let school = [ "middle", "public", "private", "private", "public", "none" ]; 
    let fragmented = bioc.SPLIT(df, school);

    expect(fragmented.middle.column("bar")).toEqual([ "chino" ]);
    expect(fragmented.public.column("bar")).toEqual([ "cocoa", "chiya" ]);
    expect(fragmented.private.column("bar")).toEqual([ "rize", "syaro" ]);
    expect(fragmented.none.column("bar")).toEqual([ "tippy" ]);
})

