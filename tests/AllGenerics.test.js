import * as bioc from "../src/index.js";

test("cloning works on the built-in types", () => {
    expect(bioc.CLONE(1)).toBe(1);
    expect(bioc.CLONE("asdasd")).toBe("asdasd");
    expect(bioc.CLONE(null)).toBeNull();

    {
        let x = [ "Y", "Z" ];
        let xclone = bioc.CLONE(x);
        expect(xclone).toEqual(x);
        xclone[0] = "A";
        expect(x[0]).toEqual("Y");
    }

    // Works for nested arrays.
    {
        let x = [ [ "A", "B" ], "Z" ];
        let xclone = bioc.CLONE(x);
        expect(xclone).toEqual(x);
        xclone[0][1] = "Y";
        expect(x[0][1]).toEqual("B");
    }

    // Works for objects.
    {
        let x = { a: "foo", b: "bar", c: [ "princess" ] };
        let xclone = bioc.CLONE(x);
        expect(xclone).toEqual(x);
        xclone.c[0] = "tatiana";
        expect(x.c[0]).toEqual("princess");
    }

    // Works for Maps.
    {
        let x = new Map;
        x.set("A", "Aaron");
        x.set("B", "Brian");

        let xclone = bioc.CLONE(x);
        expect(xclone.get("A")).toEqual("Aaron");
        expect(xclone.get("B")).toEqual("Brian");

        xclone.set("C", "Christine");
        expect(x.has("C")).toBe(false);
    }

    // Works for Sets.
    {
        let x = new Set([ "Aaron", "Brian" ]);

        let xclone = bioc.CLONE(x);
        expect(xclone.has("Aaron")).toBe(true);
        expect(xclone.has("Brian")).toBe(true);

        xclone.add("Christine");
        expect(x.has("Christine")).toBe(false);
    }
})

test("splitting works correctly for regular vectors", () => {
    let bar = [ "chino", "cocoa", "rize", "syaro", "chiya", "tippy" ]
    let school = [ "middle", "public", "private", "private", "public", "none" ]; 

    let fragmented = bioc.SPLIT(bar, school);
    expect(fragmented.middle).toEqual([ "chino" ]);
    expect(fragmented.public).toEqual([ "cocoa", "chiya" ]);
    expect(fragmented.private).toEqual([ "rize", "syaro" ]);
    expect(fragmented.none).toEqual([ "tippy" ]);

    // Respects type and presplits.
    let foo = new Int32Array([ 9, 8, 7, 6, 5, 4 ]);
    let presplit = bioc.presplitFactor(school);
    let fragmented2 = bioc.SPLIT(foo, presplit);

    expect(fragmented2.middle).toEqual(new Int32Array([ 9 ]));
    expect(fragmented2.public).toEqual(new Int32Array([ 8, 5 ]));
    expect(fragmented2.private).toEqual(new Int32Array([ 7, 6 ]));
    expect(fragmented2.none).toEqual(new Int32Array([ 4 ]));
})

