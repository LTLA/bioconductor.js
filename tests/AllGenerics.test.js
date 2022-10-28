import * as bioc from "../src/index.js";

test("cloning works on the built-in types", () => {
    expect(bioc.CLONE(1)).toBe(1);
    expect(bioc.CLONE("asdasd")).toBe("asdasd");

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
})

test("splitting works correctly", () => {
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

    let presplit = bioc.presplitFactor(school);
    let fragmented2 = bioc.SPLIT(df, presplit);

    expect(fragmented2.middle.column("foo")).toEqual([ 9 ]);
    expect(fragmented2.public.column("foo")).toEqual([ 8, 5 ]);
    expect(fragmented2.private.column("foo")).toEqual([ 7, 6 ]);
    expect(fragmented2.none.column("foo")).toEqual([ 4 ]);
})

