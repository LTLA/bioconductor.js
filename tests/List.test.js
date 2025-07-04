import * as bioc from "../src/index.js";

test("constructing a List from an array", () => {
    let ll = new bioc.List([1,2,3,4,5]);
    expect(ll.length()).toEqual(5);
    expect(ll.names()).toBeNull();
    expect(ll.values()).toEqual([1,2,3,4,5]);

    ll = new bioc.List([1,2,3,4,5], { names: ["a", "b", "c", "d", "e"] });
    expect(ll.names()).toEqual(['a', 'b', 'c', 'd', 'e']);

    expect(() => new bioc.List([1,2,3,4,5], { names: ["a", "b", "c", "d"] })).toThrow("same length");
    expect(() => new bioc.List([1,2,3,4,5], { names: ["a", 2, "b", "c", "d"] })).toThrow("array of strings");
})

test("constructing a List from a Map", () => {
    let src = new Map;
    src.set("A", 1);
    src.set("B", 2);
    src.set("C", 3);
    src.set("D", 4);
    src.set("E", 5);

    let ll = new bioc.List(src);
    expect(ll.length()).toEqual(5);
    expect(ll.names()).toEqual(["A","B","C","D","E"]);
    expect(ll.values()).toEqual([1,2,3,4,5]);

    let nonstrkey = new Map;
    nonstrkey.set(1, 2);
    expect(() => new bioc.List(nonstrkey)).toThrow("should be strings");

    ll = new bioc.List(src, { names: ["E", "D", "C", "B", "A"] });
    expect(ll.names()).toEqual(['E', 'D', 'C', 'B', 'A']);
    expect(ll.values()).toEqual([5,4,3,2,1]);

    expect(() => new bioc.List(src, { names: ["A", "B", "C", "D"] })).toThrow("equal to length");
    expect(() => new bioc.List(src, { names: ["A", 2, "B", "C", "D"] })).toThrow("array of strings");
    expect(() => new bioc.List(src, { names: ["A", "B", "C", "D", "F"] })).toThrow("missing name");
})

test("constructing a List from an object", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };

    let ll = new bioc.List(src);
    expect(ll.length()).toEqual(5);
    expect(ll.names()).toEqual(["A","B","C","D","E"]);
    expect(ll.values()).toEqual([1,2,3,4,5]);

    ll = new bioc.List(src, { names: ["E", "D", "C", "B", "A"] });
    expect(ll.names()).toEqual(['E', 'D', 'C', 'B', 'A']);
    expect(ll.values()).toEqual([5,4,3,2,1]);

    expect(() => new bioc.List(src, { names: ["A", "B", "C", "D"] })).toThrow("equal to length");
    expect(() => new bioc.List(src, { names: ["A", 2, "B", "C", "D"] })).toThrow("array of strings");
    expect(() => new bioc.List(src, { names: ["A", "B", "C", "D", "F"] })).toThrow("missing name");
})

test("List getters work as expected", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);

    expect(ll.get(0)).toEqual(1);
    expect(ll.get(4)).toEqual(5);

    expect(ll.get("A")).toEqual(1);
    expect(ll.get("E")).toEqual(5);
    expect(ll.get("C")).toEqual(3);

    expect(ll.nameToIndex("B")).toEqual(1);
    expect(ll.nameToIndex("D")).toEqual(3);

    expect(() => ll.get(5)).toThrow("out of range");
    expect(() => ll.get(-1)).toThrow("out of range");
    expect(() => ll.get("F")).toThrow("no matching name");

    let unnamed = new bioc.List([6,7,8]);
    expect(() => unnamed.get("F")).toThrow("no available names");
})

test("List setters work with indices", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);
    let unnamed = new bioc.List([1,2,3,4,5]);

    {
        let ll2 = ll.set(1, null);
        expect(ll2.get(1)).toBeNull();
        expect(ll.get(1)).toEqual(2);
    }

    {
        let ll2 = ll.set(1, null, { name: "foo" });
        expect(ll2.get(1)).toBeNull();
        expect(ll2.names()).toEqual(["A","foo","C","D","E"]);
        expect(ll.names()).toEqual(["A","B","C","D","E"]);

        let unnamed2 = unnamed.set(1, null, { name: "foo" });
        expect(unnamed2.get(1)).toBeNull();
        expect(unnamed2.names()).toEqual(["","foo","","",""]);
        expect(unnamed.names()).toBeNull();

        // Check that the lookup isn't contaminated in the existing object. 
        expect(() => ll.get("foo")).toThrow("no matching name");
        expect(ll2.get("foo")).toBeNull();
    }

    // Appends automatically.
    {
        let ll2 = ll.set(5, null);
        expect(ll2.length()).toEqual(6);
        expect(ll2.get(5)).toEqual(null);
        expect(ll2.names()).toEqual(["A","B","C","D","E",""]);
        expect(ll.length()).toEqual(5);

        expect(() => ll.set(6, null)).toThrow("out of range");
    }

    {
        let ll2 = ll.set(5, null, { name: "foo" });
        expect(ll2.length()).toEqual(6);
        expect(ll2.get(5)).toEqual(null);
        expect(ll2.names()).toEqual(["A","B","C","D","E","foo"]);
        expect(ll.length()).toEqual(5);

        let unnamed2 = unnamed.set(5, null, { name: "foo" });
        expect(unnamed2.get(5)).toBeNull();
        expect(unnamed2.names()).toEqual(["","","","","","foo"]);
        expect(unnamed.names()).toBeNull();

        expect(() => unnamed.set(5, null, { name: 2 })).toThrow("should be a string");

        // Check that the lookup isn't contaminated in the existing object. 
        expect(() => ll.get("foo")).toThrow("no matching name");
        expect(ll2.get("foo")).toBeNull();
    }

    // Works in-place as well.
    ll.set(5, null, { inPlace: true });
    expect(ll.length()).toEqual(6);
    expect(ll.get(5)).toBeNull();
    expect(ll.get("")).toBeNull();

    ll.set(5, null, { inPlace: true, name: "foo" });
    expect(ll.length()).toEqual(6);
    expect(() => ll.get("")).toThrow("no matching name"); // check that the lookup is properly reset.
    expect(ll.get("foo")).toBeNull();
})

test("List setters work with names", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);
    let unnamed = new bioc.List([1,2,3,4,5]);

    {
        let ll2 = ll.set("A", null);
        expect(ll2.get(0)).toBeNull();
        expect(ll2.get("A")).toBeNull();
    }

    {
        let ll2 = ll.set("F", null);
        expect(ll2.names()).toEqual(["A","B","C","D","E","F"]);
        expect(ll2.get(5)).toBeNull();
        expect(ll2.get("F")).toBeNull();
    }

    {
        let unnamed2 = unnamed.set("foo", null);
        expect(unnamed2.names()).toEqual(["","","","","","foo"]);
        expect(unnamed2.get(5)).toBeNull();
        expect(unnamed2.get("foo")).toBeNull();

        // Check that the lookup isn't contaminated in the existing object. 
        expect(() => unnamed.get("foo")).toThrow("no available names");
        expect(unnamed2.get("foo")).toBeNull();
    }

    // Works in-place as well.
    ll.set("A", null, { inPlace: true });
    expect(ll.get(0)).toBeNull();
})

test("name replacement works in a List", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);

    let unnamed = ll.setNames(null);
    expect(unnamed.names()).toBeNull();
    expect(ll.names()).not.toBeNull();

    let renamed = ll.setNames(['e','d','c','b','a']);
    expect(renamed.names()).toEqual(['e','d','c','b','a']);
    expect(ll.names()).toEqual(['A','B','C','D','E']);

    expect(() => ll.setNames(['A'])).toThrow("same length");
    expect(() => ll.setNames(['A', 1, 2, 3, 4])).toThrow("array of strings");

    ll.setNames(['e','d','c','b','a'], { inPlace: true });
    expect(ll.names()).toEqual(['e','d','c','b','a']);
    expect(() => ll.get("A")).toThrow("no matching name");
    expect(ll.get("a")).toEqual(5);
})
