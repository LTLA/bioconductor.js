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

    let dupped = new bioc.List([1,2,3,4,5], { names: ["a", "b", "a", "c", "c"] });
    expect(dupped.get("a")).toBe(1);
    expect(dupped.get("b")).toBe(2);
    expect(dupped.get("c")).toBe(4);
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

    // Works with duplicates.
    {
        let dupped = new bioc.List([1,2,3,4,5], { names: ["a", "b", "a", "c", "c"] });
        dupped.set(0, null, { name: "b", inPlace: true });
        expect(dupped.get("b")).toBeNull();
        expect(dupped.get("a")).toBe(3);
        dupped.set(0, null, { name: "a", inPlace: true });
        expect(dupped.get("a")).toBeNull();
        expect(dupped.get("b")).toBe(2);
        dupped.set(1, null, { name: "a", inPlace: true });
        expect(dupped.get("a")).toBeNull();
        expect(dupped.has("b")).toBe(false);
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

    // Works with duplicates.
    {
        let dupped = new bioc.List([1,2,3,4,5], { names: ["a", "b", "a", "c", "c"] });
        dupped.set("a", null, { inPlace: true });
        expect(dupped.get(0)).toBeNull();
        expect(dupped.get("a")).toBeNull();
        expect(dupped.get(2)).toBe(3);
        dupped.set("c", null, { inPlace: true });
        expect(dupped.get(3)).toBeNull();
        expect(dupped.get("c")).toBeNull();
        expect(dupped.get(4)).toBe(5);
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

test("List deleters work with indices", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);
    let unnamed = new bioc.List([1,2,3,4,5]);

    {
        let ll2 = ll.delete(2);
        expect(ll2.get(0)).toBe(1);
        expect(ll2.get(2)).toBe(4);
        expect(ll2.get("B")).toBe(2);
        expect(ll2.get("D")).toBe(4);
        expect(() => ll2.get("C")).toThrow("no matching name");
        expect(ll2.names()).toEqual(["A","B","D","E"]);
    }

    {
        let unnamed2 = unnamed.delete(1);
        expect(unnamed2.get(0)).toBe(1);
        expect(unnamed2.get(1)).toBe(3);
    }

    // Works with duplicates.
    {
        let dupped = new bioc.List([1,2,3,4,5], { names: ["a", "b", "a", "c", "c"] });
        let deleted = dupped.delete(0);
        expect(deleted.get(0)).toBe(2);
        expect(deleted.get("a")).toBe(3);
        deleted = dupped.delete(3);
        expect(deleted.get("c")).toBe(5);
    }

    // Works in place.
    expect(ll.get("D")).toEqual(4); // populating the lookup.
    ll.delete(3, { inPlace: true });
    expect(ll.names()).toEqual(["A","B","C","E"]);
    expect(() => ll.get("D")).toThrow("no matching name"); // lookup is properly reset.
})

test("List deleters work with names", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);

    // Test error code without a lookup.
    expect(() => ll.delete("F")).toThrow("no matching name");

    {
        let ll2 = ll.delete("C"); // checking deletion without a lookup.
        expect(ll2.get(0)).toBe(1);
        expect(ll2.get(2)).toBe(4);
        expect(ll2.get("B")).toBe(2);
        expect(ll2.get("D")).toBe(4);
        expect(() => ll2.get("C")).toThrow("no matching name");
        expect(ll2.names()).toEqual(["A","B","D","E"]);
    }

    {
        let res = ll.get("C"); // populate the look-up table.
        expect(res).toBe(3); // existing objects are unaffected.
        let ll2 = ll.delete("A"); // checking deletion with a lookup.
        expect(ll2.get(0)).toBe(2);
    }

    {
        let unnamed = new bioc.List([1,2,3,4,5]);
        expect(() => unnamed.delete("A")).toThrow("no available names");
    }

    // Works with duplicates.
    {
        let dupped = new bioc.List([1,2,3,4,5], { names: ["a", "b", "a", "c", "c"] });
        let deleted = dupped.delete("a");
        expect(deleted.get(0)).toBe(2);
        expect(deleted.get("a")).toBe(3);
        deleted = deleted.delete("c");
        expect(deleted.get("c")).toBe(5);
    }

    // Works in place.
    expect(ll.get("D")).toEqual(4); // populating the lookup.
    ll.delete("C", { inPlace: true });
    expect(ll.names()).toEqual(["A","B","D","E"]);
    expect(() => ll.get("C")).toThrow("no matching name"); // lookup is properly reset.
})

test("List conversion to different types", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);

    expect(ll.toObject()).toEqual(src);
    expect(ll.toArray()).toEqual([1,2,3,4,5]);
    let map = ll.toMap();
    expect(map.size).toBe(5);
    expect(map.get("A")).toBe(1);
    expect(map.get("E")).toBe(5);

    // Handles duplicates correctly.
    ll = new bioc.List([1,2,3,4,5], { names: ["a", "b", "a", "c", "b"] });
    expect(ll.toArray()).toEqual([1,2,3,4,5]);
    expect(ll.toObject()).toEqual({ a : 1, b: 2, c : 4 });
    map = ll.toMap();
    expect(map.size).toBe(3);
    expect(map.get("a")).toBe(1);
    expect(map.get("b")).toBe(2);
    expect(map.get("c")).toBe(4);
})

test("List slicing for ranges", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);

    let sliced = ll.sliceRange(1, 3);
    expect(sliced.length()).toEqual(2);
    expect(sliced.values()).toEqual([2, 3]);
    expect(sliced.names()).toEqual(["B", "C"]);

    let unnamed = new bioc.List([1,2,3,4,5]);
    let usliced = unnamed.sliceRange(0, 4);
    expect(usliced.values()).toEqual([1,2,3,4]);
    expect(usliced.names()).toBeNull();

    // Works in place.
    ll.sliceRange(2, 3, { inPlace: true });
    expect(ll.length()).toEqual(1);
    expect(ll.values()).toEqual([3]);
    expect(ll.names()).toEqual(["C"]);
})

test("List slicing for indices", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);

    {
        let sliced = ll.sliceIndices([1, 3]);
        expect(sliced.length()).toEqual(2);
        expect(sliced.values()).toEqual([2, 4]);
        expect(sliced.names()).toEqual(["B", "D"]);

        sliced = ll.sliceIndices(["D", "C", "A"]);
        expect(sliced.length()).toEqual(3);
        expect(sliced.values()).toEqual([4, 3, 1]);
        expect(sliced.names()).toEqual(["D", "C", "A"]);
    }

    {
        let unnamed = new bioc.List([1,2,3,4,5]);
        let usliced = unnamed.sliceIndices([0, 4, 3]);
        expect(usliced.values()).toEqual([1,5,4]);
        expect(usliced.names()).toBeNull();

        expect(() => unnamed.sliceIndices(["D", "C", "A"])).toThrow("available names");
    }

    // Works with duplicates.
    {
        let dupped = new bioc.List([1,2,3,4,5], { names: ["a", "b", "a", "c", "c"] });
        let sliced = dupped.sliceIndices(["a", "c", "b"]);
        expect(sliced.values()).toEqual([1,4,2]);
        expect(sliced.names()).toEqual(["a","c","b"]);
    }

    // Works in place.
    expect(ll.get("A")).toEqual(1); // populating the lookup.
    ll.sliceIndices([2, 3], { inPlace: true });
    expect(ll.length()).toEqual(2);
    expect(ll.values()).toEqual([3,4]);
    expect(() => ll.get("A")).toThrow("no matching name"); // check that the lookup is properly reset.
})

test("Lists can be iterated over", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);

    let collected = [];
    for (const x of ll) {
        collected.push(x);
    }
    expect(collected).toEqual([1,2,3,4,5]);
})

test("LENGTH works for List", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);
    expect(bioc.LENGTH(ll)).toEqual(5);
})

test("SLICE works for List", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);
    let slice = bioc.SLICE(ll, [ 0, 2, 4 ]);
    expect(slice.values()).toEqual([1, 3, 5]);
    expect(slice.names()).toEqual(["A", "C", "E"]);
})

test("CLONE works for List", () => {
    let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
    let ll = new bioc.List(src);
    let deepcopy = bioc.CLONE(ll);

    deepcopy.values()[0] = 3;
    expect(deepcopy.get(0)).toEqual(3);
    expect(ll.get(0)).toEqual(1); // original is unaffected.
})

test("COMBINE works for List", () => {
    // All named.
    {
        let src = { A: 1, B: 2, C: 3, D: 4, E: 5 };
        let ll = new bioc.List(src);

        let src2 = { a: -1, b: -2, c: -3, d: -4, e: -5 };
        let ll2 = new bioc.List(src2);
        let combined = bioc.COMBINE([ll, ll2]);
        expect(combined.values()).toEqual([1,2,3,4,5,-1,-2,-3,-4,-5]);
        expect(combined.names()).toEqual(["A", "B", "C", "D", "E", "a", "b", "c", "d", "e"]);

        let src3 = { alpha: 100 };
        let ll3 = new bioc.List(src3);
        combined = bioc.COMBINE([ll, ll2, ll3]);
        expect(combined.values()).toEqual([1,2,3,4,5,-1,-2,-3,-4,-5, 100]);
        expect(combined.names()).toEqual(["A", "B", "C", "D", "E", "a", "b", "c", "d", "e", "alpha"]);
    }

    // All unnamed.
    {
        let unnamed = new bioc.List([1,2,3,4,5]);
        let unnamed2 = new bioc.List([-1,-2,-3,-4,-5]);
        let combined = bioc.COMBINE([unnamed, unnamed2]);
        expect(combined.values()).toEqual([1,2,3,4,5,-1,-2,-3,-4,-5]);
        expect(combined.names()).toBeNull();

        let unnamed3 = new bioc.List([100]);
        combined = bioc.COMBINE([unnamed, unnamed2, unnamed3]);
        expect(combined.values()).toEqual([1,2,3,4,5,-1,-2,-3,-4,-5, 100]);
        expect(combined.names()).toBeNull();
    }

    // Mix of named and unnamed.
    {
        let unnamed = new bioc.List([1,2,3,4,5]);
        let src2 = { a: -1, b: -2, c: -3, d: -4, e: -5 };
        let ll2 = new bioc.List(src2);
        let unnamed3 = new bioc.List([100]);

        let combined = bioc.COMBINE([unnamed, ll2, unnamed3]);
        expect(combined.values()).toEqual([1,2,3,4,5,-1,-2,-3,-4,-5, 100]);
        expect(combined.names()).toEqual(["", "", "", "", "", "a", "b", "c", "d", "e", ""]);
    }

    // Try other things.
    {
        let unnamed = new bioc.List([1,2,3,4,5]);
        let combined = bioc.COMBINE([unnamed, [-1,-2], { "x": 300, "y": 400 }]);
        expect(combined.values()).toEqual([1,2,3,4,5,-1,-2,300,400]);
        expect(combined.names()).toEqual(["", "", "", "", "", "", "", "x", "y"]);
    }
})
