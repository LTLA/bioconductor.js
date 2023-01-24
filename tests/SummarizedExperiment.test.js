import * as bioc from "../src/index.js";

function spawn_random_vector(N) {
    let v = new Float64Array(N);
    v.forEach((x, i) => { v[i] = Math.random(); });
    return v;
}

function spawn_random_matrix(NR, NC) {
    let v = spawn_random_vector(NR * NC);
    return new bioc.DenseMatrix(NR, NC, v);
}

test("Basic construction of a SummarizedExperiment works", () => {
    let mat = spawn_random_matrix(10, 20);
    let out = new bioc.SummarizedExperiment({ counts: mat });
    expect(out.numberOfRows()).toEqual(10);
    expect(out.numberOfColumns()).toEqual(20);

    expect(out.assayNames()).toEqual(["counts"]);
    expect(out.assay(0).values()).toEqual(mat.values());
    expect(out.assay("counts").values()).toEqual(mat.values());

    expect(out.rowData().numberOfRows()).toBe(10);
    expect(out.rowData().numberOfColumns()).toBe(0);
    expect(out.columnData().numberOfRows()).toBe(20);
    expect(out.columnData().numberOfColumns()).toBe(0);

    expect(out.rowNames()).toBeNull()
    expect(out.columnNames()).toBeNull()
})

test("Construction of a SummarizedExperiment works with DFs", () => {
    let mat = spawn_random_matrix(10, 20);
    let rdf = new bioc.DataFrame({ foo: spawn_random_vector(10) });
    let cdf = new bioc.DataFrame({ bar: spawn_random_vector(20) });

    let se = new bioc.SummarizedExperiment({ counts: mat }, { rowData: rdf, columnData: cdf });
    expect(se.rowData().column("foo")).toEqual(rdf.column("foo"));
    expect(se.columnData().column("bar")).toEqual(cdf.column("bar"));

    // Errors.
    expect(() => new bioc.SummarizedExperiment({ counts: mat }, { rowData: rdf, columnData: rdf })).toThrow("should be equal to the number of columns");
    expect(() => new bioc.SummarizedExperiment({ counts: mat }, { rowData: cdf, columnData: cdf })).toThrow("should be equal to the number of rows");

    // DFs supply dimensions.
    let none = new bioc.SummarizedExperiment({}, { rowData: rdf, columnData: cdf });
    expect(none.numberOfRows()).toBe(10);
    expect(none.numberOfColumns()).toBe(20);
    expect(none.assayNames()).toEqual([]);

    expect(() => new bioc.SummarizedExperiment({}, { columnData: cdf })).toThrow("rowData");
    expect(() => new bioc.SummarizedExperiment({}, { rowData: rdf })).toThrow("columnData");
})

test("Construction of a SummarizedExperiment works with multiple assays", () => {
    let assays = { counts: spawn_random_matrix(10, 20), logcounts: spawn_random_matrix(10, 20) };
    let se = new bioc.SummarizedExperiment(assays);

    expect(se.numberOfRows()).toEqual(10);
    expect(se.numberOfColumns()).toEqual(20);
    expect(se.assayNames()).toEqual(["counts", "logcounts"]);
    expect(se.assay('counts').values()).toEqual(assays["counts"].values());
    expect(se.assay(1).values()).toEqual(assays["logcounts"].values());

    // Works with custom ordering.
    {
        let se = new bioc.SummarizedExperiment(assays, { assayOrder: [ "logcounts", "counts" ] });
        expect(se.assayNames()).toEqual([ "logcounts", "counts" ]);
        expect(se.assay(0).values()).toEqual(assays["logcounts"].values());
        expect(se.assay("counts").values()).toEqual(assays["counts"].values());
    }

    // Errors with mismatching assays.
    assays["logcounts"] = spawn_random_matrix(10, 30);
    expect(() => new bioc.SummarizedExperiment(assays)).toThrow("same number of rows and columns");
    assays["logcounts"] = spawn_random_matrix(30, 20);
    expect(() => new bioc.SummarizedExperiment(assays)).toThrow("same number of rows and columns");
})

test("Construction of a SummarizedExperiment works with names", () => {
    let mat = spawn_random_matrix(5, 4);
    let rnames = [ "A", "B", "C", "D", "E" ];
    let cnames = [ "a", "b", "c", "d" ];

    let out = new bioc.SummarizedExperiment({ counts: mat }, { rowNames: rnames, columnNames: cnames });
    expect(out.numberOfRows()).toEqual(5);
    expect(out.numberOfColumns()).toEqual(4);
    expect(out.rowNames()).toEqual(rnames);
    expect(out.columnNames()).toEqual(cnames);

    // Errors.
    expect(() => new bioc.SummarizedExperiment({ counts: mat }, { rowNames: cnames, columnNames: cnames })).toThrow("number of rows");
    expect(() => new bioc.SummarizedExperiment({ counts: mat }, { rowNames: rnames, columnNames: rnames })).toThrow("number of columns");
})

test("Construction of a SummarizedExperiment works with metadata", () => {
    let mat = spawn_random_matrix(5, 4);
    let se = new bioc.SummarizedExperiment({ counts: mat }, { metadata: { foo: 1 } });
    expect(se.metadata().foo).toEqual(1);
})

