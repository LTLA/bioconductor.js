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

test("setting/removing of the assays works as expected", () => {
    let mat = spawn_random_matrix(10, 20);
    let out = new bioc.SummarizedExperiment({ counts: mat });
    expect(out.assayNames()).toEqual(["counts"]);
    expect(out.assay(0).values()).toEqual(mat.values());

    let mat2 = spawn_random_matrix(10, 20);
    out.$setAssay(0, mat2);
    expect(out.assayNames()).toEqual(["counts"]);
    expect(out.assay(0).values()).toEqual(mat2.values());

    out.$setAssay("logcounts", mat);
    expect(out.assayNames()).toEqual(["counts", "logcounts"]);
    expect(out.assay(1).values()).toEqual(mat.values());

    out.$removeAssay(0);
    expect(out.assayNames()).toEqual(["logcounts"]);
    expect(out.assay(0).values()).toEqual(mat.values());

    out.$removeAssay("logcounts");
    expect(out.assayNames()).toEqual([]);

    // Errors.
    expect(() => { out.$setAssay('foo', spawn_random_matrix(10, 10)) }).toThrow("same dimensions");
    expect(() => { out.$setAssay('foo', spawn_random_matrix(20, 20)) }).toThrow("same dimensions");
    expect(() => { out.$removeAssay('foo') }).toThrow("no assay");
})

test("setting/removing of the DFs works as expected", () => {
    let mat = spawn_random_matrix(10, 20);
    let out = new bioc.SummarizedExperiment({ counts: mat });

    let rdf = new bioc.DataFrame({ foo: spawn_random_vector(10) });
    out.$setRowData(rdf);
    expect(out.rowData().column("foo")).toEqual(rdf.column("foo"));

    let cdf = new bioc.DataFrame({ bar: spawn_random_vector(20) });
    out.$setColumnData(cdf);
    expect(out.columnData().column("bar")).toEqual(cdf.column("bar"));

    // Errors.
    expect(() => { out.$setRowData(bioc.SLICE(rdf, { start: 5, end: 10 })) }).toThrow("same number of rows");
    expect(() => { out.$setColumnData(bioc.SLICE(rdf, { start: 5, end: 10 })) }).toThrow("same number of rows");
})

test("setting/removing of the names works as expected", () => {
    let mat = spawn_random_matrix(5, 4);
    let out = new bioc.SummarizedExperiment({ counts: mat });

    let rnames = [ "A", "B", "C", "D", "E" ];
    out.$setRowNames(rnames);
    expect(out.rowNames()).toEqual(rnames);

    let cnames = [ "a", "b", "c", "d" ];
    out.$setColumnNames(cnames);
    expect(out.columnNames()).toEqual(cnames);

    out.$setRowNames(null);
    out.$setColumnNames(null);
    expect(out.rowNames()).toBeNull();
    expect(out.columnNames()).toBeNull();

    // Errors.
    expect(() => { out.$setRowNames(["a", "b"]) }).toThrow("numberOfRows");
    expect(() => { out.$setColumnNames(["a", "b"]) }).toThrow("numberOfColumns");
})

test("NUMBER_OF generics work as expected", () => {
    let mat = spawn_random_matrix(11, 13);
    let out = new bioc.SummarizedExperiment({ counts: mat });
    expect(bioc.NUMBER_OF_ROWS(out)).toBe(11);
    expect(bioc.NUMBER_OF_COLUMNS(out)).toBe(13);
})

test("SLICE_2D generic works as expected", () => {
    let mat = spawn_random_matrix(11, 13);
    let rdf = new bioc.DataFrame({ foo: spawn_random_vector(11) });
    let cdf = new bioc.DataFrame({ bar: spawn_random_vector(13) });
    let out = new bioc.SummarizedExperiment({ counts: mat }, { rowData: rdf, columnData: cdf });

    // No-op.
    {
        let sliced = bioc.SLICE_2D(out, null, null);

        expect(sliced.numberOfRows()).toBe(out.numberOfRows());
        expect(sliced.rowData().column("foo")).toBe(rdf.column("foo"));
        expect(sliced.numberOfColumns()).toBe(out.numberOfColumns());
        expect(sliced.columnData().column("bar")).toBe(cdf.column("bar"));

        expect(sliced.assayNames()).toEqual(out.assayNames());
        expect(sliced.assay(0).values()).toEqual(mat.values());

        expect(sliced.rowNames()).toBeNull();
        expect(sliced.columnNames()).toBeNull();
    }

    // Rows only.
    {
        let ridx = [ 1, 1, 2, 3, 5, 8 ];
        let sliced = bioc.SLICE_2D(out, ridx, null);

        expect(sliced.numberOfRows()).toBe(ridx.length);
        expect(sliced.rowData().column("foo")).toEqual(bioc.SLICE(rdf.column("foo"), ridx));
        expect(sliced.numberOfColumns()).toBe(out.numberOfColumns());
        expect(sliced.columnData().column("bar")).toEqual(cdf.column("bar"));

        expect(sliced.assay(0).numberOfRows()).toEqual(ridx.length);
        expect(sliced.assay("counts").numberOfColumns()).toEqual(out.numberOfColumns());

        expect(sliced.rowNames()).toBeNull();
        expect(sliced.columnNames()).toBeNull();
    }

    // Columns only.
    {
        let cidx = [ 1, 1, 2, 3, 5, 8 ];
        let sliced = bioc.SLICE_2D(out, null, cidx);

        expect(sliced.numberOfRows()).toBe(out.numberOfRows());
        expect(sliced.rowData().column("foo")).toEqual(rdf.column("foo"));
        expect(sliced.numberOfColumns()).toBe(cidx.length);
        expect(sliced.columnData().column("bar")).toEqual(bioc.SLICE(cdf.column("bar"), cidx));

        expect(sliced.assay("counts").numberOfRows()).toEqual(out.numberOfRows());
        expect(sliced.assay(0).numberOfColumns()).toEqual(cidx.length);

        expect(sliced.rowNames()).toBeNull();
        expect(sliced.columnNames()).toBeNull();
    }

    // Both at the same time.
    {
        let sliced = bioc.SLICE_2D(out, { start: 5, end: 9 }, { start: 0, end: 10 });

        expect(sliced.numberOfRows()).toBe(4);
        expect(sliced.rowData().column("foo")).toEqual(rdf.column("foo").slice(5, 9));
        expect(sliced.numberOfColumns()).toBe(10);
        expect(sliced.columnData().column("bar")).toEqual(cdf.column("bar").slice(0, 10));

        expect(sliced.assay("counts").numberOfRows()).toEqual(4);
        expect(sliced.assay(0).numberOfColumns()).toEqual(10);

        expect(sliced.rowNames()).toBeNull();
        expect(sliced.columnNames()).toBeNull();
    }

    // Handles names and metadata properly.
    {
        let mat = spawn_random_matrix(6, 5);
        let se = new bioc.SummarizedExperiment({ foo: mat }, { 
            rowNames: [ "A", "B", "C", "D", "E", "F" ], 
            columnNames: [ "a", "b", "c", "d", "e" ],
            metadata: { wailord: 999 }
        });

        let sliced = bioc.SLICE_2D(se, [ 5, 3, 1 ], [ 0, 2, 4 ]);
        expect(sliced.rowNames()).toEqual(["F", "D", "B"]);
        expect(sliced.columnNames()).toEqual(["a", "c", "e"]);
        expect(sliced.metadata().wailord).toBe(999);
    }
})

function mock_names(n, prefix="X") {
    let output = [];
    for (var i = 0; i < n; i++) {
        output.push(prefix + String(i));
    }
    return output;
}

test("COMBINE_ROWS generic works as expected", () => {
    let NC = 16;

    let NR1 = 11;
    let mat1 = spawn_random_matrix(NR1, NC);
    let rdf1 = new bioc.DataFrame({ foo: spawn_random_vector(NR1) });
    let cdf1 = new bioc.DataFrame({ bar: spawn_random_vector(NC) });
    let se1 = new bioc.SummarizedExperiment({ counts: mat1 }, { rowData: rdf1, columnData: cdf1, metadata: { bob: "builder" } });

    let NR2 = 9;
    let mat2 = spawn_random_matrix(NR2, NC);
    let rdf2 = new bioc.DataFrame({ foo: spawn_random_vector(NR2) });
    let cdf2 = new bioc.DataFrame({ weyland: spawn_random_vector(NC) });
    let se2 = new bioc.SummarizedExperiment({ counts: mat2 }, { rowData: rdf2, columnData: cdf2, metadata: { second: "yutani" }});

    // Basic handling.
    {
        let combined = bioc.COMBINE_ROWS([se1, se2]);
        expect(combined.numberOfRows()).toEqual(NR1 + NR2);
        expect(combined.numberOfColumns()).toEqual(NC);

        expect(combined.rowData().column("foo")).toEqual(bioc.COMBINE([rdf1.column("foo"), rdf2.column("foo")]));
        expect(combined.columnData().column("bar")).toEqual(cdf1.column("bar"));
        expect(combined.columnData().hasColumn("weyland")).toBe(false);

        expect(combined.assayNames()).toEqual(["counts"]);
        expect(combined.assay("counts").values()).toEqual(bioc.COMBINE_ROWS([mat1, mat2]).values());

        expect(combined.metadata()).toEqual({ bob: "builder" });
    }

    // Multiple assays.
    {
        let lmat1 = spawn_random_matrix(NR1, NC);
        let lmat2 = spawn_random_matrix(NR2, NC);
        let se1 = new bioc.SummarizedExperiment({ counts: mat1, logcounts: lmat1 });
        let se2 = new bioc.SummarizedExperiment({ counts: mat2, logcounts: lmat2 });

        let combined = bioc.COMBINE_ROWS([se1, se2]);
        expect(combined.assayNames()).toEqual(["counts", "logcounts"]);
        expect(combined.assay(0).values()).toEqual(bioc.COMBINE_ROWS([mat1, mat2]).values());
        expect(combined.assay(1).values()).toEqual(bioc.COMBINE_ROWS([lmat1, lmat2]).values());
    }

    // With names.
    {
        let se1 = new bioc.SummarizedExperiment({ counts: mat1 }, { rowNames: mock_names(NR1), columnNames: mock_names(NC, "Y") });
        let se2 = new bioc.SummarizedExperiment({ counts: mat2 }, { rowNames: mock_names(NR2), columnNames: mock_names(NC, "Y") });

        let combined = bioc.COMBINE_ROWS([se1, se2]);
        expect(combined.rowNames()).toEqual(bioc.COMBINE([se1.rowNames(), se2.rowNames()]));
        expect(combined.columnNames()).toEqual(se1.columnNames());
    }
})

test("COMBINE_COLUMNS generic works as expected", () => {
    let NR = 9

    let NC1 = 7;
    let mat1 = spawn_random_matrix(NR, NC1);
    let rdf1 = new bioc.DataFrame({ foo: spawn_random_vector(NR) });
    let cdf1 = new bioc.DataFrame({ bar: spawn_random_vector(NC1) });
    let se1 = new bioc.SummarizedExperiment({ counts: mat1 }, { rowData: rdf1, columnData: cdf1, metadata: { bob: "builder" } });

    let NC2 = 8;
    let mat2 = spawn_random_matrix(NR, NC2);
    let rdf2 = new bioc.DataFrame({ weyland: spawn_random_vector(NR) });
    let cdf2 = new bioc.DataFrame({ bar: spawn_random_vector(NC2) });
    let se2 = new bioc.SummarizedExperiment({ counts: mat2 }, { rowData: rdf2, columnData: cdf2, metadata: { second: "yutani" }});

    // Basic handling.
    {
        let combined = bioc.COMBINE_COLUMNS([se1, se2]);
        expect(combined.numberOfRows()).toEqual(NR);
        expect(combined.numberOfColumns()).toEqual(NC1 + NC2);

        expect(combined.rowData().column("foo")).toEqual(rdf1.column("foo"));
        expect(combined.rowData().hasColumn("weyland")).toBe(false);
        expect(combined.columnData().column("bar")).toEqual(bioc.COMBINE([cdf1.column("bar"), cdf2.column("bar")]));

        expect(combined.assayNames()).toEqual(["counts"]);
        expect(combined.assay("counts").values()).toEqual(bioc.COMBINE_COLUMNS([mat1, mat2]).values());

        expect(combined.metadata()).toEqual({ bob: "builder" });
    }

    // Multiple assays.
    {
        let lmat1 = spawn_random_matrix(NR, NC1);
        let lmat2 = spawn_random_matrix(NR, NC2);
        let se1 = new bioc.SummarizedExperiment({ counts: mat1, logcounts: lmat1 });
        let se2 = new bioc.SummarizedExperiment({ counts: mat2, logcounts: lmat2 });

        let combined = bioc.COMBINE_COLUMNS([se1, se2]);
        expect(combined.assayNames()).toEqual(["counts", "logcounts"]);
        expect(combined.assay(0).values()).toEqual(bioc.COMBINE_COLUMNS([mat1, mat2]).values());
        expect(combined.assay(1).values()).toEqual(bioc.COMBINE_COLUMNS([lmat1, lmat2]).values());
    }

    // With names.
    {
        let se1 = new bioc.SummarizedExperiment({ counts: mat1 }, { rowNames: mock_names(NR), columnNames: mock_names(NC1, "Y") });
        let se2 = new bioc.SummarizedExperiment({ counts: mat2 }, { rowNames: mock_names(NR), columnNames: mock_names(NC2, "Y") });

        let combined = bioc.COMBINE_COLUMNS([se1, se2]);
        expect(combined.rowNames()).toEqual(se1.rowNames());
        expect(combined.columnNames()).toEqual(bioc.COMBINE([se1.columnNames(), se2.columnNames()]));
    }
})

test("CLONE generic works as expected", () => {
    let mat = spawn_random_matrix(11, 13);
    let rdf = new bioc.DataFrame({ foo: spawn_random_vector(11) });
    let cdf = new bioc.DataFrame({ bar: spawn_random_vector(13) });
    let out = new bioc.SummarizedExperiment({ counts: mat }, { rowData: rdf, columnData: cdf });

    // Deep copy
    let deep = bioc.CLONE(out);
    deep.assay(0).values()[0] = -100;
    expect(deep.assay(0).values()[0]).toEqual(-100);
    expect(out.assay(0).values()[0]).not.toEqual(-100);

    // Shallow copy.
    let shallow = bioc.CLONE(out, { deepCopy: false });
    shallow.assay(0).values()[0] = -100;
    expect(shallow.assay(0).values()[0]).toEqual(-100);
    expect(out.assay(0).values()[0]).toEqual(-100);
})
