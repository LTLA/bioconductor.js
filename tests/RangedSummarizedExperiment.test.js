import * as bioc from "../src/index.js";

function spawn_random_vector(N) {
    let v = new Float64Array(N);
    v.forEach((x, i) => { v[i] = Math.random(); });
    return v;
}

function spawn_random_GRanges(N) {
    let s = new Int32Array(N);
    s.forEach((x, i) => { s[i] = Math.random() * 1000; });

    let w = new Int32Array(N);
    w.forEach((x, i) => { w[i] = Math.random() * 100; });

    let n = new Array(N);
    for (var i = 0; i < N; ++i) {
        n[i] = "chr" + String(Math.random() * 3 + 1);
    }

    return new bioc.GRanges(n, new bioc.IRanges(s, w));
}

function spawn_random_matrix(NR, NC) {
    let v = spawn_random_vector(NR * NC);
    return new bioc.DenseMatrix(NR, NC, v);
}

test("Basic construction of a RangedSummarizedExperiment works", () => {
    let mat = spawn_random_matrix(10, 20);
    let gr = spawn_random_GRanges(10);
    let rse = new bioc.RangedSummarizedExperiment({ counts: mat }, gr);

    expect(rse.numberOfRows()).toBe(10);
    expect(rse.numberOfColumns()).toBe(20);
    expect(rse.rowRanges().start()).toEqual(gr.start());

    // Errors.
    expect(() => new bioc.RangedSummarizedExperiment({ counts: mat }, bioc.SLICE(gr, [1,2,3]))).toThrow("number of rows")
    expect(() => new bioc.RangedSummarizedExperiment({ counts: mat }, [1,2,3])).toThrow("GRanges")
})

test("Construction of a RangedSummarizedExperiment works with all bits and pieces", () => {
    let mat = spawn_random_matrix(5, 4);
    let mat2 = spawn_random_matrix(5, 4);
    let gr = spawn_random_GRanges(5);
    let rdf = new bioc.DataFrame({ foo: spawn_random_vector(5) });
    let cdf = new bioc.DataFrame({ bar: spawn_random_vector(4) });
    let rnames = [ "A", "B", "C", "D", "E" ];
    let cnames = [ "a", "b", "c", "d" ];

    let rse = new bioc.RangedSummarizedExperiment(
        { counts: mat, logcounts: mat2 }, 
        gr, 
        { 
            assayOrder: [ "logcounts", "counts" ],
            rowData: rdf, 
            columnData: cdf, 
            rowNames: rnames, 
            columnNames: cnames, 
            metadata: { bob: 1 }
        }
    );

    expect(rse.assayNames()).toEqual(["logcounts", "counts"]);
    expect(rse.rowData().column("foo")).toEqual(rdf.column("foo"));
    expect(rse.columnData().column("bar")).toEqual(cdf.column("bar"));
    expect(rse.rowNames()).toEqual(rnames);
    expect(rse.columnNames()).toEqual(cnames);
    expect(rse.metadata().bob).toBe(1);
})

test("setting/removing of the rowRanges works as expected", () => {
    let mat = spawn_random_matrix(10, 20);
    let gr = spawn_random_GRanges(10);
    let rse = new bioc.RangedSummarizedExperiment({ counts: mat }, gr);

    let gr2 = spawn_random_GRanges(10);
    rse.$setRowRanges(gr2);
    expect(rse.rowRanges().start()).toEqual(gr2.start());
    expect(rse.rowRanges().start()).not.toEqual(gr.start());

    // Errors.
    expect(() => rse.$setRowRanges(bioc.SLICE(gr2, [1,2,3]))).toThrow("number of rows");
    expect(() => rse.$setRowRanges([1,2,3])).toThrow("GRanges");
})

test("SLICE_2D generic works as expected", () => {
    let mat = spawn_random_matrix(10, 20);
    let gr = spawn_random_GRanges(10);
    let rse = new bioc.RangedSummarizedExperiment({ counts: mat }, gr);

    // No-op.
    {
        let sliced = bioc.SLICE_2D(rse, null, null);
        expect(sliced.rowRanges().start()).toEqual(gr.start());
    }

    // Rows only.
    {
        let ridx = [ 1, 1, 2, 3, 5, 8 ];
        let sliced = bioc.SLICE_2D(rse, ridx, null);
        expect(sliced.rowRanges().start()).toEqual(bioc.SLICE(gr.start(), ridx));
    }

    // Columns only.
    {
        let cidx = [ 1, 1, 2, 3, 5, 8 ];
        let sliced = bioc.SLICE_2D(rse, null, cidx);
        expect(sliced.rowRanges().start()).toEqual(gr.start());
    }
})

test("COMBINE_ROWS generic works as expected", () => {
    let NC = 16;

    let NR1 = 11;
    let mat1 = spawn_random_matrix(NR1, NC);
    let gr1 = spawn_random_GRanges(NR1);
    let se1 = new bioc.RangedSummarizedExperiment({ counts: mat1 }, gr1);

    let NR2 = 9;
    let mat2 = spawn_random_matrix(NR2, NC);
    let gr2 = spawn_random_GRanges(NR2);
    let se2 = new bioc.RangedSummarizedExperiment({ counts: mat2 }, gr2);

    let combined = bioc.COMBINE_ROWS([se1, se2]);
    expect(combined.numberOfRows()).toEqual(NR1 + NR2);
    expect(combined.numberOfColumns()).toEqual(NC);
    expect(combined.rowRanges().start()).toEqual(bioc.COMBINE([gr1.start(), gr2.start()]));
})

test("COMBINE_COLUMNS generic works as expected", () => {
    let NR = 9

    let NC1 = 7;
    let mat1 = spawn_random_matrix(NR, NC1);
    let gr1 = spawn_random_GRanges(NR);
    let se1 = new bioc.RangedSummarizedExperiment({ counts: mat1 }, gr1);

    let NC2 = 8;
    let mat2 = spawn_random_matrix(NR, NC2);
    let gr2 = spawn_random_GRanges(NR);
    let se2 = new bioc.SummarizedExperiment({ counts: mat2 }, gr2);

    let combined = bioc.COMBINE_COLUMNS([se1, se2]);
    expect(combined.numberOfRows()).toEqual(NR);
    expect(combined.numberOfColumns()).toEqual(NC1 + NC2);
    expect(combined.rowRanges().start()).toEqual(gr1.start());
})

test("CLONE generic works as expected", () => {
    let mat = spawn_random_matrix(11, 13);
    let gr = spawn_random_GRanges(11);
    let out = new bioc.RangedSummarizedExperiment({ counts: mat }, gr);

    // Deep copy
    let deep = bioc.CLONE(out);
    deep.rowRanges().start()[0] = 100;
    expect(deep.rowRanges().start()[0]).toEqual(100);
    expect(out.rowRanges().start()[0]).not.toEqual(100);

    // Shallow copy.
    let shallow = bioc.CLONE(out, { deepCopy: false });
    shallow.rowRanges().start()[0] = 100;
    expect(shallow.rowRanges().start()[0]).toEqual(100);
    expect(out.rowRanges().start()[0]).toEqual(100);
})
