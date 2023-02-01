import * as bioc from "../src/index.js";
import * as utils from "./utils.js";

test("Basic construction of a RangedSummarizedExperiment works", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let gr = utils.spawn_random_GRanges(10);
    let rse = new bioc.RangedSummarizedExperiment({ counts: mat }, gr);

    expect(rse.numberOfRows()).toBe(10);
    expect(rse.numberOfColumns()).toBe(20);
    expect(rse.rowRanges().start()).toEqual(gr.start());

    // Errors.
    expect(() => new bioc.RangedSummarizedExperiment({ counts: mat }, bioc.SLICE(gr, [1,2,3]))).toThrow("number of rows")
    expect(() => new bioc.RangedSummarizedExperiment({ counts: mat }, [1,2,3])).toThrow("GRanges")
})

test("Construction of a RangedSummarizedExperiment works with all bits and pieces", () => {
    let mat = utils.spawn_random_matrix(5, 4);
    let mat2 = utils.spawn_random_matrix(5, 4);
    let gr = utils.spawn_random_GRanges(5);
    let rdf = new bioc.DataFrame({ foo: utils.spawn_random_vector(5) });
    let cdf = new bioc.DataFrame({ bar: utils.spawn_random_vector(4) });
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
    expect(rse.metadata().get("bob")).toBe(1);
})

test("Construction of a RangedSummarizedExperiment works with null rowRanges", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let rse = new bioc.RangedSummarizedExperiment({ counts: mat }, null);

    expect(rse.rowRanges() instanceof bioc.GroupedGRanges).toBe(true);
    expect(bioc.LENGTH(rse.rowRanges())).toEqual(10);
    expect(rse.rowRanges().rangeStarts()[0]).toBe(0);
    expect(rse.rowRanges().rangeStarts()[9]).toBe(0);
})

test("setting/removing of the rowRanges works as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let gr = utils.spawn_random_GRanges(10);
    let rse = new bioc.RangedSummarizedExperiment({ counts: mat }, gr);

    let gr2 = utils.spawn_random_GRanges(10);
    let rse2 = rse.setRowRanges(gr2);
    expect(rse2.rowRanges().start()).toEqual(gr2.start());
    expect(rse.rowRanges().start()).toEqual(gr.start()); // original is still preserved.

    rse.$setRowRanges(gr2);
    expect(rse.rowRanges().start()).toEqual(gr2.start());
    expect(rse.rowRanges().start()).not.toEqual(gr.start());

    // Works with a GroupedGRanges.
    rse.$setRowRanges(bioc.GroupedGRanges.empty(10));
    expect(rse.rowRanges() instanceof bioc.GroupedGRanges).toBe(true);

    // Errors.
    expect(() => rse.$setRowRanges(bioc.SLICE(gr2, [1,2,3]))).toThrow("number of rows");
    expect(() => rse.$setRowRanges([1,2,3])).toThrow("GRanges");
})

test("SLICE_2D generic works as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let gr = utils.spawn_random_GRanges(10);
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

test("COMBINE_ROWS generic works as expected with GRanges", () => {
    let NC = 16;

    let NR1 = 11;
    let mat1 = utils.spawn_random_matrix(NR1, NC);
    let gr1 = utils.spawn_random_GRanges(NR1);
    let se1 = new bioc.RangedSummarizedExperiment({ counts: mat1 }, gr1);

    let NR2 = 9;
    let mat2 = utils.spawn_random_matrix(NR2, NC);
    let gr2 = utils.spawn_random_GRanges(NR2);
    let se2 = new bioc.RangedSummarizedExperiment({ counts: mat2 }, gr2);

    let combined = bioc.COMBINE_ROWS([se1, se2]);
    expect(combined.numberOfRows()).toEqual(NR1 + NR2);
    expect(combined.numberOfColumns()).toEqual(NC);
    expect(combined.rowRanges().start()).toEqual(bioc.COMBINE([gr1.start(), gr2.start()]));
})

test("COMBINE_ROWS generic works as expected with mixed GRanges and GroupedGRanges", () => {
    let NC = 16;

    let NR1 = 11;
    let mat1 = utils.spawn_random_matrix(NR1, NC);
    let gr1 = utils.spawn_random_GRanges(NR1);
    let se1 = new bioc.RangedSummarizedExperiment({ counts: mat1 }, gr1);

    let NR2 = 9;
    let mat2 = utils.spawn_random_matrix(NR2, NC);
    let width = new Int32Array(NR2);
    width.fill(2);
    let ggr2 = new bioc.GroupedGRanges(utils.spawn_random_GRanges(NR2*2), { rangeLengths: width });
    let se2 = new bioc.RangedSummarizedExperiment({ counts: mat2 }, ggr2);

    let combined = bioc.COMBINE_ROWS([se1, se2]);
    let comranges = combined.rowRanges();
    expect(comranges instanceof bioc.GroupedGRanges).toBe(true);
    expect(comranges.numberOfGroups()).toBe(NR1 + NR2);

    let first = bioc.SLICE(comranges, { start: 0, end: NR1 });
    expect(first.ranges().start()).toEqual(gr1.start());
    let filling = new Int32Array(NR1);
    filling.fill(1);
    expect(first.rangeLengths()).toEqual(filling);

    let second = bioc.SLICE(comranges, { start: NR1, end: NR1 + NR2 });
    expect(second.ranges().start()).toEqual(ggr2.ranges().start());
    expect(second.rangeLengths()).toEqual(ggr2.rangeLengths());
})

test("COMBINE_ROWS generic works as expected with mixed GRanges and empty objects", () => {
    let NC = 12;

    let NR1 = 5;
    let mat1 = utils.spawn_random_matrix(NR1, NC);
    let gr1 = utils.spawn_random_GRanges(NR1);
    let se1 = new bioc.RangedSummarizedExperiment({ counts: mat1 }, gr1);

    let NR2 = 9;
    let mat2 = utils.spawn_random_matrix(NR2, NC);
    let se2 = new bioc.SummarizedExperiment({ counts: mat2 });

    let combined = bioc.COMBINE_ROWS([se1, se2]);
    let comranges = combined.rowRanges();
    expect(comranges instanceof bioc.GroupedGRanges).toBe(true);
    expect(comranges.numberOfGroups()).toBe(NR1 + NR2);

    let first = bioc.SLICE(comranges, { start: 0, end: NR1 });
    expect(first.ranges().start()).toEqual(gr1.start());
    let filling = new Int32Array(NR1);
    filling.fill(1);
    expect(first.rangeLengths()).toEqual(filling);

    let second = bioc.SLICE(comranges, { start: NR1, end: NR1 + NR2 });
    expect(bioc.LENGTH(second.ranges())).toEqual(0);
})

test("COMBINE_COLUMNS generic works as expected", () => {
    let NR = 9

    let NC1 = 7;
    let mat1 = utils.spawn_random_matrix(NR, NC1);
    let gr1 = utils.spawn_random_GRanges(NR);
    let se1 = new bioc.RangedSummarizedExperiment({ counts: mat1 }, gr1);

    let NC2 = 8;
    let mat2 = utils.spawn_random_matrix(NR, NC2);
    let gr2 = utils.spawn_random_GRanges(NR);
    let se2 = new bioc.RangedSummarizedExperiment({ counts: mat2 }, gr2);

    let combined = bioc.COMBINE_COLUMNS([se1, se2]);
    expect(combined.numberOfRows()).toEqual(NR);
    expect(combined.numberOfColumns()).toEqual(NC1 + NC2);
    expect(combined.rowRanges().start()).toEqual(gr1.start());
})

test("CLONE generic works as expected", () => {
    let mat = utils.spawn_random_matrix(11, 13);
    let gr = utils.spawn_random_GRanges(11);
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
