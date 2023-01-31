import * as bioc from "../src/index.js";
import * as utils from "./utils.js";

test("Construction of a SingleCellExperiment works with reduced dimensions", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let rd1 = utils.spawn_random_matrix(20, 5);
    let rd2 = utils.spawn_random_matrix(20, 2);
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1, MDS: rd2 } });

    expect(sce.numberOfRows()).toBe(10);
    expect(sce.numberOfColumns()).toBe(20);

    expect(sce.reducedDimensionNames()).toEqual(["PCA", "MDS"]);
    expect(sce.reducedDimension(0).column(0)).toEqual(rd1.column(0));
    expect(sce.reducedDimension("MDS").column(1)).toEqual(rd2.column(1));

    // Flipping the order.
    let sce2 = new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1, MDS: rd2 }, reducedDimensionOrder: [ "MDS", "PCA"] });
    expect(sce2.reducedDimensionNames()).toEqual(["MDS", "PCA"]);
    expect(sce2.reducedDimension(0).column(0)).toEqual(rd2.column(0));

    // Errors.
    expect(() => new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1, MDS: bioc.SLICE_2D(rd2, [1,2,3], null) } })).toThrow("number of rows");
    expect(() => new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1, MDS: rd2 }, reducedDimensionOrder: [] })).toThrow("reduced dimension list");
})

test("Construction of a SingleCellExperiment works with alternative experiments", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let alt1 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) });
    let alt2 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(3, 20) });
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { ADT: alt1, ERCC: alt2 } });

    expect(sce.numberOfRows()).toBe(10);
    expect(sce.numberOfColumns()).toBe(20);

    expect(sce.alternativeExperimentNames()).toEqual(["ADT", "ERCC"]);
    expect(sce.alternativeExperiment(0).assay(0).column(0)).toEqual(alt1.assay(0).column(0));
    expect(sce.alternativeExperiment("ERCC").assay(0).column(0)).toEqual(alt2.assay(0).column(0));

    // Flipping the order.
    let sce2 = new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { ADT: alt1, ERCC: alt2 }, alternativeExperimentOrder: [ "ERCC", "ADT" ] });
    expect(sce2.alternativeExperimentNames()).toEqual(["ERCC", "ADT"]);
    expect(sce2.alternativeExperiment(0).assay(0).column(0)).toEqual(alt2.assay(0).column(0));

    // Errors.
    expect(() => new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { ADT: alt1, ERCC: 2 } })).toThrow("SummarizedExperiment");
    expect(() => new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { ADT: alt1, ERCC: bioc.SLICE_2D(alt2, null, [1,2,3]) } })).toThrow("number of columns");
    expect(() => new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { ADT: alt1, ERCC: alt2 }, alternativeExperimentOrder: [] })).toThrow ("alternative experiment list");
})

test("Construction of a SingleCellExperiment works with all bits and pieces", () => {
    let mat = utils.spawn_random_matrix(5, 4);
    let mat2 = utils.spawn_random_matrix(5, 4);
    let gr = utils.spawn_random_GRanges(5);
    let rdf = new bioc.DataFrame({ foo: utils.spawn_random_vector(5) });
    let cdf = new bioc.DataFrame({ bar: utils.spawn_random_vector(4) });
    let rnames = [ "A", "B", "C", "D", "E" ];
    let cnames = [ "a", "b", "c", "d" ];

    let sce = new bioc.SingleCellExperiment(
        { counts: mat, logcounts: mat2 }, 
        { 
            assayOrder: [ "logcounts", "counts" ],
            rowRanges: gr, 
            rowData: rdf, 
            columnData: cdf, 
            rowNames: rnames, 
            columnNames: cnames, 
            metadata: { bob: 1 }
        }
    );

    expect(sce.assayNames()).toEqual(["logcounts", "counts"]);
    expect(sce.rowData().column("foo")).toEqual(rdf.column("foo"));
    expect(sce.columnData().column("bar")).toEqual(cdf.column("bar"));
    expect(sce.rowNames()).toEqual(rnames);
    expect(sce.columnNames()).toEqual(cnames);
    expect(sce.metadata().get("bob")).toBe(1);
})

test("Reduced dimension setters work as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let rd1 = utils.spawn_random_matrix(20, 5);
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1 } });
    expect(sce.reducedDimensionNames()).toEqual(["PCA"]);

    let rd2 = utils.spawn_random_matrix(20, 5);
    
    {
        // Originals are unaffected.
        let sce2 = sce.setReducedDimension("MDS", rd2);
        expect(sce2.reducedDimensionNames()).toEqual(["PCA", "MDS"]);
        expect(sce.reducedDimensionNames()).toEqual(["PCA"]); 

        sce2 = sce.removeReducedDimension("PCA");
        expect(sce2.reducedDimensionNames()).toEqual([]);
        expect(sce.reducedDimensionNames()).toEqual(["PCA"]); 
    }

    sce.$setReducedDimension("MDS", rd2);
    expect(sce.reducedDimensionNames()).toEqual(["PCA", "MDS"]);
    expect(sce.reducedDimension(1).numberOfColumns()).toEqual(5);

    let rd3 = utils.spawn_random_matrix(20, 2);
    sce.$setReducedDimension(1, rd3);
    expect(sce.reducedDimension(1).numberOfColumns()).toEqual(2);
    sce.$setReducedDimension("MDS", rd2);
    expect(sce.reducedDimension(1).numberOfColumns()).toEqual(5);

    sce.$removeReducedDimension(0);
    expect(sce.reducedDimensionNames()).toEqual(["MDS"]);
    sce.$removeReducedDimension("MDS");
    expect(sce.reducedDimensionNames()).toEqual([]);

    // Errors.
    expect(() => sce.$setReducedDimension(0, rd1)).toThrow("out of range");
    expect(() => sce.$removeReducedDimension(0)).toThrow("out of range");
    expect(() => sce.$removeReducedDimension("FOO")).toThrow("failed to remove");
    expect(() => sce.$setReducedDimension("PCA", utils.spawn_random_matrix(10, 1))).toThrow("number of rows");
})

test("Reduced dimension renaming works as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let rd1 = utils.spawn_random_matrix(20, 5);
    let rd2 = utils.spawn_random_matrix(20, 3);
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1, MDS: rd2 } });

    let out = sce.setReducedDimensionNames(["TSNE", "UMAP"]);
    expect(out.reducedDimensionNames()).toEqual(["TSNE", "UMAP"]);
    expect(sce.reducedDimensionNames()).toEqual(["PCA", "MDS"]);

    sce.$setReducedDimensionNames(["pca", "mds"]);
    expect(sce.reducedDimensionNames()).toEqual(["pca", "mds"]);
    expect(sce.reducedDimension(0).column(0)).toEqual(rd1.column(0));
    expect(sce.reducedDimension(1).column(0)).toEqual(rd2.column(0));

    expect(() => sce.$setReducedDimensionNames([1,2])).toThrow("failed to set the reduced dimension names");
})

test("Reduced dimension slicing works as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let rd1 = utils.spawn_random_matrix(20, 5);
    let rd2 = utils.spawn_random_matrix(20, 3);
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1, MDS: rd2 } });

    let out = sce.sliceReducedDimensions([1, 0]);
    expect(out.reducedDimensionNames()).toEqual(["MDS", "PCA"]);
    expect(out.reducedDimension(0).column(0)).toEqual(rd2.column(0));
    expect(out.reducedDimension(1).column(0)).toEqual(rd1.column(0));
    expect(sce.reducedDimensionNames()).toEqual(["PCA", "MDS"]); // original is unaffected.

    sce.$sliceReducedDimensions(["MDS"]);
    expect(sce.reducedDimensionNames()).toEqual(["MDS"]);
    expect(sce.reducedDimension(0).column(0)).toEqual(rd2.column(0));

    expect(() => sce.$sliceReducedDimensions([1,2])).toThrow("failed to slice the reduced dimensions");
})

test("Alternative experiment setters work as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let alt1 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) });
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { ADT: alt1 } });

    let alt2 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(3, 20) });

    {
        // Originals are unaffected.
        let sce2 = sce.setAlternativeExperiment("XXXX", alt2);
        expect(sce2.alternativeExperimentNames()).toEqual(["ADT", "XXXX"]);
        expect(sce.alternativeExperimentNames()).toEqual(["ADT"]); 

        sce2 = sce.removeAlternativeExperiment("ADT");
        expect(sce2.alternativeExperimentNames()).toEqual([]);
        expect(sce.alternativeExperimentNames()).toEqual(["ADT"]); 
    }

    sce.$setAlternativeExperiment("SPIKE", alt2);
    expect(sce.alternativeExperimentNames()).toEqual(["ADT", "SPIKE"]);
    expect(sce.alternativeExperiment(1).numberOfRows()).toEqual(3);

    let alt3 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(11, 20) });
    sce.$setAlternativeExperiment(1, alt3);
    expect(sce.alternativeExperiment(1).numberOfRows()).toEqual(11);
    sce.$setAlternativeExperiment("SPIKE", alt2);
    expect(sce.alternativeExperiment(1).numberOfRows()).toEqual(3);

    sce.$removeAlternativeExperiment(0);
    expect(sce.alternativeExperimentNames()).toEqual(["SPIKE"]);
    sce.$removeAlternativeExperiment("SPIKE");
    expect(sce.alternativeExperimentNames()).toEqual([]);

    // Errors.
    expect(() => sce.$setAlternativeExperiment(0, alt1)).toThrow("out of range");
    expect(() => sce.$removeAlternativeExperiment(0)).toThrow("out of range");
    expect(() => sce.$removeAlternativeExperiment("FOO")).toThrow("failed to remove");
    expect(() => sce.$setAlternativeExperiment("X", utils.spawn_random_matrix(10, 1))).toThrow("SummarizedExperiment");
    expect(() => sce.$setAlternativeExperiment("X", bioc.SLICE_2D(alt1, null, [1,2,3]))).toThrow("number of columns");
})

test("Alternative experiment renaming works as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let alt1 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) });
    let alt2 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(2, 20) });
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { foo: alt1, bar: alt2 } });

    let out = sce.setAlternativeExperimentNames(["WHEE", "stuff"]);
    expect(out.alternativeExperimentNames()).toEqual(["WHEE", "stuff"]);
    expect(sce.alternativeExperimentNames()).toEqual(["foo", "bar"]);

    sce.$setAlternativeExperimentNames(["blah", "whee"]);
    expect(sce.alternativeExperimentNames()).toEqual(["blah", "whee"]);
    expect(out.alternativeExperiment(0).assay(0).column(0)).toEqual(alt1.assay(0).column(0));
    expect(out.alternativeExperiment(1).assay(0).column(0)).toEqual(alt2.assay(0).column(0));

    expect(() => sce.$setAlternativeExperimentNames([1,2])).toThrow("failed to set the alternative experiment names");
})

test("Alternative experiment slicing works as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let alt1 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) });
    let alt2 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(2, 20) });
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { foo: alt1, bar: alt2 } });

    let out = sce.sliceAlternativeExperiments([1, 0]);
    expect(out.alternativeExperimentNames()).toEqual(["bar", "foo"]);
    expect(out.alternativeExperiment(0).assay(0).column(0)).toEqual(alt2.assay(0).column(0));
    expect(out.alternativeExperiment(1).assay(0).column(0)).toEqual(alt1.assay(0).column(0));
    expect(sce.alternativeExperimentNames()).toEqual(["foo", "bar"]);

    sce.$sliceAlternativeExperiments(["bar"]);
    expect(sce.alternativeExperimentNames()).toEqual(["bar"]);
    expect(sce.alternativeExperiment(0).assay(0).column(0)).toEqual(alt2.assay(0).column(0));

    expect(() => sce.$sliceAlternativeExperiments([1,2])).toThrow("failed to slice the alternative experiments");
})

test("SLICE_2D generic works as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let alt = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) });
    let rd = utils.spawn_random_matrix(20, 5);
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd }, alternativeExperiments: { STUFF: alt } });

    // By column.
    {
        let sliced = bioc.SLICE_2D(sce, null, [1,2,3,4]);
        expect(sliced.reducedDimension(0).numberOfRows()).toBe(4);
        expect(sliced.alternativeExperiment(0).numberOfColumns()).toBe(4);
    }

    // By row (no effect).
    {
        let sliced = bioc.SLICE_2D(sce, [1,2,3,4], null);
        expect(sliced.reducedDimension(0).numberOfRows()).toBe(20);
        expect(sliced.alternativeExperiment(0).numberOfColumns()).toBe(20);
    }
})

test("COMBINE_ROWS generic works as expected", () => {
    let mat1 = utils.spawn_random_matrix(10, 20);
    let alt1 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) });
    let rd1 = utils.spawn_random_matrix(20, 5);
    let sce1 = new bioc.SingleCellExperiment({ counts: mat1 }, { reducedDimensions: { PCA: rd1 }, alternativeExperiments: { STUFF: alt1 } });

    let se2 = new bioc.RangedSummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) }, null);
    let combined = bioc.COMBINE_ROWS([sce1, se2]);
    expect(combined.numberOfRows()).toEqual(15);

    // No-op on the reduced dimensions.
    expect(combined.reducedDimensionNames()).toEqual(["PCA"]);
    expect(combined.reducedDimension(0).column(0)).toEqual(rd1.column(0));
    expect(combined.alternativeExperimentNames()).toEqual(["STUFF"]);
    expect(combined.alternativeExperiment(0).assay(0).column(0)).toEqual(alt1.assay(0).column(0));
})

test("COMBINE_COLUMNS generic works as expected", () => {
    let mat1 = utils.spawn_random_matrix(10, 20);
    let alt1 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) });
    let rd1 = utils.spawn_random_matrix(20, 5);
    let sce1 = new bioc.SingleCellExperiment({ counts: mat1 }, { reducedDimensions: { PCA: rd1 }, alternativeExperiments: { STUFF: alt1 } });

    let mat2 = utils.spawn_random_matrix(10, 30);
    let alt2 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 30) });
    let rd2 = utils.spawn_random_matrix(30, 5);
    let sce2 = new bioc.SingleCellExperiment({ counts: mat2 }, { reducedDimensions: { PCA: rd2 }, alternativeExperiments: { STUFF: alt2 } });

    let combined = bioc.COMBINE_COLUMNS([ sce1, sce2 ]);
    expect(combined.numberOfColumns()).toEqual(50);
    expect(combined.reducedDimensionNames()).toEqual(["PCA"]);
    expect(combined.reducedDimension(0).column(0)).toEqual(bioc.COMBINE([rd1.column(0), rd2.column(0)]));
    expect(combined.alternativeExperimentNames()).toEqual(["STUFF"]);
    expect(combined.alternativeExperiment(0).assay(0).row(0)).toEqual(bioc.COMBINE([alt1.assay(0).row(0), alt2.assay(0).row(0)]));

    // Errors.
    let sce3 = new bioc.SingleCellExperiment({ counts: mat1 }, { reducedDimensions: { PCA: rd1 } });
    expect(() => bioc.COMBINE_COLUMNS([ sce1, sce3 ])).toThrow("failed to combine alternative experiments");
    sce3.$removeReducedDimension(0);
    expect(() => bioc.COMBINE_COLUMNS([ sce1, sce3 ])).toThrow("failed to combine reduced dimensions");
})

