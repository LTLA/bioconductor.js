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
    expect(() => new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1, MDS: rd2 }, reducedDimensionOrder: [] })).toThrow("reducedDimensionOrder");
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
    expect(() => new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { ADT: alt1, ERCC: alt2 }, alternativeExperimentOrder: [] })).toThrow ("alternativeExperimentOrder");
})

test("Reduced dimension setters work as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let rd1 = utils.spawn_random_matrix(20, 5);
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { reducedDimensions: { PCA: rd1 } });
    expect(sce.reducedDimensionNames()).toEqual(["PCA"]);

    let rd2 = utils.spawn_random_matrix(20, 5);
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
    expect(() => sce.$removeReducedDimension("FOO")).toThrow("no reduced dimension");
    expect(() => sce.$setReducedDimension("PCA", utils.spawn_random_matrix(10, 1))).toThrow("number of rows");
})

test("Alternative experiment setters work as expected", () => {
    let mat = utils.spawn_random_matrix(10, 20);
    let alt1 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(5, 20) });
    let sce = new bioc.SingleCellExperiment({ counts: mat }, { alternativeExperiments: { ADT: alt1 } });

    let alt2 = new bioc.SummarizedExperiment({ counts: utils.spawn_random_matrix(3, 20) });
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
    expect(() => sce.$removeAlternativeExperiment("FOO")).toThrow("no alternative experiment");
    expect(() => sce.$setAlternativeExperiment("X", utils.spawn_random_matrix(10, 1))).toThrow("SummarizedExperiment");
    expect(() => sce.$setAlternativeExperiment("X", bioc.SLICE_2D(alt1, null, [1,2,3]))).toThrow("number of columns");
})

