import * as bioc from "../src/index.js";
import * as utils from "./utils.js";

test("Basic construction of a SingleCellExperiment works with reduced dimensions", () => {
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

test("Basic construction of a SingleCellExperiment works with alternative experiments", () => {
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
