import * as generics from "./AllGenerics.js";
import * as rse from "./RangedSummarizedExperiment.js";
import * as utils from "./utils.js";

export class SingleCellExperiment extends rse.RangedSummarizedExperiment {
    static #check_reduced_dims(reducedDimensions, numberOfColumns) {
        for (const [k, v] of Object.entries(reducedDimensions)) {
            if (bioc.NUMBER_OF_ROWS(v) !== numberOfColumns) {
                throw new Error("number of rows for reduced dimension '" + k + "' is not equal to number of columns for this SingleCellExperiment");
            }
        }
    }

    static #check_alternative_experiments(alternativeExperiments, numberOfColumns) {
        for (const [k, v] of Object.entries(reducedDimensions)) {
            if (!(v instanceof se.SummarizedExperiment)) {
                throw new Error("alternative experiment '" + k + "' is not a SummarizedExperiment");
            }
            if (bioc.numberOfColumns(v) !== numberOfColumns) {
                throw new Error("number of columns for alternative experiment '" + k + "' is not equal to number of columns for this SingleCellExperiment");
            }
        }
    }

    /**
     * @param {Object} assays - Object where keys are the assay names and values are multi-dimensional arrays of experimental data.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?(GRanges|GroupedGRanges)} rowRanges - Genomic ranges corresponding to each row, see the {@linkplain RangedSummarizedExperiment} constructor.
     * @param {?Array} [options.assayOrder=null] - Array of strings specifying the ordering of the assays.
     * If non-`null`, this should have the same values as the keys of `assays`.
     * If `null`, an arbitrary ordering is obtained from `assays`.
     * @param {?DataFrame} [options.rowData=null] - Data frame of row annotations.
     * If non-`null`, this should have a number of rows equal to the number of rows in each entry of `assays`.
     * If `null`, an empty {@linkplain DataFrame} is automatically created.
     * @param {?DataFrame} [options.columnData=null] - Data frame of column annotations.
     * If non-`null`, this should have a number of columns equal to the number of columns in each entry of `assays`.
     * If `null`, an empty {@linkplain DataFrame} is automatically created.
     * @param {?Array} [options.rowNames=null] - Array of strings of length equal to the number of rows in the `assays`, containing row names.
     * Alternatively `null`, if no row names are present.
     * @param {?Array} [options.columnNames=null] - Array of strings of length equal to the number of columns in the `assays`, containing column names.
     * Alternatively `null`, if no column names are present.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(assays, { 
        reducedDimensions = null,
        reducedDimensionOrder = null,
        alternativeExperiments = null,
        alternativeExperimentOrder = null,
        rowRanges = null, 
        assayOrder = null, 
        rowData = null, 
        columnData = null, 
        rowNames = null, 
        columnNames = null, 
        metadata = {} 
    } = {}) {
        super(assays, rowRanges, { assayOrder, rowData, columnData, rowNames, columnNames, metadata });

        if (reducedDimensions === null) {
            reducedDimensions = {};
        }
        SingleCellExperiment.#check_reduced_dimensions(reducedDimensions, this.numberOfColumns());
        this._reducedDimensions = {
            entries: reducedDimensions,
            order: utils.checkEntryOrder(reducedDimensions, reducedDimensionOrder, "reducedDimension")
        };

        if (alternativeExperiments === null) {
            alternativeExperiments = {};
        }
        SingleCellExperiment.#check_alternative_experiments(alternativeExperiments, this.numberOfColumns());
        this._alternativeExperiments = {
            entries: alternativeExperiments,
            order: utils.checkEntryOrder(alternativeExperiments, alternativeExperimentOrder, "alternativeExperiment")
        };

        return;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    reducedDimensionNames() {
        return this._reducedDimensions.order;
    }

    reducedDimension(i) {
        return utils.retrieveSingleEntry(this._reducedDimensions.entries, this._reducedDimensions.order, i, "reduced dimension", "SingleCellExperiment");
    }

    alternativeExperimentNames() {
        return this._alternativeExperiments.order;
    }

    alternativeExperiment(i) {
        return utils.retrieveSingleEntry(this._alternativeExperiments.entries, this._alternativeExperiments.order, i, "alternative experiment", "SingleCellExperiment");
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    $removeReducedDimension(i) {
        utils.removeSingleEntry(this._reducedDimensions.entries, this._reducedDimensions.order, i, "reduced dimension", "SingleCellExperiment");
        return this;
    }

    $setReducedDimension(i, value) {
        utils.setSingleEntry(this._reducedDimensions.entries, this._reducedDimensions.order, i, value, "reduced dimension", "SingleCellExperiment");
        return this;
    }

    $removeAlternativeExperiment(i) {
        utils.removeSingleEntry(this._alternativeExperiments.entries, this._alternativeExperiments.order, i, "alternative experiment", "SingleCellExperiment");
        return this;
    }

    $setAlternativeExperiment(i, value) {
        utils.setSingleEntry(this._alternativeExperiments.entries, this._alternativeExperiments.order, i, value, "alternative experiment", "SingleCellExperiment");
        return this;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_SLICE_2D(output, rows, columns, { allowView = false }) {
        super._bioconductor_SLICE_2D(output, rows, columns, { allowView });

        output.reducedDimensions = utils.shallowCloneEntries(this._reducedDimensions);
        output.alternativeExperiments = utils.shallowCloneEntries(this.alternativeExperiments);

        if (columns !== null) {
            for (const [k, v] of Object.entries(output.reducedDimensions.entries)) {
                output.reducedDimensions.entries[k] = generics.SLICE_2D(v, columns, null, { allowView });
            }
            for (const [k, v] of Object.entries(output.alternativeExperiments.entries)) {
                output.alternativeExperiments.entries[k] = generics.SLICE_2D(v, null, columns, { allowView });
            }
        }
    }

    _bioconductor_COMBINE_ROWS(output, objects) {
        super._bioconductor_COMBINE_ROWS(output, objects);

        output.reducedDimensions = utils.shallowCloneEntries(this._reducedDimensions);
        output.alternativeExperiments = utils.shallowCloneEntries(this.alternativeExperiments);

        return;
    }

    _bioconductor_COMBINE_COLUMNS(output, objects) {
        super._bioconductor_COMBINE_COLUMNS(output, objects);

        output._reducedDimensions = utils.combineEntries(objects.map(x => x._reducedDimensions), generics.COMBINE_ROWS, "reducedDimNames" , "SingleCellExperiment"); 
        output._alternativeExperiments = utils.combineEntries(objects.map(x => x._alternativeExperiments), generics.COMBINE_COLUMNS, "alternativeExperimentNames" , "SingleCellExperiment"); 

        return;
    }

    _bioconductor_CLONE(output, { deepCopy }) {
        super._bioconductor_CLONE(output, { deepCopy });

        let FUN = (deepCopy ? generics.CLONE : utils.shallowCloneEntries);
        output.reducedDimensions = FUN(this._reducedDimensions);
        output.alternativeExperiments = FUN(this.alternativeExperiments);

        return;
    }

    

}
