import * as generics from "./AllGenerics.js";
import * as rse from "./RangedSummarizedExperiment.js";
import * as se from "./SummarizedExperiment.js";
import * as utils from "./utils.js";

/**
 * A SingleCellExperiment is a {@linkplain RangedSummarizedExperiment} subclass that contains additional fields for storing reduced dimensions and alternative experiments.
 * It supports the same set of generics as the {@linkplain SummarizedExperiment}.
 *
 * Each reduced dimension instance should have number of rows equal to the number of columns of the SingleCellExperiment.
 * Each instance is expected to provide methods for the following generics:
 *
 * - {@linkcode NUMBER_OF_ROWS}
 * - {@linkcode SLICE_2D}
 * - {@linkcode COMBINE_ROWS}
 * - {@linkcode CLONE}
 *
 * Each alternative experiment should be a {@linkplain SummarizedExperiment} with number of columns equal to that of the SingleCellExperiment.
 *
 * @extends RangedSummarizedExperiment
 */
export class SingleCellExperiment extends rse.RangedSummarizedExperiment {
    static #check_reduced_dimensions(reducedDimensions, numberOfColumns) {
        for (const [k, v] of Object.entries(reducedDimensions)) {
            if (generics.NUMBER_OF_ROWS(v) !== numberOfColumns) {
                throw new Error("number of rows for reduced dimension '" + k + "' is not equal to number of columns for this SingleCellExperiment");
            }
        }
    }

    static #check_alternative_experiments(alternativeExperiments, numberOfColumns) {
        for (const [k, v] of Object.entries(alternativeExperiments)) {
            if (!(v instanceof se.SummarizedExperiment)) {
                throw new Error("alternative experiment '" + k + "' is not a SummarizedExperiment");
            }
            if (v.numberOfColumns(v) !== numberOfColumns) {
                throw new Error("number of columns for alternative experiment '" + k + "' is not equal to number of columns for this SingleCellExperiment");
            }
        }
    }

    /**
     * @param {Object} assays - Object where keys are the assay names and values are multi-dimensional arrays of experimental data.
     * @param {Object} [options={}] - Optional parameters.
     * @param {Object} [reducedDimensions={}] - Object containing named reduced dimensions.
     * Each value should be a 2-dimensional object with number of rows equal to the number of columns of the assays.
     * @param {?Array} [reducedDimensionOrder=null] - Array containing the order of the reduced dimensions.
     * This should have the same values as the keys of `reducedDimensions`, and defaults to those keys if `null`.
     * @param {Object} [alternativeExperiments={}] - Object containing named alternative experiments.
     * Each value should be a 2-dimensional object with number of columns equal to that of the assays.
     * @param {?Array} [alternativeExperimentOrder=null] - Array containing the order of the alternative experiments.
     * This should have the same values as the keys of `alternativeExperiments`, and defaults to those keys if `null`.
     * @param {?(GRanges|GroupedGRanges)} [rowRanges=null] - Genomic ranges corresponding to each row, see the {@linkplain RangedSummarizedExperiment} constructor.
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
        reducedDimensions = {},
        reducedDimensionOrder = null,
        alternativeExperiments = {},
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

        SingleCellExperiment.#check_reduced_dimensions(reducedDimensions, this.numberOfColumns());
        this._reducedDimensions = {
            entries: reducedDimensions,
            order: utils.checkEntryOrder(reducedDimensions, reducedDimensionOrder, "reducedDimension")
        };

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

    /**
     * @return {Array} Array of strings containing the (ordered) names of the reduced dimensions.
     */
    reducedDimensionNames() {
        return this._reducedDimensions.order;
    }

    /**
     * @param {string|number} i - Reduced dimension to retrieve, either by name or index.
     * @return {*} The contents of reduced dimension `i` as an multi-dimensional array-like object.
     */
    reducedDimension(i) {
        return utils.retrieveSingleEntry(this._reducedDimensions.entries, this._reducedDimensions.order, i, "reduced dimension", "SingleCellExperiment");
    }

    /**
     * @return {Array} Array of strings containing the (ordered) names of the alternative experiments.
     */
    alternativeExperimentNames() {
        return this._alternativeExperiments.order;
    }

    /**
     * @param {string|number} i - Alternative experiment to retrieve, either by name or index.
     * @return {SummarizedExperiment} The specified alternative experiment `i`. 
     */
    alternativeExperiment(i) {
        return utils.retrieveSingleEntry(this._alternativeExperiments.entries, this._alternativeExperiments.order, i, "alternative experiment", "SingleCellExperiment");
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {string|number} i - Identity of the reduced dimension to remove, either by name or index.
     * @return {SingleCellExperiment} Reference to this SingleCellExperiment after removing the specified assay.
     */
    $removeReducedDimension(i) {
        utils.removeSingleEntry(this._reducedDimensions.entries, this._reducedDimensions.order, i, "reduced dimension", "SingleCellExperiment");
        return this;
    }

    /**
     * @param {string|number} i - Identity of the reduced dimension to add, either by name or index.
     * - If `i` is a number, the reduced dimension at the specified index is replaced.
     *   `i` should be non-negative and less than the number of reduced dimensions.
     * - If `i` is a string, any reduced dimension with the same name is replaced.
     *   If no such reduced dimension exists, a new reduced dimension is appended to the list of reduced dimensions.
     * @param {*} value - Multi-dimensional array-like object to set/add as the reduced dimension.
     * @return {SingleCellExperiment} Reference to this SingleCellExperiment with modified reduced dimensions.
     */
    $setReducedDimension(i, value) {
        if (generics.NUMBER_OF_ROWS(value) != this.numberOfColumns()) {
            throw new Error("number of rows of 'value' should be the same as the number of columns of this SingleCellExperiment");
        }
        utils.setSingleEntry(this._reducedDimensions.entries, this._reducedDimensions.order, i, value, "reduced dimension", "SingleCellExperiment");
        return this;
    }

    /**
     * @param {string|number} i - Identity of the reduced dimension to remove, either by name or index.
     * @return {SingleCellExperiment} Reference to this SingleCellExperiment after removing the specified assay.
     */
    $removeAlternativeExperiment(i) {
        utils.removeSingleEntry(this._alternativeExperiments.entries, this._alternativeExperiments.order, i, "alternative experiment", "SingleCellExperiment");
        return this;
    }

    /**
     * @param {string|number} i - Identity of the alternative experiment to add, either by name or index.
     * - If `i` is a number, the alternative experiment at the specified index is replaced.
     *   `i` should be non-negative and less than the number of alternative experiments.
     * - If `i` is a string, any alternative experiment with the same name is replaced.
     *   If no such alternative experiment exists, a new alternative experiment is appended to the list of alternative experiments.
     * @param {*} value - Multi-dimensional array-like object to set/add as the alternative experiment.
     * @return {SingleCellExperiment} Reference to this SingleCellExperiment with modified alternative experiments.
     */
    $setAlternativeExperiment(i, value) {
        if (!(value instanceof se.SummarizedExperiment) || generics.NUMBER_OF_COLUMNS(value) != this.numberOfColumns()) {
            throw new Error("'value' should be a SummarizedExperiment with the same number of columns as this SingleCellExperiment");
        }
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
