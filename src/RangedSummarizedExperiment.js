import * as generics from "./AllGenerics.js";
import * as gr from "./GRanges.js";
import * as ggr from "./GroupedGRanges.js";
import * as se from "./SummarizedExperiment.js";
import * as utils from "./utils.js";

/**
 * A RangedSummarizedExperiment is a {@linkplain SummarizedExperiment} subclass where each row represents a genomic interval.
 * As such, it stores an additional {@linkplain GRanges} of length equal to the number of rows,
 * where each element represents the genomic range for the corresponding row of the SummarizedExperiment.
 * It supports the same set of generics as the {@linkplain SummarizedExperiment}.
 *
 * @extends SummarizedExperiment
 */
export class RangedSummarizedExperiment extends se.SummarizedExperiment {
    #check_rowRanges(x) {
        if (!(x instanceof gr.GRanges) && !(x instanceof ggr.GroupedGRanges)) {
            throw new Error("'rowRanges' should be a 'GRanges' or 'GroupedGRanges' instance");
        }
        if (generics.LENGTH(x) !== this._rowData.numberOfRows()) {
            throw utils.formatLengthError("'rowRanges'", "the number of rows");
        }
    }

    /**
     * @param {Object} assays - Object where keys are the assay names and values are multi-dimensional arrays of experimental data.
     * All arrays should have the same number of rows and columns.
     * @param {?(GRanges|GroupedGRanges)} rowRanges - Genomic ranges corresponding to each row.
     *
     * Alternatively, each row may correspond to a group of genomic ranges.
     *
     * If `null`, a {@linkplain GroupedGRanges} is constructed where each row corresponds to one group of ranges of zero length.
     * @param {Object} [options={}] - Optional parameters.
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
    constructor(assays, rowRanges, { assayOrder = null, rowData = null, columnData = null, rowNames = null, columnNames = null, metadata = {} } = {}) {
        if (arguments.length == 0) {
            super();
            return;
        }

        super(assays, { assayOrder, rowData, columnData, rowNames, columnNames, metadata });

        if (rowRanges === null) {
            rowRanges = ggr.GroupedGRanges.empty(this.numberOfRows());
        } else {
            this.#check_rowRanges(rowRanges);
        }
        this._rowRanges = rowRanges;

        return;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {GRanges} Genomic ranges corresponding to each row.
     */
    rowRanges() {
        return this._rowRanges;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {GRanges} value - Genomic ranges corresponding to each row.
     * This should have length equal to the number of rows in this RangedSummarizedExperiment.
     * @return {RangedSummarizedExperiment} A reference to this RangedSummarizedExperiment after modifying its `rowRanges`.
     */
    $setRowRanges(value) {
        this.#check_rowRanges(value);
        this._rowRanges = value;
        return this;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_SLICE_2D(output, rows, columns, { allowView = false }) {
        super._bioconductor_SLICE_2D(output, rows, columns, { allowView });
        if (rows !== null) {
            output._rowRanges = generics.SLICE(this._rowRanges, rows);
        } else {
            output._rowRanges = this._rowRanges;
        }
    }

    _bioconductor_COMBINE_ROWS(output, objects) {
        super._bioconductor_COMBINE_ROWS(output, objects);

        let collected = objects.map(x => x._rowRanges);
        output._rowRanges = generics.COMBINE(collected);

        return;
    }

    _bioconductor_COMBINE_COLUMNS(output, objects) {
        super._bioconductor_COMBINE_COLUMNS(output, objects);

        output._rowRanges = objects[0]._rowRanges;

        return;
    }

    _bioconductor_CLONE(output, { deepCopy }) {
        super._bioconductor_CLONE(output, { deepCopy });

        output._rowRanges = generics.CLONE(this._rowRanges, { deepCopy });

        return;
    }
}
