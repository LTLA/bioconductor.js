import * as generics from "./AllGenerics.js";
import * as gr from "./GRanges.js";
import * as ggr from "./GroupedGRanges.js";
import * as se from "./SummarizedExperiment.js";
import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";

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
     * @param {Object} [options={}] - Optional parameters, including those used in the {@linkplain SummarizedExperiment} constructor.
     */
    constructor(assays, rowRanges, options = {}) {
        if (arguments.length == 0) {
            super();
            return;
        }

        super(assays, options);

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
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this Annotated instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {RangedSummarizedExperiment} The RangedSummarizedExperiment after modifying its `rowRanges`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setRowRanges(value, { inPlace = false } = {}) {
        this.#check_rowRanges(value);
        let target = cutils.setterTarget(this, inPlace);
        target._rowRanges = value;
        return target;
    }

    /**
     * @param {GRanges} value - Genomic ranges corresponding to each row.
     * This should have length equal to the number of rows in this RangedSummarizedExperiment.
     * @return {RangedSummarizedExperiment} A reference to this RangedSummarizedExperiment after modifying its `rowRanges`.
     */
    $setRowRanges(value) {
        return this.setRowRanges(value, { inPlace: true });
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

        output._rowRanges = cutils.cloneField(this._rowRanges, deepCopy);

        return;
    }
}
