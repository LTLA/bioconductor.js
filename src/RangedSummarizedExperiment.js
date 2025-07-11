import * as generics from "./AllGenerics.js";
import * as gr from "./GRanges.js";
import * as ggr from "./GroupedGRanges.js";
import * as se from "./SummarizedExperiment.js";
import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";

/**
 * A RangedSummarizedExperiment is a {@linkplain SummarizedExperiment} subclass where each row represents a genomic interval.
 * As such, it stores an additional {@linkplain GRanges} or {@linkplain GroupedGRanges} of length equal to the number of rows,
 * where each element represents the genomic range(s) for the corresponding row of the SummarizedExperiment.
 *
 * The RangedSummarizedExperiment supports the same set of generics as the {@linkplain SummarizedExperiment}.
 * Each method will call the base method, with the following extensions:
 *
 * - {@linkcode SLICE_2D} will additionally slice the supplied genomic ranges by the desired `rows`.
 * - {@linkcode COMBINE_ROWS} will combine genomic ranges across objects.
 *   If some objects contain a GroupedGRanges and other objects contain GRanges, the latter will be coerced to a GroupedGRanges (where each group contains one range) before combining.
 *   If any object is a base SummarizedExperiment, a GroupedGRanges containing zero-length groups will be automatically constructed to attempt combining.
 * - {@linkcode COMBINE_COLUMNS} will use the genomic ranges from the first object.
 *
 * Constructors of RangedSummarizedExperiment subclasses should be callable with no arguments, possibly creating an empty object with no properties.
 * This will be used by the `_bioconductor_CLONE`, `_bioconductor_COMBINE_ROWS`, `_bioconductor_COMBINE_COLUMNS` and `_bioconductor_SLICE_2D` methods to return an instance of the subclass.
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

    $setRowRanges(value) {
        return this.setRowRanges(value, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_SLICE_2D(rows, columns, { allowView = false }) {
        let output = super._bioconductor_SLICE_2D(rows, columns, { allowView });
        if (rows !== null) {
            output._rowRanges = generics.SLICE(this._rowRanges, rows);
        } else {
            output._rowRanges = this._rowRanges;
        }
        return output;
    }

    _bioconductor_COMBINE_ROWS(objects) {
        let output = super._bioconductor_COMBINE_ROWS(objects);

        let collected = [this._rowRanges];
        let has_empty = false;
        let has_ggr = (this.rowRanges instanceof ggr.GroupedGRanges);
        for (const x of objects) {
            if (x instanceof RangedSummarizedExperiment) {
                let y = x._rowRanges;
                if (y instanceof ggr.GroupedGRanges) {
                    has_ggr = true;
                }
                collected.push(y);
            } else if (x instanceof se.SummarizedExperiment) {
                has_empty = true;
                collected.push(null);
            } else {
                throw new Error("objects to be combined must be SummarizedExperiments (failing for object " + String(i) + ")");
            }
        }

        // Promoting nulls and GRanges to GroupedGRanges, if necessary.
        if (has_empty || has_ggr) {
            for (var i = 0; i < collected.length; i++) {
                let current = collected[i];

                if (current instanceof gr.GRanges) {
                    let widths = new Int32Array(generics.LENGTH(current));
                    widths.fill(1);

                    let options = { 
                        rangeLengths: widths,
                        names: current.names(),
                        elementMetadata: current.elementMetadata(),
                        metadata: current.metadata()
                    };

                    if (options.names !== null) {
                        current = current.setNames(null);
                    } 

                    if (options.elementMetadata.metadata().size > 0 || options.elementMetadata.numberOfColumns() > 0) {
                        current = current.setElementMetadata(null);
                    }

                    if (options.metadata.size > 0) {
                        current = current.setMetadata(new Map);
                    }

                    collected[i] = new ggr.GroupedGRanges(current, options);

                } else if (current === null){
                    const x = (i == 0 ? this : objects[i - 1]);
                    collected[i] = ggr.GroupedGRanges.empty(x.numberOfRows());
                }
            }
        }

        output._rowRanges = generics.COMBINE(collected);
        return output;
    }

    _bioconductor_COMBINE_COLUMNS(objects) {
        let output = super._bioconductor_COMBINE_COLUMNS(objects);
        output._rowRanges = this._rowRanges;
        return output;
    }

    _bioconductor_CLONE({ deepCopy }) {
        let output = super._bioconductor_CLONE({ deepCopy });
        output._rowRanges = cutils.cloneField(this._rowRanges, deepCopy);
        return output;
    }
}
