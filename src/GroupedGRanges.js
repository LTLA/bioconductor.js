import * as vec from "./Vector.js";
import * as gr from "./GRanges.js";
import * as utils from "./utils.js";
import * as generics from "./AllGenerics.js";

/**
 * A GroupedGRanges object is a collection of groups of genomic ranges, inspired by the `GRangesList` class from the Bioconductor ecosystem.
 * Each group consists of a {@linkplain GRanges} object of arbitrary length, which is most often used to represent a multi-exon gene.
 * The GroupedGRanges can be considered a vector of groups, and defines methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 *
 * @extends Vector
 */
export class GroupedGRanges extends vec.Vector {
    static #computeStarts(lengths) {
        let starts = new Int32Array(lengths.length);
        let last = 0;
        for (var i = 0; i < lengths.length; i++) {
            starts[i] = last;
            last += lengths[i];
        }
        return { starts: starts, total: last };
    }

    /**
     * @param {Array|GRanges} ranges - An array of {@linkplain GRanges} objects, where each element represents a group of genomic ranges.
     * 
     * Alternatively, a single GRanges containing a concatenation of ranges from all groups.
     * In this case, `rangeLengths` must be supplied.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?(TypedArray|Array)} [rangeLengths=null] - Length of the ranges within each group.
     * This should be coercible to an Int32Array, contain non-negative values, and have a sum equal to the length of `ranges`.
     * Only used if `ranges` is a single {@linkplain GRanges} object, where each group's ranges are assumed to form contiguous intervals along `ranges`.
     * @param {?Array} [options.names=null] - Array of strings of length equal to `start`, containing names for each genomic range.
     * Alternatively `null`, in which case the ranges are assumed to be unnamed.
     * @param {?DataFrame} [options.elementMetadata=null] - A {@linkplain DataFrame} with number of rows equal to the length of `start`, containing arbitrary per-range annotations.
     * Alternatively `null`, in which case a zero-column DataFrame is automatically constructed.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(ranges, { rangeLengths = null, names = null, elementMetadata = null, metadata = {} } = {}) {
        if (ranges.constructor == Array) {
            super(ranges.length, { names, elementMetadata, metadata });
            rangeLengths = new Int32Array(ranges.length);
            for (var i = 0; i < rangeLengths.length; i++) {
                if (!(ranges[i] instanceof gr.GRanges)) {
                    throw new Error("'ranges' must either be a 'GRanges' or an array of 'GRanges'");
                }
                rangeLengths[i] = generics.LENGTH(ranges[i]);
            }
            ranges = generics.COMBINE(ranges);

        } else {
            if (!(ranges instanceof gr.GRanges)) {
                throw new Error("'ranges' must either be a 'GRanges' or an array of 'GRanges'");
            }
            if (rangeLengths == null) {
                throw new Error("'rangeLengths' must be specified when 'ranges' is a 'GRanges'");
            }
            super(rangeLengths.length, { names, elementMetadata, metadata });
            rangeLengths = utils.convertToInt32Array(rangeLengths);
            utils.checkNonNegative(rangeLengths);
        }

        this._ranges = ranges;
        this._rangeLengths = rangeLengths;

        let accumulated = GroupedGRanges.#computeStarts(rangeLengths);
        this._rangeStarts = accumulated.starts;

        if (accumulated.total !== generics.LENGTH(ranges)) {
            throw new Error("sum of 'rangeLengths' must be equal to the length of 'ranges'");
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {GRanges} The concatenated set of ranges across all groups. 
     */
    ranges() {
        return this._ranges;
    }

    /**
     * @return {Int32Array} The start indices for each group's ranges along the concatenated set of ranges returned by {@linkcode GroupedGRanges#ranges ranges}.
     */
    rangeStarts() {
        return this._rangeStarts;
    }

    /**
     * @return {Int32Array} The length of each group's ranges along the concatenated set of ranges returned by {@linkcode GroupedGRanges#ranges ranges}.
     */
    rangeLengths() {
        return this._rangeLengths;
    }

    /**
     * @param {number} i - Index of the group of interest.
     * @param {boolean} [options.allowView=false] - Whether a view can be created in any internal slicing operations.
     *
     * @return {GRanges} The genomic ranges for group `i`.
     */
    group(i, { allowView = false } = {}) {
        let s = this._rangeStarts[i];
        return generics.SLICE(this._ranges, { start: s, end: s + this._rangeLengths[i] }, { allowView });
    }

    /**
     * @return {number} Number of groups in this object.
     */
    numberOfGroups() {
        return this._rangeStarts.length;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {GRanges} ranges - Genomic ranges of length equal to the concatenated set of ranges returned by {@linkcode GroupedGRanges#ranges ranges}.
     * @return {GroupedGRanges} A reference to this GroupedGRanges object after modifying the internal ranges.
     */
    $setRanges(ranges) {
        if (!(ranges instanceof gr.GRanges)) {
            throw new Error("'ranges' must be a 'GRanges'");
        }
        if (generics.LENGTH(ranges) !== generics.LENGTH(this._ranges)) {
            throw utils.formatLengthError("'ranges'", "number of ranges");
        }
        this._ranges = ranges;
        return this;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {Object} [options={}] - Optional parameters.
     * @param {?(Array|Set)} [options.restrictToSeqnames=null] - Array or Set containing the sequence names to use in the index.
     * If `null`, all available sequence names are used.
     * @param {?(Array|Set)} [options.restrictToStrand=null] - Array or Set containing the strands to use in the index.
     * If `null`, all available strands are used.
     *
     * @return {GroupedGRangesOverlapIndex} A pre-built index for computing overlaps with other {@linkplain GRanges} instances.
     */
    buildOverlapIndex({ restrictToSeqnames = null, restrictToStrand = null } = {}) {
        return new GroupedGRangesOverlapIndex(
            this._ranges.buildOverlapIndex({ restrictToSeqnames, restrictToStrand }),
            this._rangeStarts,
            this._rangeLengths
        );
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_LENGTH() {
        return this._rangeStarts.length;
    }

    _bioconductor_SLICE(output, i, { allowView = false } = {}) {
        super._bioconductor_SLICE(output, i, { allowView });

        output._rangeLengths = generics.SLICE(this._rangeLengths, i, { allowView });
        let accumulated = GroupedGRanges.#computeStarts(output._rangeLengths);
        output._rangeStarts = accumulated.starts;

        if (i.constructor == Object) {
            // Handle this specially for optimizing allowView = true.
            let s = this._rangeStarts[i.start];
            output._ranges = generics.SLICE(this._ranges, { start: s, end: s + accumulated.total }, { allowView });
        } else {
            let keep = new Int32Array(accumulated.total);

            let counter = 0;
            for (const j of i) {
                let start = this._rangeStarts[j];
                let end = start + this._rangeLengths[j];
                for (var k = start; k < end; k++) {
                    keep[counter] = k;
                    counter++;
                }
            }

            output._ranges = generics.SLICE(this._ranges, keep, { allowView });
        }
        return;
    }

    _bioconductor_COMBINE(output, objects) {
        super._bioconductor_COMBINE(output, objects);

        output._rangeLengths = generics.COMBINE(objects.map(x => x._rangeLengths));
        let accumulated = GroupedGRanges.#computeStarts(output._rangeLengths);
        output._rangeStarts = accumulated.starts;
        output._ranges = generics.COMBINE(objects.map(x => x._ranges));

        return;
    }

    _bioconductor_CLONE(output, { deepCopy = true }) {
        super._bioconductor_CLONE(output, { deepCopy });
        output._rangeLengths = generics.CLONE(this._rangeLengths, { deepCopy });
        output._rangeStarts = generics.CLONE(this._rangeStarts, { deepCopy });
        output._ranges = generics.CLONE(this._ranges, { deepCopy });
        return;
    }
}

/**
 * Pre-built index for overlapping {@linkplain GroupedGRanges} objects.
 * This is typically constructed using the {@linkcode GroupedGRanges#buildOverlapIndex GroupedGRanges.buildOverlapIndex} method for a "reference" object,
 * and can be applied to different query GroupedGRanges or {@linkplain GRanges} to identify overlaps with the reference.
 *
 * @hideconstructor
 */
export class GroupedGRangesOverlapIndex {
    constructor(index, rangeStarts, rangeLengths) {
        this._index = index;
        this._rangeStarts = rangeStarts;
        this._rangeLengths = rangeLengths;
    }

    /**
     * @param {GroupedGRanges|GRanges} query - The query object, containing ranges to be overlapped with those in the reference GroupedGRanges (that was used to construct this GroupedGRangesOverlapIndex object).
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.ignoreStrand=true] - Whether to ignore differences in strandedness between the ranges in `query` and the reference object.
     *
     * @return {Array} An array of length equal to the number of ranges or groups in `query`,
     * where each element is an array containing the indices of the overlapping ranges in the reference {@linkplain GRanges} object.
     */
    overlap(query, { ignoreStrand = true } = {}) {
        let output = new Array(this._rangeStarts.length);

        if (query instanceof GroupedGRanges) {
            let overlaps = this._index.overlap(query._ranges);

            let rev_map = new Int32Array(query._ranges.length);
            for (var i = 0; i < query._rangeStarts; i++) {
                let start = query._rangeStarts[i];
                let end = start + query._rangeLengths[i];
                for (var s = start; s < end; s++) {
                    rev_map[s] = i;
                }
            }

            for (var i = 0; i < this._rangeStarts; i++) {
                let start = this._rangeStarts[i];
                let end = start + this._rangeLengths[i];

                let results = new Set;
                for (var s = start; s < end; s++) {
                    overlaps[s].forEach(x => results.add(rev_map[x]));
                }
                output[i] = Array.from(results);
            }

        } else {
            let overlaps = this._index.overlap(query);

            for (var i = 0; i < this._rangeStarts; i++) {
                let start = this._rangeStarts[i];
                let end = start + this._rangeLengths[i];

                let results = new Set;
                for (var s = start; s < end; s++) {
                    overlaps[s].forEach(x => results.add(x));
                }
                output[i] = Array.from(results);
            }
        }

        return output;
    }

}
