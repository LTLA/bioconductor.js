import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";
import * as ir from "./IRanges.js";
import * as vec from "./Vector.js";

/**
 * A GRanges object is a collection of genomic ranges, inspired by the class of the same name from the Bioconductor ecosystem.
 * Each range consists of a sequence name, a start position on that sequence, and a width.
 * Each range may also be associated with arbitrary range-level metadata in a {@linkplain DataFrame}.
 * The GRanges defines methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 *
 * @extends Vector
 */
export class GRanges extends vec.Vector {
    /**
     * @param {Array} seqnames - Array of strings containing the sequence names for each genomic range.
     * @param {IRanges} ranges - Position and width of the range on its specified sequence.
     * This should have the same length as `seqnames`.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?Array} [options.names=null] - Array of strings of length equal to `start`, containing names for each genomic range.
     * Alternatively `null`, in which case the ranges are assumed to be unnamed.
     * @param {?DataFrame} [options.elementMetadata=null] - A {@linkplain DataFrame} with number of rows equal to the length of `start`, containing arbitrary per-range annotations.
     * Alternatively `null`, in which case a zero-column DataFrame is automatically constructed.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(seqnames, ranges, { names = null, elementMetadata = null, metadata = {} } = {}) {
        super(seqnames.length, { names, elementMetadata, metadata });

        utils.checkStringArray(seqnames, "seqnames");
        this._seqnames = seqnames;

        let n = seqnames.length;
        if (n !== generics.LENGTH(ranges)) {
            throw utils.formatLengthError("'ranges'", "'seqnames'");
        }
        this._ranges = ranges;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {Int32Array} Array of integers containing the start position for each genomic range.
     */
    start() {
        return this._ranges.start();
    }

    /**
     * @return {Int32Array} Array of integers containing the end position (specifically, one-past-the-end) for each genomic range.
     */
    end() {
        return this._ranges.end();
    }

    /**
     * @return {Int32Array} Array of integers containing the width of each genomic range.
     */
    width() {
        return this._ranges.width();
    }

    /**
     * @return {Array} Array of strings containing the sequence name for each genomic range.
     */
    seqnames() {
        return this._seqnames;
    }

    /**
     * @return {IRanges} Start positions and widths for all ranges on their specified sequence names.
     */
    ranges() {
        return this._ranges;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {Array} seqnames - Array of strings containing the sequence names for each genomic range.
     * @return {GRanges} A reference to this GRanges object, after setting the sequence names to `seqnames`.
     */
    $setSeqnames(seqnames) {
        utils.checkNamesArray(seqnames, "replacement 'seqnames'", generics.LENGTH(this), "'LENGTH(<GRanges>)'");
        this._seqnames = seqnames;
        return this;
    }

    /**
     * @param {Array|TypedArray} start - Array of start positions for each genomic range.
     * This should have length equal to the number of ranges and be coercible into an Int32Array.
     * @return {GRanges} A reference to this GRanges object, after setting the start positions to `start`.
     */
    $setStart(start) {
        this._ranges.$setStart(start);
        return this;
    }

    /**
     * @param {Array|TypedArray} width - Array of widths for each genomic range.
     * This should have length equal to the number of ranges and be coercible into an Int32Array.
     * @return {GRanges} A reference to this GRanges object, after setting the widths to `width`.
     */
    $setWidth(width) {
        this._ranges.$setWidth(width);
        return this;
    }

    /**
     * @param {IRanges} ranges - Start positions and widths for each genomic range.
     * This should have length equal to the number of ranges. 
     * @return {GRanges} A reference to this GRanges object, after setting the ranges to `ranges`.
     */
    $setRanges(ranges) {
        if (!(ranges instanceof ir.IRanges)) {
            throw new Error("'ranges' should be an IRanges object");
        }
        if (generics.LENGTH(ranges) !== generics.LENGTH(this._ranges)) {
            throw utils.formatLengthError("replacement 'ranges'", "'LENGTH(<GRanges>)'");
        }
        this._ranges = ranges;
        return this;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_LENGTH() {
        return this._seqnames.length;
    }

    _bioconductor_SLICE(output, i, { allowView = false }) {
        super._bioconductor_SLICE(output, i, { allowView });
        output._seqnames = generics.SLICE(this._seqnames, i, { allowView });
        output._ranges = generics.SLICE(this._ranges, i, { allowView });
        return;
    }

    _bioconductor_COMBINE(output, objects) {
        super._bioconductor_COMBINE(output, objects);

        let all_s = [];
        let all_rr = [];
        for (const x of objects) {
            all_s.push(x._seqnames);
            all_rr.push(x._ranges);
        }

        output._seqnames = generics.COMBINE(all_s);
        output._ranges = generics.COMBINE(all_rr);
        return;
    }

    _bioconductor_CLONE(output, { deepcopy = true }) {
        let options = { deepcopy };
        super._bioconductor_CLONE(output, options);
        output._seqnames = generics.CLONE(this._seqnames, options);
        output._ranges = generics.CLONE(this._ranges, options);
        return;
    }
}

