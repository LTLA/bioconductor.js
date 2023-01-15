import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";
import * as df from "./DataFrame.js";
import * as ir from "./IRanges.js";

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
 */
export class GRanges {
    /**
     * @param {Array} seqnames - Array of strings containing the sequence names for each genomic range.
     * @param {IRanges} ranges - Position and width of the range on its specified sequence.
     * This should have the same length as `seqnames`.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?Array} [options.names=null] - Array of strings of length equal to `start`, containing names for each genomic range.
     * Alternatively `null`, in which case the ranges are assumed to be unnamed.
     * @param {?DataFrame} [options.rangeMetadata=null] - A {@linkplain DataFrame} with number of rows equal to the length of `start`, containing arbitrary per-range annotations.
     * Alternatively `null`, in which case a zero-column DataFrame is automatically constructed.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(seqnames, ranges, { names = null, rangeMetadata = null, metadata = {} } = {}) {
        utils.checkStringArray(seqnames, "seqnames");
        this._seqnames = seqnames;

        let n = seqnames.length;
        if (n !== generics.LENGTH(ranges)) {
            throw utils.formatLengthError("'ranges'", "'seqnames'");
        }
        this._ranges = ranges;

        this._rangeMetadata = ir.verifyRangeMetadata(rangeMetadata, n, "seqnames.length");

        if (names !== null) {
            utils.checkNamesArray(names, "'names'", n, "'seqnames.length'");
        }
        this._names = names;

        // Storing this separately.
        this._metadata = metadata;
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

    /**
     * @return {?Array} Array of strings containing the name of each genomic range, or `null` if no names are available.
     */
    names() {
        return this._names;
    }

    /**
     * @return {DataFrame} A {@linkplain DataFrame} with one row corresponding to each genomic range, containing arbitrary per-range metadata.
     */
    rangeMetadata() {
        return this._rangeMetadata;
    }

    /**
     * @return {Object} Object containing arbitrary metadata.
     */
    metadata() {
        return this._metadata;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {Array} seqnames - Array of strings containing the sequence names for each genomic range.
     * @return {GRanges} A reference to this GRanges object, after setting the sequence names to `seqnames`.
     */
    $setSeqnames(seqnames) {
        utils.checkNamesArray(value, "replacement 'seqnames'", generics.LENGTH(this), "'LENGTH(<GRanges>)'");
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
        if (generics.LENGTH(ranges) !== generics.LENGTH(this._ranges)) {
            throw utils.formatLengthError("replacement 'ranges'", "'LENGTH(<GRanges>)'");
        }
        this._ranges = ranges;
        return this;
    }

    /**
     * @param {?Array} names - Array of strings containing a name for each genomic range.
     * This should have length equal to the number of ranges.
     * Alternatively `null`, if no names are present.
     * @return {GRanges} A reference to this GRanges object, after setting the names to `names`.
     */
    $setNames(names) {
        if (names !== null) {
            utils.checkNamesArray(names, "replacement 'names'", generics.LENGTH(this), "'LENGTH(<GRanges>)'");
        }
        this._names = names;
        return this;
    }

    /**
     * @param {?DataFrame} value - Arbitrary metadata for each genomic range.
     * This should have number of rows equal to the number of ranges.
     * Alternatively `null`, in which case all existing per-range metadata is removed.
     * @return {GRanges} A reference to this GRanges object, after setting the widths to `value`.
     */
    $setRangeMetadata(value) {
        this._rangeMetadata = ir.verifyRangeMetadata(rangeMetadata, n, "seqnames.length");
        return this;
    }

    /**
     * @param {Object} value - Object containing the metadata.
     *
     * @return {IRanges} Reference to this DataFrame after replacing the metadata.
     */
    $setMetadata(value) {
        this._metadata = value;
        return this;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_LENGTH() {
        return this._start.length;
    }

    _bioconductor_SLICE(i, { allowInPlace = false, allowView = false }) {
        let options = { allowInPlace, allowView };
        let s = generics.SLICE(this._start, i, options);
        let w = generics.SLICE(this._width, i, options);
        let r = generics.SLICE(this._rangeMetadata, i, options);

        if (allowInPlace) {
            this._start = s;
            this._width = w;
            this._rangeMetadata = r;
            return this;
        } else {
            let output = Object.create(this.constructor.prototype); // avoid validity checks.
            output._start = s;
            output._width = w;
            output._rangeMetadata = r;
            output._metadata = this._metadata;
            return output;
        }
    }

    _bioconductor_COMBINE(objects, { allowAppend = false }) {
        let all_s = [];
        let all_w = [];
        let all_r = [];

        for (const x of objects) {
            all_s.push(x.start());
            all_w.push(x.width());
            all_r.push(x.rangeMetadata());
        }

        let combined_s = generics.COMBINE(all_s);
        let combined_w = generics.COMBINE(all_w);
        let combined_r = generics.COMBINE(all_r);

        if (allowAppend) {
            this._start = combined_s;
            this._width = combined_w;
            this._rangeMetadata = combined_r;
            return this;
        } else {
            let output = Object.create(this.constructor.prototype); // avoid validity checks.
            output._start = combined_s;
            output._width = combined_w;
            output._rangeMetadata = combined_r;
            output._metadata = this._metadata;
            return output;
        }
    }

    _bioconductor_CLONE({ deepcopy = true }) {
        let options = { deepcopy };
        let output = Object.create(this.constructor.prototype); // avoid validity checks.
        output._start = generics.CLONE(this._start, options);
        output._width = generics.CLONE(this._width, options);
        output._rangeMetadata = generics.CLONE(this._rangeMetadata, options);
        output._metadata = generics.CLONE(this._metadata, options);
        return output;
    }
}

