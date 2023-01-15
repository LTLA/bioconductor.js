import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";
import * as df from "./DataFrame.js";

/**
 * An IRanges object is a collection of integer ranges, inspired by the class of the same name from the Bioconductor ecosystem.
 * Each range consists of a start position and a width, and may be associated with arbitrary range-level metadata in a {@linkplain DataFrame}.
 * The IRanges defines methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 */
export class IRanges {
    static #convertToInt32Array(x) {
        if (x instanceof Int32Array) {
            return x;
        } else {
            return new Int32Array(x);
        }
    }

    static #checkNonNegative(x, msg) {
        for (const y of x) {
            if (y < 0) {
                throw new Error("detected a negative entry in '" + msg + "'");
            }
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {Array|TypedArray} start - Array of start positions for each range.
     * This should be coercible into an Int32Array.
     * @param {Array|TypedArray} width - Array of widths for each range.
     * This should be coercible into an Int32Array.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?Array} [options.names=null] - Array of strings of length equal to `start`, containing names for each range.
     * Alternatively `null`, in which case the ranges are assumed to be unnamed.
     * @param {?DataFrame} [options.rangeMetadata=null] - A {@linkplain DataFrame} with number of rows equal to the length of `start`, containing arbitrary per-range annotations.
     * Alternatively `null`, in which case a zero-column DataFrame is automatically constructed.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(start, width, { names = null, rangeMetadata = null, metadata = {} } = {}) {
        this._start = IRanges.#convertToInt32Array(start);
        IRanges.#checkNonNegative(this._start, "start");

        this._width = IRanges.#convertToInt32Array(width);
        IRanges.#checkNonNegative(this._width, "width");

        let n = this._start.length;
        if (n !== this._width.length) {
            throw new Error("'start' and 'width' should have the same length");
        }

        if (rangeMetadata !== null) {
            if (!(rangeMetadata instanceof df.DataFrame)) {
                throw new Error("'rangeMetadata' should be a DataFrame");
            }
            if (generics.LENGTH(rangeMetadata) !== n) {
                throw new Error("'rangeMetadata' should have the same number of rows as 'start.length'");
            }
        } else {
            rangeMetadata = new df.DataFrame({}, { numberOfRows: n });
        }
        this._rangeMetadata = rangeMetadata;

        if (names !== null) {
            utils.checkNamesArray(names);
            if (names.length != n) {
                throw new Error("'start' and 'names' should have the same length");
            }
        }
        this._rangeMetadata.$setRowNames(names);

        this._metadata = metadata;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {Int32Array} Array of integers containing the start position for each range.
     */
    start() {
        return this._start;
    }

    /**
     * @return {Int32Array} Array of integers containing the end position (specifically, one-past-the-end) for each range.
     */
    end() {
        return this._start.map((x, i) => x + this._width[i]);
    }

    /**
     * @return {Int32Array} Array of integers containing the width of each range.
     */
    width() {
        return this._width;
    }

    /**
     * @return {?Array} Array of strings containing the name of each range, or `null` if no names are available.
     */
    names() {
        return this._rangeMetadata.rowNames();
    }

    /**
     * @return {DataFrame} A {@linkplain DataFrame} with one row corresponding to each range, containing arbitrary per-range metadata.
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
     * @param {Array|TypedArray} value - Array of start positions for each range.
     * This should have length equal to the number of ranges and be coercible into an Int32Array.
     * @return {IRanges} A reference to this IRanges object, after setting the start positions to `value`.
     */
    $setStart(value) {
        let candidate = IRanges.#convertToInt32Array(value);
        if (candidate.length !== generics.LENGTH(this)) {
            throw new Error("'start' should be replaced by array of the same length");
        }
        IRanges.#checkNonNegative(candidate, "start");
        this._start = candidate;
        return this;
    }

    /**
     * @param {Array|TypedArray} value - Array of widths for each range.
     * This should have length equal to the number of ranges and be coercible into an Int32Array.
     * @return {IRanges} A reference to this IRanges object, after setting the widths to `value`.
     */
    $setWidth(value) {
        let candidate = IRanges.#convertToInt32Array(value);
        if (candidate.length !== generics.LENGTH(this)) {
            throw new Error("'width' should be replaced by array of the same length");
        }
        IRanges.#checkNonNegative(candidate, "width");
        this._width = candidate;
        return this;
    }

    /**
     * @param {?Array} value - Array of strings containing a name for each range.
     * This should have length equal to the number of ranges.
     * Alternatively `null`, if no names are present.
     * @return {IRanges} A reference to this IRanges object, after setting the names to `value`.
     */
    $setNames(value) {
        this._rangeMetadata.$setRowNames(value);
        return this;
    }

    /**
     * @param {?DataFrame} value - Arbitrary metadata for each range.
     * This should have number of rows equal to the number of ranges.
     * Alternatively `null`, in which case all existing per-range metadata is removed.
     * @return {IRanges} A reference to this IRanges object, after setting the widths to `value`.
     */
    $setRangeMetadata(value) {
        if (value !== null) {
            if (!(value instanceof df.DataFrame) || generics.LENGTH(value) !== generics.LENGTH(this)) {
                throw new Error("'rangeMetadata' should be replaced by a DataFrame with the same number of rows");
            }
        } else {
            let existing = this._rangeMetadata;
            value = new df.DataFrame({}, { rowNames: existing.rowNames(), numberOfRows: existing.numberOfRows() });
        }
        this._rangeMetadata = value;
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

