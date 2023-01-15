import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";
import * as df from "./DataFrame.js";
import * as vec from "./Vector.js";

/**
 * An IRanges object is a collection of integer ranges, inspired by the class of the same name from the Bioconductor ecosystem.
 * Each range consists of a start position and a width, and may be associated with arbitrary range-level metadata in a {@linkplain DataFrame}.
 * The IRanges defines methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 *
 * @extends Vector
 */
export class IRanges extends vec.Vector {
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
     * @param {?DataFrame} [options.elementMetadata=null] - A {@linkplain DataFrame} with number of rows equal to the length of `start`, containing arbitrary per-range annotations.
     * Alternatively `null`, in which case a zero-column DataFrame is automatically constructed.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(start, width, { names = null, elementMetadata = null, metadata = {} } = {}) {
        super(start.length, { names, elementMetadata, metadata });

        this._start = IRanges.#convertToInt32Array(start);
        IRanges.#checkNonNegative(this._start, "start");

        this._width = IRanges.#convertToInt32Array(width);
        IRanges.#checkNonNegative(this._width, "width");

        let n = this._start.length;
        if (n !== this._width.length) {
            throw new Error("'start' and 'width' should have the same length");
        }
    }

    static className = "IRanges";

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

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_LENGTH() {
        return this._start.length;
    }

    _bioconductor_SLICE(output, i, { allowView = false }) {
        super._bioconductor_SLICE(output, i, { allowView });
        output._start = generics.SLICE(this._start, i, { allowView });
        output._width = generics.SLICE(this._width, i, { allowView });
        return;
    }

    _bioconductor_COMBINE(output, objects) {
        super._bioconductor_COMBINE(output, objects);

        let all_s = [];
        let all_w = [];
        for (const x of objects) {
            all_s.push(x._start);
            all_w.push(x._width);
        }

        output._start = generics.COMBINE(all_s);
        output._width = generics.COMBINE(all_w);
        return;
    }

    _bioconductor_CLONE(output, { deepcopy = true }) {
        super._bioconductor_CLONE(output, { deepcopy });
        output._start = generics.CLONE(this._start, { deepcopy });
        output._width = generics.CLONE(this._width, { deepcopy });
        return;
    }
}

