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

    buildOverlapIndex({ slice = null } = {}) {
        let n = (slice == null ? this._start.length : slice.length);
        let positions = new Int32Array(n * 2);
        let add = new Uint8Array(n * 2);
        let indices = new Int32Array(n * 2);

        let counter = 0;
        let fillIndex = i => {
            let at = counter * 2;
            let next = at + 1;
            positions[at] = this._start[i];
            positions[next] = this._start[i] + this._width[i];
            add[at] = 1;
            add[next] = 0;
            indices[at] = i;
            indices[next] = 0;
            counter++;
        };

        if (slice === null) {
            for (var i = 0; i < n; i++) {
                fillIndex(i);                                
            }
        } else {
            for (const i of slice) {
                fillIndex(i);
            }
        }

        // Sort by position, then by add-action. The latter is necessary to
        // ensure that add actions get processed before the corresponding
        // deletion action, otherwise there'd be nothing to delete!
        let order = utils.createSequence(positions.length);
        order.sort((i, j) => {
            let p = positions[i] - positions[j];
            return (p == 0 ? add[i] - add[j] : p); 
        });

        return new iro.IRangesOverlapIndex(
            generics.SLICE(positions, order), 
            generics.SLICE(add, order), 
            generics.SLICE(indices, order)
        );
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

export class IRangesOverlapIndex {
    constructor(position, add, index) {
        this._position = position;
        this._add = add;
        this._index = index;
    }

    overlap(queryIndex) {
        let self_position = this._position;
        let self_add = this._add;
        let self_index = this._index;
        let query_position = queryIndex._position;
        let query_add = queryIndex._add;
        let query_index = queryIndex._index;

        // Iterating through both indices and collecting overlaps.
        let overlaps = [];
        let self_active = new Set;
        let query_active = new Set;
        let i = 0;
        let j = 0;

        while (i < self_position.length && j < query_position.length) {
            if (self_position[i] < query_position[j]) {
                let sx = self_index[i];
                if (self_add[i]) {
                    self_active.add(sx);
                    query_active.forEach(q => {
                        overlaps.push([sx, q]);
                    });
                } else {
                    self_active.delete(sx);
                }
                i++;

            } else if (self_position[i] > query_position[j]) {
                let qx = query_index[j];
                if (query_add[j]) {
                    query_active.add(qx);
                    self_active.forEach(s => {
                        overlaps.push([s, qx]);
                    });
                } else {
                    query_active.delete(qx);
                }
                j++;

            } else {
                let current = self_position[i];
                let self_just_added = new Set;
                let query_just_added = new Set;

                // First, we fully process all additions/deletions on this
                // position.  Specifically, we need to process the ends of
                // ranges so that they don't get overlapped with new ranges
                // that are starting at this position.
                while (i < self_position.length && self_position[i] == current) {
                    let sx = self_index[i];
                    if (self_add[i]) {
                        self_just_added.add(sx);
                        self_active.add(sx);
                    } else {
                        self_active.delete(sx);
                    }
                    i++;
                }

                while (j < query_position.length && query_position[j] == current) {
                    let qx = query_index[j];
                    if (query_add[j]) {
                        query_just_added.add(qx);
                        query_active.add(qx);
                    } else {
                        query_active.delete(qx);
                    }
                    j++;
                }

                // Then we add the overlaps for all the newly added ranges on
                // this position. Note that this includes zero-length ranges
                // that would have already been removed from self/query_active.
                self_just_added.forEach(s => {
                    query_just_added.forEach(q => {
                        overlaps.push([s, q]);
                    });
                    query_active.forEach(q => {
                        overlaps.push([s, q]);
                    });
                });

                query_just_added.forEach(q => {
                    self_active.forEach(s => {
                        overlaps.push([s, q]);
                    });
                });
            }
        }

        return overlaps;
    }
}
