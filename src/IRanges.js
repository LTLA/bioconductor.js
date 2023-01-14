import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";
import * as df from "./DataFrame.js";

/**
 * An IRanges object is a collection of integer ranges.
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

    start() {
        return this._start;
    }

    end() {
        return this._start.map((x, i) => x + this._width[i]);
    }

    width() {
        return this._width;
    }

    names() {
        return this._rangeMetadata.rowNames();
    }

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

    $setStart(value) {
        let candidate = IRanges.#convertToInt32Array(value);
        if (candidate.length !== generics.LENGTH(this)) {
            throw new Error("'start' should be replaced by array of the same length");
        }
        IRanges.#checkNonNegative(candidate, "start");
        this._start = candidate;
        return;
    }

    $setWidth(value) {
        let candidate = IRanges.#convertToInt32Array(value);
        if (candidate.length !== generics.LENGTH(this)) {
            throw new Error("'width' should be replaced by array of the same length");
        }
        IRanges.#checkNonNegative(candidate, "width");
        this._width = candidate;
        return;
    }

    $setNames(value) {
        this._rangeMetadata.$setRowNames(value);
        return;
    }

    $setRangeMetadata(value) {
        if (!(value instanceof df.DataFrame) || generics.LENGTH(value) !== generics.LENGTH(this)) {
            throw new Error("'rangeMetadata' should be replaced by a DataFrame with the same number of rows");
        }
        this._rangeMetadata = value;
        return;
    }

    $setMetadata(value) {
        this._metadata = value;
        return;
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

