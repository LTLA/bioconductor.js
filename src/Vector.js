import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as df from "./DataFrame.js";
import * as ann from "./Annotated.js";

function verifyElementMetadata(elementMetadata, numExpected, className) {
    if (elementMetadata !== null) {
        if (!(elementMetadata instanceof df.DataFrame)) {
            throw new Error("'elementMetadata' should be a DataFrame");
        }
        if (generics.LENGTH(elementMetadata) !== numExpected) {
            throw new Error("'elementMetadata' should have the same number of rows as 'LENGTH(<" + className + ">)'");
        }
    } else {
        elementMetadata = new df.DataFrame({}, { numberOfRows: numExpected });
    }
    return elementMetadata;
}

/**
 * The Vector class implements a store for arbitrary per-element metadata and per-element names.
 * It is intended as a base class for other structures that have a concept of "vector-ness".
 * It should not be constructed directly.
 *
 * @augments Annotated
 */
export class Vector extends ann.Annotated {
    /**
     * @param {number} length - Number of elements in this vector-like object.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?DataFrame} [options.elementMetadata=null] - A {@linkplain DataFrame} with number of rows equal to the length of `start`, containing arbitrary per-element annotations.
     * Alternatively `null`, in which case a zero-column DataFrame is automatically constructed.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(length, { names = null, elementMetadata = null, metadata = {} } = {}) {
        if (arguments.length == 0) {
            super();
            return;
        }

        super(metadata);

        this._elementMetadata = verifyElementMetadata(elementMetadata, length, this.constructor.className);

        if (names !== null) {
            utils.checkNamesArray(names, "'names'", length, "'LENGTH(<" + this.constructor.className + ">)'");
        }
        this._names = names;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {DataFrame} A {@linkplain DataFrame} with one row corresponding to each vector element, containing arbitrary per-element metadata.
     */
    elementMetadata() {
        return this._elementMetadata;
    }

    /**
     * @return {?Array} Array of strings containing the name of each range, or `null` if no names are available.
     */
    names() {
        return this._names;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {?DataFrame} elementMetadata - Arbitrary metadata for each vector element.
     * This should have number of rows equal to the vector length.
     * Alternatively `null`, in which case all existing per-element metadata is removed.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this Vector instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {Vector} The Vector object after setting the element metadata to `value`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setElementMetadata(elementMetadata, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        target._elementMetadata = verifyElementMetadata(elementMetadata, generics.LENGTH(target), target.constructor.className);
        return target;
    }

    $setElementMetadata(elementMetadata) {
        return this.setElementMetadata(elementMetadata, { inPlace: true });
    }

    /**
     * @param {?Array} names - Array of strings containing a name for each range.
     * This should have length equal to the number of ranges.
     * Alternatively `null`, if no names are present.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this Vector instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {Vector} The Vector object after setting the names to `value`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setNames(names, { inPlace = false } = {}) {
        if (names !== null) {
            utils.checkNamesArray(names, "replacement 'names'", generics.LENGTH(this), "'LENGTH(<" + this.constructor.className + ">)'");
        } 
        let target = cutils.setterTarget(this, inPlace);
        target._names = names;
        return target;
    }

    $setNames(names) {
        return this.setNames(names, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_SLICE(output, i, { allowView = false }) {
        output._elementMetadata = generics.SLICE(this._elementMetadata, i, { allowView });
        output._names = (this._names === null ? null : generics.SLICE(this._names, i, { allowView }));
        output._metadata = this._metadata;
        return;
    }

    _bioconductor_COMBINE(output, objects) {
        let all_em = [];
        let all_n = [];
        let all_l = [];

        for (const x of objects) {
            all_em.push(x._elementMetadata);
            all_n.push(x._names);
            all_l.push(generics.LENGTH(x));
        }

        output._elementMetadata = generics.COMBINE(all_em);
        output._names = utils.combineNames(all_n, all_l);
        return;
    }

    _bioconductor_CLONE(output, { deepCopy = true }) {
        super._bioconductor_CLONE(output, { deepCopy });
        output._elementMetadata = cutils.cloneField(this._elementMetadata, deepCopy);
        output._names = cutils.cloneField(this._names, deepCopy);
        return;
    }
}
