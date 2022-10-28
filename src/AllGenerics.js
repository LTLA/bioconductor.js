import * as utils from "./utils.js";

/**
 * Compute the length of a vector-like object.
 * For built-in Array and TypedArrays, this just returns the `length` property directly.
 * Custom classes should provide a `_bioconductor_LENGTH` method to describe their length.
 *
 * @param {*} x - Some vector-like object.
 * @return {number} Length of the object.
 */
export function LENGTH(x) {
    if ("_bioconductor_LENGTH" in x) {
        return x._bioconductor_LENGTH();
    }

    if (!utils.isArrayLike(x)) {
        throw new Error("no method for 'LENGTH' in '" + x.constructor.name + "' instance");
    }

    return x.length;
}

/**
 * Slice a vector-like object.
 * For built-in Array and TypedArrays, this just uses `slice()` or `subarray()`.
 * Custom classes should provide a `_bioconductor_SLICE` method to create a slice.
 *
 * @param {*} x - Some vector-like object.
 * @param {Object|Array|TypedArray} i - An Array or TypedArray of integer indices specifying the slice of `x` to retain.
 *
 * Alternatively, an object containing `start` and `end`, where the slice is defined as the sequence of consecutive integers in `[start, end)`.
 * @param {Object} [options={}] - Optional parameters.
 * @param {boolean} [options.allowInPlace=false] - Whether the slice is allowed to be performed in-place, thus modifying `x` by reference.
 * Whether this is actually done depends on the method, but may improve efficiency by avoiding unnecessary copies.
 * Setting `allowInPlace=true` should only be used if the input `x` is no longer needed.
 * @param {boolean} [options.allowView=false] - Whether a view can be created to mimic the slice operation.
 * Whether this is actually done depends on the method, but may improve efficiency by avoiding unnecessary copies.
 *
 * @return {*} A vector-like object, typically of the same class as `x`, containing data for the specified slice.
 *
 * If `allowInPlace = true`, `x` _may_ be modified in place, and the return value _may_ be a reference to `x`. 
 */
export function SLICE(x, i, { allowInPlace = false, allowView = false } = {}) {
    if ("_bioconductor_SLICE" in x) {
        return x._bioconductor_SLICE(i, { allowInPlace, allowView });
    }

    if (!utils.isArrayLike(x)) {
        throw new Error("no method for 'SLICE' in '" + x.constructor.name + "' instance");
    }

    if (i.constructor == Object) {
        if (allowView && ArrayBuffer.isView(x)) {
            return x.subarray(i.start, i.end);
        } else {
            return x.slice(i.start, i.end);
        }
    } else {
        let output = new x.constructor(i.length);
        i.forEach((y, j) => {
            output[j] = x[y];
        });
        return output;
    }
}

/**
 * Combine multiple vector-like objects.
 * Custom classes should provide a `_bioconductor_COMBINE` method to define the combining operation.
 *
 * @param {*} x - Some vector-like object.
 * @param {Array} y - Array of vector-like objects that are compatible with `x`.
 * @param {Object} [options={}] - Optional parameters.
 * @param {boolean} [options.allowAppend=false] - Whether the method can append elements of `y` onto `x`, thus modifying `x` in place.
 * Whether this is actually done depends on the method, but may improve efficiency by avoiding unnecessary copies.
 * Setting `allowAppend=true` should only be used if the input `x` is no longer needed.
 *
 * @return {*} A vector-like object containing the concatenated data from the input objects.
 * - If `x` is an instance of a custom class, the return value should be of the same class.
 * - If `x` and all elements of `y` are TypedArrays of the same class, the return value will be a TypedArray of that class.
 * - If any of `x` or the elements of `y` are Arrays, the return value will be an Array.
 * - If any of `x` or `y` are 64-bit TypedArrays of different classes, the return value will be an Array.
 * - Otherwise, for any other classes of TypedArrays in `x` or `y`, the return value will be a Float64Array.
 *
 * If `allowAppend=true`, the return value _may_ be a reference to `x` that is modified in-place.
 */
export function COMBINE(x, y, { allowAppend = false } = {}) {
    if ("_bioconductor_COMBINE" in x) {
        return x._bioconductor_COMBINE(y, { allowAppend });
    }

    if (!utils.isArrayLike(x)) {
        throw new Error("no method for 'COMBINE' in '" + x.constructor.name + "' instance");
    }

    // It is assumed that every 'y' is of some compatible Array-like type as well.
    let everything = [x, ...y];
    let total_LENGTH = 0;
    let constructor = x.constructor;

    for (const obj of everything) {
        total_LENGTH += obj.length;
        constructor = utils.chooseArrayConstructors(constructor, obj.constructor);
    }

    let output = new constructor(total_LENGTH);
    let position = 0;
    for (const obj of everything) {
        if ("set" in output) {
            output.set(obj, position);
            position += obj.length;
        } else {
            obj.forEach(x => {
                output[position] = x;
                position++;
            });
        }
    }

    return output;
}

/**
 * Clone a vector-like object.
 * Custom classes should provide a `_bioconductor_CLONE` method to define the cloning operation.
 *
 * @param {*} x - Some vector-like object.
 * @param {Object} [options={}] - Optional parameters.
 * @param {boolean} [options.deepCopy=true] - Whether to create a deep copy.
 * The exact interpretation of `deepCopy=false` is left to each method, but generally speaking, 
 * any setter (`$`-marked) functions operating on the copy should not alter `x`.
 *
 * @return {*} A clone of `x`, i.e., the return value and `x` should not compare equal.
 * If `deepCopy=true`, all internal components are also cloned.
 */
export function CLONE(x, { deepCopy = true } = {}) {
    if ("_bioconductor_CLONE" in x) {
        return x._bioconductor_CLONE({ deepCopy });
    }

    if (utils.isArrayLike(x)) {
        if (x.constructor == Array) {
            return x.slice();
        } else if (deepCopy) {
            return x.slice();
        } else {
            return x.subarray();
        }
    }

    if (x.constructor == Object) {
        if (deepCopy) {
            let output = {};
            for (const [k, v] of Object.entries(x)) {
                output[k] = CLONE(v);
            }
            return output;
        } else {
            return { ...x };
        }
    }

    // Immutable atomics should be all that's left.
    return x;
}
