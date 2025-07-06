import * as utils from "./utils.js";
import * as misc from "./miscellaneous.js";

/**
 * Compute the length of a vector-like object.
 *
 * For Array and TypedArrays, this just returns the `length` property directly.
 *
 * Custom classes should implement a `_bioconductor_LENGTH` method to return their length via this generic.
 * This method should accept no arguments. 
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
 *
 * For Array and TypedArrays, this just uses `slice()` or `subarray()`.
 *
 * Custom classes should implement a `_bioconductor_SLICE` method to create a slice via this generic.
 * This method should accept `i` and `options`, and return a appropriately sliced instance of the same class as `x`.
 *
 * @param {*} x - Some vector-like object.
 * @param {Object|Array|TypedArray} i - An Array or TypedArray of integer indices specifying the slice of `x` to retain.
 *
 * Alternatively, an object containing `start` and `end`, where the slice is defined as the sequence of consecutive integers in `[start, end)`.
 * @param {Object} [options={}] - Optional parameters.
 * @param {boolean} [options.allowView=false] - Whether a view can be created to mimic the slice operation.
 * Whether this is actually done depends on the method, but may improve efficiency by avoiding unnecessary copies.
 *
 * @return {*} A vector-like object of the same class as `x`, containing data for the specified slice.
 *
 * If `allowInPlace = true`, `x` _may_ be modified in place, and the return value _may_ be a reference to `x`. 
 */
export function SLICE(x, i, { allowView = false } = {}) {
    if ("_bioconductor_SLICE" in x) {
        return x._bioconductor_SLICE(i, { allowView });
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
 *
 * For Array and TypedArrays, the combined array is of a class that avoids information loss.
 *
 * Custom classes should implement a `_bioconductor_COMBINE` method to combine objects via this generic. 
 * This method will be invoked from the first instance of `objects` and should accept an array of remaining objects (i.e., `objects.slice(1)`) to be combined to the first instance. 
 * It return an appropriately-combined instance of the same class as the first instance.
 *
 * @param {Array} objects - Array of vector-like objects to be combined.
 * It is assumed that the objects are of the same class, or at least compatible with each other -
 * for custom classes, the definition of "compatibility" depends on the `_bioconductor_COMBINE` method of the first element of `objects`.
 *
 * @return {*} A vector-like object containing the concatenated data from the input objects.
 * - If the first entry of `objects` is an instance of a custom class, the return value should be of the same class.
 * - If all `objects` are TypedArrays of the same class, the return value will be a TypedArray of that class.
 * - If any of the `objects` are Arrays, the return value will be an Array.
 * - If any of the `objects` are 64-bit TypedArrays of different classes, the return value will be an Array.
 * - Otherwise, for any other classes of TypedArrays in `objects`, the return value will be a Float64Array.
 */
export function COMBINE(objects) {
    let x = objects[0];
    if ("_bioconductor_COMBINE" in x) {
        return x._bioconductor_COMBINE(objects.slice(1));
    }

    if (!utils.isArrayLike(x)) {
        throw new Error("no method for 'COMBINE' in '" + x.constructor.name + "' instance");
    }

    // It is assumed that every 'y' is of some compatible Array-like type as well.
    let total_LENGTH = 0;
    let constructor = x.constructor;

    for (const obj of objects) {
        total_LENGTH += obj.length;
        constructor = utils.chooseArrayConstructors(constructor, obj.constructor);
    }

    let output = new constructor(total_LENGTH);
    let position = 0;
    for (const obj of objects) {
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
 * 
 * For TypedArrays, this just uses `slice()`.
 * For Arrays, this creates a copy and runs `CLONE` on each element in the copy.
 *
 * Custom classes should provide a `_bioconductor_CLONE` method to define the cloning operation.
 * This method should accept `options` and should return an appropriately cloned instance of the same class as `x`.
 *
 * @param {*} x - Some vector-like object.
 * @param {Object} [options={}] - Optional parameters.
 * @param {boolean} [options.deepCopy=true] - Whether to create a deep copy.
 * The exact interpretation of `deepCopy=false` is left to each method.
 *
 * @return {*} A clone of `x`, i.e., the return value and `x` should not compare equal.
 * If `deepCopy=true`, all internal components are also cloned.
 */
export function CLONE(x, { deepCopy = true } = {}) {
    if (x instanceof Object) {
        let options = { deepCopy };
        if ("_bioconductor_CLONE" in x) {
            return x._bioconductor_CLONE(options);
        }

        if (utils.isArrayLike(x)) {
            if (x.constructor == Array) {
                return x.map(y => CLONE(y, options));
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

        if (x.constructor == Map) {
            let output = new Map;
            for (const [k, v] of x) {
                output.set(k, deepCopy ? CLONE(v) : v);
            }
            return output;
        }

        if (x.constructor == Set) {
            let output = new Set;
            for (const k of x) {
                output.add(deepCopy ? CLONE(k) : k);
            }
            return output;
        }


        throw new Error("unknown CLONE operation for instance of class '" + x.constructor.name + "'");
    }

    // Immutable atomics should be all that's left.
    return x;
}

/**
 * Split a vector-like object along its length according to the levels of a factor of the same length.
 * This works automatically for all classes for which there is a {@linkcode SLICE} method.
 * Custom classes may also choose to define their own `_bioconductor_SPLIT` method,
 * which should accept `factor` and return an array of sliced instances of the same class as `x`.
 *
 * @param {*} x - Some vector-like object.
 * @param {Array|TypedArray} factor - Array containing the factor to use for splitting.
 * This should have the same length as `x`.
 *
 * Alternatively, the output of {@linkcode presplitFactor} can be supplied.
 *
 * @return {Object} An object containing one key per level of `factor`,
 * where the value is the slice of `x` corresponding to the indices of that level in `factor`.
 */
export function SPLIT(x, factor) {
    if (factor.constructor != Object) {
        factor = misc.presplitFactor(factor);
    }

    if ("_bioconductor_SPLIT" in x) {
        return x._bioconductor_SPLIT(factor);
    }

    let output = {};
    for (const [k, v] of Object.entries(factor)) {
        output[k] = SLICE(x, v);
    }

    return output;
}

/**
 * Return the number of rows for a two-dimensional object.
 * Custom classes should implement a `_bioconductor_NUMBER_OF_ROWS` method, accepting no arguments.
 *
 * @param {*} x - Some two-dimensional object.
 * @return {number} Number of rows.
 */
export function NUMBER_OF_ROWS(x) {
    if (!("_bioconductor_NUMBER_OF_ROWS" in x)) {
        throw new Error("no 'NUMBER_OF_ROWS' method available for '" + x.constructor.name + "' instance");
    }
    return x._bioconductor_NUMBER_OF_ROWS();
}

/**
 * Return the number of columns for a two-dimensional object.
 * Custom classes should implement a `_bioconductor_NUMBER_OF_COLUMNS` method, accepting no arguments.
 *
 * @param {*} x - Some two-dimensional object.
 * @return {number} Number of columns.
 */
export function NUMBER_OF_COLUMNS(x) {
    if (!("_bioconductor_NUMBER_OF_COLUMNS" in x)) {
        throw new Error("no 'NUMBER_OF_COLUMNS' method available for '" + x.constructor.name + "' instance");
    }
    return x._bioconductor_NUMBER_OF_COLUMNS();
}

/**
 * Slice a two-dimensional object by its rows and/or columns.
 *
 * Custom classes should implement a `_bioconductor_SLICE_2D` method to slice rows/columns via this generic.
 * This method should accept `rows`, `columns` and `options`, and should return a sliced instance of the same class as `x`.
 *
 * @param {*} x - Some two-dimensional object.
 * @param {?(Object|Array|TypedArray)} rows - An Array or TypedArray of integer indices specifying the row-wise slice of `x` to retain.
 *
 * Alternatively, an object containing `start` and `end`, where the slice is defined as the sequence of consecutive integers in `[start, end)`.
 * 
 * Alternatively `null`, to indicate that no slicing is to be performed on the rows.
 * @param {?(Object|Array|TypedArray)} columns - An Array or TypedArray of integer indices specifying the column-wise slice of `x` to retain.
 *
 * Alternatively, an object containing `start` and `end`, where the slice is defined as the sequence of consecutive integers in `[start, end)`.
 *
 * Alternatively `null`, to indicate that no slicing is to be performed on the columns.
 * @param {Object} [options={}] - Optional parameters.
 * @param {boolean} [options.allowView=false] - Whether a view can be created to mimic the slice operation.
 * Whether this is actually done depends on the method, but may improve efficiency by avoiding unnecessary copies.
 *
 * @return {*} A two-dimensional object, typically of the same class as `x`, containing data for the specified slice.
 */
export function SLICE_2D(x, rows, columns, { allowView = false } = {}) {
    if (!("_bioconductor_SLICE_2D" in x)) {
        throw new Error("no 'SLICE_2D' method available for '" + x.constructor.name + "' instance");
    }
    return x._bioconductor_SLICE_2D(rows, columns, { allowView });
}

/**
 * Combine multiple two-dimensional objects by row.
 *
 * Custom classes should implement a `_bioconductor_COMBINE_ROWS` method to combine objects via this generic.
 * This method will be invoked from the first instance of `objects` and should accept an array of remaining objects (i.e., `objects.slice(1)`) to be combined to the first instance. 
 * It return an appropriately-combined instance of the same class as the first instance.
 *
 * @param {Array} objects - Array of two-dimensional objects to be combined by row.
 * It is assumed that the objects are of the same class, or at least compatible with each other -
 * for custom classes, the definition of "compatibility" depends on the `_bioconductor_COMBINE_ROWS` method of the first element of `objects`.
 *
 * @return {*} A two-dimensional object of the same class as the first entry of `objects`, ontaining the row-wise concatenated data from all `objects`.
 */
export function COMBINE_ROWS(objects) {
    let x = objects[0];
    if (!("_bioconductor_COMBINE_ROWS" in x)) {
        throw new Error("no 'COMBINE_ROWS' method available for '" + x.constructor.name + "' instance");
    }
    return x._bioconductor_COMBINE_ROWS(objects.slice(1));
}

/**
 * Combine multiple two-dimensional objects by column.
 *
 * Custom classes should implement a `_bioconductor_COMBINE_COLUMNS` method to combine objects via this generic.
 * This method will be invoked from the first instance of `objects` and should accept an array of remaining objects (i.e., `objects.slice(1)`) to be combined to the first instance. 
 * It return an appropriately-combined instance of the same class as the first instance.
 *
 * @param {Array} objects - Array of two-dimensional objects to be combined by column.
 * It is assumed that the objects are of the same class, or at least compatible with each other -
 * for custom classes, the definition of "compatibility" depends on the `_bioconductor_COMBINE_COLUMNS` method of the first element of `objects`.
 *
 * @return {*} A two-dimensional object of the same class as the first entry of `objects`, ontaining the column-wise concatenated data from all `objects`.
 */
export function COMBINE_COLUMNS(objects) {
    let x = objects[0];
    if (!("_bioconductor_COMBINE_COLUMNS" in x)) {
        throw new Error("no 'COMBINE_COLUMNS' method available for '" + x.constructor.name + "' instance");
    }
    return x._bioconductor_COMBINE_COLUMNS(objects.slice(1));
}
