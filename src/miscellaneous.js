/**
 * @param {Array|TypedArray} x - Array of values to be interpreted as truthy or falsey.
 * @param {Object} [options={}] - Optional parameters.
 * @param {boolean} [options.not=false] - Whether to select the entries of `x` that are falsey.
 *
 * @return {Array} Array of indices of the entries of `x` that are truthy (if `not=false`) or falsey (if `not=true`).
 * This array is guaranteed to be sorted in ascending order.
 */
export function which(x, { not = false } = {}) {
    let output = [];
    x.forEach((y, i) => {
        if ((!y) == not) {
            output.push(i);
        }
    });
    return output;
}

/**
 * Given a factor, return the indices corresponding to each level.
 * This can be used in subsequent {@linkcode splitRows} calls.
 *
 * @param {Array|TypedArray} factor - Array containing the factor of interest.
 *
 * @return {Object} Object where each key is a factor level and each value is an array containing the indices corresponding to that level in `factor`.
 */
export function presplitFactor(factor) {
    let by = {};
    factor.forEach((x, i) => {
        if (!(x in by)) {
            by[x] = [];
        }
        by[x].push(i);
    });
    return by;
}
