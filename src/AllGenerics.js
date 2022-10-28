import * as utils from "./utils.js";

export function LENGTH(x) {
    if ("_bioconductor_LENGTH" in x) {
        return x._bioconductor_LENGTH();
    }

    if (!utils.isArrayLike(x)) {
        throw new Error("no method for 'LENGTH' in '" + x.constructor.name + "' instance");
    }

    return x.length;
}

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
