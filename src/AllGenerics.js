export function _length(x) {
    if ("_bioconductor_length" in x) {
        return x._bioconductor_length();
    }

    if ("length" in x) {
        let out = x.length;
        let type = typeof out;
        if (type == "function") {
            return out();
        } else if (type == "number") {
            return out;
        }
    }
    
    throw new Error("no method for 'length' in '" + x.constructor.name + "' instance");
}

export function _slice(x, i, { allowInPlace = false, allowView = false } = {}) {
    if ("_bioconductor_slice" in x) {
        return x._bioconductor_slice(i, { allowInPlace, allowView });
    }

    // Otherwise, 'x' is probably some TypedArray or Array.
    if (i.constructor == Object) {
        if (i.end >= i.start) {
            if (allowView) {
                return x.subarray(i.start, i.end);
            } else {
                return x.slice(i.start, i.end);
            }
        } else {
            let output = x.slice(i.end, i.start);
            return output.reverse();
        }
    } else {
        let output = new x.constructor(i.length);
        i.forEach((y, j) => {
            output[j] = x[y];
        });
        return output;
    }
}

export function _combine(x, y, { allowAppend = false } = {}) {
    if ("_bioconductor_combine" in x) {
        return x._biocondutor_combine(y, { allowAppend });
    }

    // Otherwise, 'x' is probably some TypedArray or Array. It is assumed that
    // every 'y' is of some compatible Array-like type as well.
    let everything = [x, ...y];
    let total_length = 0;
    for (const obj of everything) {
        total_length += obj.length;
    }

    let output = new x.constructor(total_length);
    let position = 0;
    for (const obj of everything) {
        output.set(obj, position);
        position += obj.length;
    }

    return output;
}
