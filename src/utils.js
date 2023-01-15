export function areArraysEqual(x, y) {
    if (x.length !== y.length) {
        return false;
    }

    for (var i = 0; i < x.length; i++) {
        if (x[i] != y[i]) {
            return false;
        }
    }

    return true;
}

export function isArrayLike(x) {
    return x.constructor == Array || ArrayBuffer.isView(x);
}

export function chooseArrayConstructors(con1, con2) {
    if (con1 == con2) {
        return con1;
    }

    if (con1 == Array || con2 == Array) {
        return Array;
    }

    if (con1 == BigInt64Array || con2 == BigInt64Array || con1 == BigUint64Array || con2 == BigUint64Array) {
        return Array;
    }

    return Float64Array;
}

export function formatLengthError(left, right) {
    return new Error(left + " should have length equal to " + right);
}

export function checkNamesArray(names, typeMessage, numExpected, lengthMessage) {
    for (const x of names) {
        if (typeof x !== "string") {
            throw new Error(typeMessage + " array should only contain strings");
        }
    }
    if (names.length != numExpected) {
        throw formatLengthError(typeMessage + " array ", lengthMessage);
    }
}

export function sum(y) {
    let total = 0;
    y.forEach(x => { total += x; });
    return total;
}

export function combineNames(all_names, all_lengths, total_n = null) {
    let all_null = true;
    for (var i = 0; i < all_names.length; i++) {
        if (all_names[i] !== null) {
            all_null = false;
        }
    }

    if (all_null) {
        return null;
    }

    if (total_n === null) {
        total_n = sum(all_lengths);
    }

    let output = new Array(total_n);
    let counter = 0;
    for (var i = 0; i < all_names.length; i++) {
        let n = all_names[i];
        if (n === null) {
            output.fill("", counter, counter + all_lengths[i]);
            counter += all_lengths[i];
        } else {
            n.forEach(x => {
                output[counter] = x;
                counter++;
            });
        }
    }

    return output;
}
