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

