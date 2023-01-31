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

export function checkStringArray(names, typeMessage) {
    for (const x of names) {
        if (typeof x !== "string") {
            throw new Error(typeMessage + " array should only contain strings");
        }
    }
}

export function checkNamesArray(names, typeMessage, numExpected, lengthMessage) {
    checkStringArray(names, typeMessage);
    if (names.length != numExpected) {
        throw formatLengthError(typeMessage + " array", lengthMessage);
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

export function createSequence(n) {
    let output = new Int32Array(n);
    for (var i = 0; i < n; i++) {
        output[i] = i;
    }
    return output;
}

export function isSorted(n, cmp) {
    for (var i = 1; i < n; ++i) {
        if (cmp(i-1, i) > 0) {
            return false;
        }
    }
    return true;
}

export function convertToInt32Array(x) {
    if (x instanceof Int32Array) {
        return x;
    } else {
        return new Int32Array(x);
    }
}

export function checkNonNegative(x, msg) {
    for (const y of x) {
        if (y < 0) {
            throw new Error("detected a negative entry in '" + msg + "'");
        }
    }
}

export function checkEntryOrder(entries, order, name) {
    let expected = Object.keys(entries);
    if (order == null) {
        return expected;
    }

    checkNamesArray(order, "'" + name + "Order'", expected.length, "the number of entries in '" + name + "s'");
    let observed = order.slice().sort();
    expected.sort();

    if (!areArraysEqual(observed, expected)) {
        throw new Error("values of '" + name + "Order' should be the same as the keys of '" + name + "s'");
    }

    return order;
}

export function check_entry_index(order, i, fieldName, className) {
    if (i < 0 || i >= order.length) {
        throw new Error(fieldName + " index '" + String(i) + "' out of range for this " + className);
    }
}

export function retrieveSingleEntry(entries, order, i, fieldName, className) {
    if (typeof i == "string") {
        if (!(i in entries)) {
            throw new Error("no " + fieldName + "'" + i + "' present in this " + className);
        }
        return entries[i];
    } else {
        check_entry_index(order, i, fieldName, className);
        return entries[order[i]];
    }
}

export function removeSingleEntry(entries, order, i, fieldName, className) {
    if (typeof i == "string") {
        let ii = order.indexOf(i);
        if (ii < 0) {
            throw new Error("no " + fieldName + " '" + i + "' present in this " + className);
        }
        order.splice(ii, 1); // modified by reference.
        delete entries[i];
    } else {
        check_entry_index(order, i, fieldName, className);
        let n = order[i];
        order.splice(i, 1);
        delete entries[n];
    }
}

export function setSingleEntry(entries, order, i, value, fieldName, className) {
    if (typeof i == "string") {
        if (!(i in entries)) {
            order.push(i);
        }
        entries[i] = value;
    } else {
        check_entry_index(order, i, fieldName, className);
        entries[order[i]] = value;
    }
}

export function combineEntries(dumps, combiner, orderName, className) {
    let first_order = dumps[0].order;
    for (var i = 1; i < dumps.length; i++) {
        if (!areArraysEqual(first_order, dumps[i].order)) {
            throw new Error("mismatching '" + orderName + "' for " + className + " " + String(i) + " to be combined");
        }
    }

    let new_entries = {};
    for (const k of first_order) {
        let found = dumps.map(x => x.entries[k]);
        new_entries[k] = combiner(found);
    }

    return { entries: new_entries, order: first_order };
}
