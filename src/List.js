import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as generics from "./AllGenerics.js";

class IndexedNames {
    constructor(names) {
        this._names = names;
        this._lookup = new Map;
    }

    names() {
        return this._names;
    }

    nameToIndexUncached(name) {
        if (this._lookup.has(name)) {
            return this._lookup.get(name);
        }
        return this._names.indexOf(name);
    }

    nameToIndex(name, { error = true } = {}) {
        if (this._lookup.has(name)) {
            return this._lookup.get(name);
        }

        for (var i = this._lookup.size; i < this._names.length; i++) {
            const current = this._names[i];
            if (this._lookup.has(current)) {
                continue; // only keep the first instance of a duplicated name.
            }
            this._lookup.set(current, i);
            if (this._names[i] == name) {
                return i;
            }
        }

        if (error) {
            throw new Error("no matching name for '" + name + "' in this " + this.constructor.className);
        } else {
            return -1;
        }
    }

    indexToName(i) {
        return this._names[i];
    }

    append(name, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._names = target._names.slice();
            target._lookup = new Map; // always making a new map to avoid sharing lookup tables between instances with different names.
        }
        target._names.push(name);
        return target;
    }

    set(i, name, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._names = target._names.slice();
            target._lookup = new Map; // always making a new map to avoid sharing lookup tables between instances with different names.
        } else {
            if (target._lookup.size > i) { // if the lookup never got to that point, we don't have to wipe it.
                target._lookup = new Map;
            }
        }
        target._names[i] = name;
        return target;
    }

    delete(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._names = target._names.slice();
            target._lookup = new Map; // always making a new map to avoid sharing lookup tables between instances with different names.
        } else {
            if (target._lookup.size > i) { // i.e., if the lookup never got to that point, we don't have to wipe it.
                target._lookup = new Map;
            }
        }
        target._names.splice(i, 1);
        return target;
    }

    _bioconductor_CLONE({ deepCopy = true }) {
        let output = new this.constructor;
        output._names = cutils.cloneField(this._names, deepCopy);
        output._lookup = cutils.cloneField(this._lookup, deepCopy);
        return output;
    }
}

/**
 * An R-style list with optional names.
 * Callers can get/set individual elements by positional index or name.
 * Operations like slicing and combining will apply to both the values and names.
 *
 * The List defines methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 *
 * We explicitly allow duplicates in the names to avoid errors when slicing or combining.
 * Otherwise, it would be impossible to construct a slice with duplicate indices or to combine multiple `List` instances with shared names.
 *
 * Constructors of List subclasses should be callable with no arguments, possibly creating an empty object with no properties.
 * This will be used by the `_bioconductor_CLONE`, `_bioconductor_SLICE` and `_bioconductor_COMBINE` methods to return an instance of the subclass.
 */
export class List {
    /**
     * @param {Array|Map|Object} values - Elements of the List.
     * For Maps or objects, the values (in order of iteration) are used as the List elements.
     * If no argument is supplied, it defaults to an empty Array.
     * @param {Object} [options={}] - Further options.
     * @param {?Array} [options.names=null] - An array of strings containing the names of the List elements.
     * If provided, this should be of the same length as `values`.
     * If `values` is a Map or object, `names` should have the same keys.
     * If `values` is an array, the names may contain duplicate strings.
     * If `null` and `values` is an array, the List will be unnamed.
     */
    constructor(values, { names = null } = {}) {
        if (arguments.length == 0) {
            this._values = [];
            this._names = null;
            return;
        }

        if (values instanceof Array) {
            if (names !== null) {
                if (names.length != values.length) {
                    throw new Error("'names' and 'values' should have the same length");
                }
                for (const n of names) {
                    if (typeof n != "string") {
                        throw new Error("'names' should be an array of strings");
                    }
                }
                names = new IndexedNames(names);
            }

            this._values = values;
            this._names = names;

        } else if (values instanceof Map) {
            let arr = [];
            if (names == null) {
                names = [];
                for (const [n, v] of values) {
                    if (typeof n != "string") {
                        throw new Error("keys of 'values' should be strings");
                    }
                    names.push(n);
                    arr.push(v);
                }

            } else {
                if (names.length != values.size) {
                    throw new Error("size of 'values' should be equal to length of 'names'");
                }
                for (const n of names) {
                    if (typeof n != "string") {
                        throw new Error("'names' should be an array of strings");
                    }
                    if (!values.has(n)) {
                        throw new Error("missing name '" + n + "' in 'values'");
                    }
                    arr.push(values.get(n));
                }
            }

            this._values = arr;
            this._names = new IndexedNames(names);

        } else {
            let arr = [];
            if (names == null) {
                names = [];
                for (const [n, v] of Object.entries(values)) {
                    names.push(n);
                    arr.push(v);
                }

            } else {
                if (names.length != Object.keys(values).length) {
                    throw new Error("size of 'values' should be equal to length of 'names'");
                }
                for (const n of names) {
                    if (typeof n != "string") {
                        throw new Error("'names' should be an array of strings");
                    }
                    if (!(n in values)) {
                        throw new Error("missing name '" + n + "' in 'values'");
                    }
                    arr.push(values[n]);
                }
            }

            this._values = arr;
            this._names = new IndexedNames(names);
        }
    }

    /**
     * @return {?Array} Array of names of the List elements, or `null` if the List is unnamed.
     */
    names() {
        if (this._names == null) {
            return null;
        } else {
            return this._names.names();
        }
    }

    /**
     * @return {Array} Array containing the List elements.
     */
    values() {
        return this._values;
    }

    /**
     * @return {number} Length of the list.
     */
    length() {
        return this._values.length;
    }

    static className = "List";

    /***********************************************/

    #check_index(i) {
        if (i < 0 || i >= this._values.length) {
            throw new Error(" index '" + String(i) + "' out of range for this " + this.constructor.className);
        }
    }

    /**
     * @param {number} i - Index of the List element to retrieve.
     * This should be non-negative and less than {@linkcode List#length length}.
     * @return The `i`-th List element.
     */
    getByIndex(i) {
        this.#check_index(i);
        return this._values[i];
    }

    /**
     * @param {string} name - Name of the List element to retrieve.
     * This should be present in {@linkcode List#names names}.
     * @return The List element corresponding to `name`.
     * If duplicates of `name` are present in the list, the first occurrence is returned.
     */
    getByName(name) {
        if (this._names === null) {
            throw new Error("no available names in this " + this.constructor.className);
        }
        let candidate = this._names.nameToIndex(name);
        return this._values[candidate];
    }

    /**
     * @param {string|number} i - Index or name of the List element to retrieve.
     * Numbers are passed to {@linkcode List#getByIndex getByIndex} and strings are passed to {@linkcode List#getByName getByName}.
     * @return The List element at/for `i`.
     */
    get(i) {
        if (typeof i == "number") {
            return this.getByIndex(i);
        } else {
            return this.getByName(i);
        }
    }

    /**
     * @param {string} name - Name of a List element.
     * @return {boolean} Whether the name exists in this List.
     */
    has(name) {
        return this.nameToIndex(name) >= 0;
    }

    /**
     * @param {string} name - Name of a List element.
     * @return {number} Index of the name in {@linkcode List#names names}.
     * If duplicate names are present, the first occurrence is returned.
     * If the name is not present, -1 is returned.
     */
    nameToIndex(name) {
        return this._names.nameToIndex(name, { error: false });
    }

    /***********************************************/

    /**
     * @return {Array} Array of values, equivalent to {@linkcode List#values values}.
     */
    toArray() {
        return this._values;
    }

    /**
     * @return {Map} Map of name-value pairs.
     * If duplicate names are present, only the value for the first occurrence is reported.
     * If the List is unnamed, an error is thrown.
     */
    toMap() {
        if (this._names == null) {
            throw new Error("no available names in this '" + this.constructor.className + "'");
        }
        let output = new Map;
        let names = this._names.names();
        for (var i = 0; i < this._values.length; i++) {
            const curname = names[i];
            if (output.has(curname)) {
                continue;
            }
            output.set(curname, this._values[i]);
        }
        return output;
    }

    /**
     * @return {Object} Object of name-value pairs.
     * If duplicate names are present, only the value for the first occurrence is reported.
     * If the List is unnamed, an error is thrown.
     */
    toObject() {
        if (this._names == null) {
            throw new Error("no available names in this '" + this.constructor.className + "'");
        }
        let output = {};
        let names = this._names.names();
        for (var i = 0; i < this._values.length; i++) {
            const curname = names[i];
            if (curname in output) {
                continue;
            }
            output[curname] = this._values[i];
        }
        return output;
    }

    /***********************************************/

    /**
     * @param {number} i - Index of the List element to set.
     * This should be non-negative and no greater than {@linkcode List#length length}.
     * If `i` is less than `length`, the `i`-th element is replaced by `x`.
     * If `i` is equal to `length`, `x` is appended to the end of the list.
     * @param {*} x - Value of a List element.
     * @param {Object} [options={}] - Further options.
     * @param {?string} [options.name=null] - Name for the List element at `i`.
     * If `i` is less than `length`, the name of the `i`-th element is replaced by `name`.
     * If `i` is equal to `length`, the name of the newly-appended element is set to `name`.
     * If the List did not previously have any names, the names of all other elements are set to an empty string.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} The List after setting the `i`-th element to `x`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setByIndex(i, x, { name = null, inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._values = target._values.slice();
        }

        if (i < 0 || i > this._values.length) {
            throw new Error(" index '" + String(i) + "' out of range for this " + this.constructor.className);
        }

        if (i == target._values.length) {
            target._values.push(x);
            if (name == null) {
                if (target._names != null) {
                    target._names = target._names.append("", { inPlace });
                }

            } else {
                if (typeof name != "string") {
                    throw new Error("'name' should be a string");
                }
                if (target._names == null) {
                    const new_names = new Array(target._values.length).fill("");
                    new_names[i] = name;
                    target._names = new IndexedNames(new_names);
                } else {
                    target._names = target._names.append(name, { inPlace });
                }
            }

        } else {
            target._values[i] = x;
            if (name !== null) {
                if (target._names === null) {
                    const new_names = new Array(target._values.length).fill("");
                    new_names[i] = name;
                    target._names = new IndexedNames(new_names);
                } else {
                    target._names = target._names.set(i, name, { inPlace });
                }
            }
        }

        return target;
    }

    /**
     * @param {string} name - Name of the List element to set.
     * If this already exists in {@linkcode List#names names}, the corresponding element is replaced by `x`.
     * Otherwise, `x` is appended to the List with the name `name`.
     * If the List did not previously have any names, the names of all other elements are set to an empty string.
     * @param {*} x - Value of a List element.
     * @param {Object} [options={}] - Further options.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} The List after setting the `name`d entry to `x`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setByName(name, x, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._values = target._values.slice();
        }

        if (target._names !== null) {
            let candidate = target._names.nameToIndex(name, { error: false });
            if (candidate < 0) {
                target._values.push(x);
                target._names = target._names.append(name, { inPlace });
            } else {
                target._values[candidate] = x;
            }
        } else {
            const new_names = new Array(target._values.length).fill("");
            new_names.push(name);
            target._names = new IndexedNames(new_names);
            target._values.push(x);
        }

        return target;
    }

    /**
     * @param {string|number} i - Index or name of the list element to set.
     * Numbers are passed to {@linkcode List#setByIndex setByIndex} and strings are passed to {@linkcode List#setByName setByName}.
     * @param {*} x - Value of a List element.
     * @param {Object} [options={}] - Further options.
     * @param {?string} [options.name=null] - See the argument of the same name in {@linkcode List#setByName setByName}.
     * Only used if `i` is a number.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} The List after setting the `i`-th element to `x`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    set(i, x, { name = null, inPlace = false } = {}) {
        if (typeof i == "number") {
            return this.setByIndex(i, x, { name, inPlace });
        } else {
            return this.setByName(i, x, { inPlace });
        }
    }

    /**
     * @param {?Array} names - Array of strings of length equal to {@linkcode List#length length}.
     * This may contain duplicates.
     * Alternatively `null`, to remove existing names.
     * @param {Object} [options={}] - Further options.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} The List after replacing the names with `names`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setNames(names, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);

        if (names !== null) {
            if (names.length != this._values.length) {
                throw new Error("'names' and 'values' should have the same length");
            }
            for (const n of names) {
                if (typeof n != "string") {
                    throw new Error("'names' should be an array of strings");
                }
            }
        }

        target._names = new IndexedNames(names);
        return target;
    }

    /***********************************************/

    /**
     * @param {number} i - Index of the List element to delete.
     * This should be non-negative and no less than {@linkcode List#length length}.
     * @param {Object} [options={}] - Further options.
     * @param {?string} [options.name=null] - See the argument of the same name in {@linkcode List#setByName setByName}.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} The List after deleting the `i`-th element.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    deleteByIndex(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._values = target._values.slice();
        }

        this.#check_index(i);
        target._values.splice(i, 1);
        if (target._names !== null) {
            target._names = target._names.delete(i, { inPlace });
        }

        return target;
    }

    /**
     * @param {string} name - Name of the List element to delete.
     * This should already exist in {@linkcode List#names names}.
     * @param {?string} [options.name=null] - See the argument of the same name in {@linkcode List#setByName setByName}.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} The List after deleting the `name`d element.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    deleteByName(name, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._values = target._values.slice();
        }

        if (target._names == null) {
            throw new Error("no available names in this " + this.constructor.className);
        }

        // Don't cache as we're going to reset the lookup immediately, so it would be needlessly inefficient.
        let candidate = this._names.nameToIndexUncached(name);
        if (candidate < 0) {
            throw new Error("no matching name for '" + name + "' in this " + this.constructor.className);
        }

        target._values.splice(candidate, 1);
        if (target._names !== null) {
            target._names = target._names.delete(candidate, { inPlace });
        }

        return target;
    }

    /**
     * @param {string|number} i - Index or name of the List element to delete.
     * Numbers are passed to {@linkcode List#deleteByIndex deleteByIndex} and strings are passed to {@linkcode List#deleteByName deleteByName}.
     * @param {Object} [options={}] - Further options.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} The List after deleting the `i`-th element.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    delete(i, { inPlace = false } = {}) {
        if (typeof i == "number") {
            return this.deleteByIndex(i, { inPlace });
        } else {
            return this.deleteByName(i, { inPlace });
        }
    }

    /***********************************************/

    /**
     * @param {number} start - Index of the first element in the slice.
     * @param {number} end - Index past the last element in the slice.
     * @param {Object} [options={}] - Further options.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} A List that is sliced to `[start, end)`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    sliceRange(start, end, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);

        target._values = target._values.slice(start, end);
        if (this._names !== null) {
            target._names = new IndexedNames(target._names.names().slice(start, end));
        }

        return target;
    }

    /**
     * @param {Array} indices - Array of numbers or strings specifying the List elements to retain in the slice.
     * Numbers are interpreted as positional indices while strings are interpreted as names.
     * @param {Object} [options={}] - Further options.
     * @param {boolean} [options.inPlace=false] - Whether to modify this List instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {List} A List containing the specified elements in `indices`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    sliceIndices(indices, { inPlace = false } = {}) {
        let new_names = [];
        let new_values = [];

        for (let i of indices) {
            if (typeof i == "string") {
                if (this._names == null) {
                    throw new Error("no available names in this " + this.constructor.className);
                }
                i = this._names.nameToIndex(i);
            } else {
                this.#check_index(i);
            }

            new_values.push(this._values[i]);
            if (this._names !== null) {
                new_names.push(this._names.indexToName(i));
            }
        }

        let target = cutils.setterTarget(this, inPlace);
        target._values = new_values;
        if (this._names !== null) {
            target._names = new IndexedNames(new_names);
        }

        target._lookup = new Map;
        return target;
    }

    /***********************************************/

    /**
     * @return {iterator} An iterable iterator that can be used in, e.g., `for...of` constructs to loop over the List.
     * The list values are directly returned during iteration, i.e., names are ignored.
     */
    [Symbol.iterator]() {
        let counter = 0;
        let all_values = this._values;
        return {
            next: function() {
                if (counter < all_values.length) {
                    let val = all_values[counter];
                    counter++;
                    return { done: false, value: val };
                } else {
                    return { done: true };
                }
            },
            [Symbol.iterator]() {
                return this;
            },
        }
    }

    /***********************************************/

    _bioconductor_LENGTH() {
        return this.length();
    }

    _bioconductor_SLICE(i, { allowView = false }) {
        let sliced = this.sliceIndices(i);
        let output = new this.constructor;
        output._values = sliced._values;
        output._names = sliced._names;
        return output;
    }

    _bioconductor_CLONE({ deepCopy = true }) {
        let output = new this.constructor;
        output._values = cutils.cloneField(this._values, deepCopy);
        output._names = cutils.cloneField(this._names, deepCopy);
        return output;
    }

    _bioconductor_COMBINE(objects) {
        let all_values = this._values.slice();
        let all_names = null;
        if (this._names !== null) {
            all_names = this.names().slice();
        }

        for (let x of objects) {
            if (!(x instanceof List)) {
                x = new List(x);
            }

            const xvals = x.values();
            for (const y of xvals) {
                all_values.push(y);
            }

            const xnames = x.names();
            if (xnames === null) {
                if (all_names !== null) {
                    for (const y of xvals) {
                        all_names.push("");
                    }
                }
            } else {
                if (all_names === null) {
                    all_names = new Array(all_values.length - xvals.length).fill("");
                }
                for (const yn of xnames) {
                    all_names.push(yn);
                }
            }
        }

        let output = new this.constructor;
        output._values = all_values;
        output._names = new IndexedNames(all_names);
        return output;
    }
}

/**
 * Subclass of a {@linkplain List} that only contains integers or `null`s.
 * If a `null` is present, it should be treated as a missing value.
 * @extends List
 */
export class IntegerList extends List {
    #sanitize(x) {
        if (x !== null && !Number.isInteger(x)) {
            throw new Error("only integers or nulls can be stored in an IntegerList");
        }
        return x;
    }

    /**
     * @param {Array|Map|Object} values - Elements of the List.
     * This should only contain integers or `null`s.
     * @param {Object} [options={}] - Further options.
     */
    constructor(values, options = {}) {
        if (arguments.length == 0) {
            super();
        } else {
            super(values, options);
            for (const x of this._values) {
                this.#sanitize(x);
            }
        }
    }

    /**
     * @param {number} i - Index of the List element to set, see {@linkcode List#setByIndex List.setByIndex} for details.
     * @param {?number} x - Value of a List element as an integer or `null`.
     * @param {Object} [options={}] - Further options, see {@linkcode List#setByIndex List.setByIndex} for details.
     *
     * @return {List} The List after setting the `i`-th element to `x`, see {@linkcode List#setByIndex List.setByIndex} for details.
     */
    setByIndex(i, x, options = {}) {
        return super.setByIndex(i, this.#sanitize(x), options);
    }

    /**
     * @param {string} name - Name of the List element to set, see {@linkcode List#setByName List.setByName} for details.
     * @param {?number} x - Value of a List element as an integer or `null`.
     * @param {Object} [options={}] - Further options, see {@linkcode List#setByName List.setByName} for details.
     *
     * @return {List} The List after setting the `name`d entry to `x`, see {@linkcode List#setByName List.setByName} for details.
     */
    setByName(name, x, options = {}) {
        return super.setByName(name, this.#sanitize(x), options);
    }
}

/**
 * Subclass of a {@linkplain List} that only contains numbers or `null`s.
 * If a `null` is present, it should be treated as a missing value.
 * @extends List
 */
export class NumberList extends List {
    #sanitize(x) {
        if (x !== null && typeof x !== "number") {
            throw new Error("only numbers or nulls can be stored in a NumberList");
        }
        return x;
    }

    /**
     * @param {Array|Map|Object} values - Elements of the List.
     * This should only contain numbers or `null`s.
     * @param {Object} [options={}] - Further options.
     */
    constructor(values, options = {}) {
        if (arguments.length == 0) {
            super();
        } else {
            super(values, options);
            for (const x of this._values) {
                this.#sanitize(x);
            }
        }
    }

    /**
     * @param {number} i - Index of the List element to set, see {@linkcode List#setByIndex List.setByIndex} for details.
     * @param {?number} x - Value of a List element as a number or `null`.
     * @param {Object} [options={}] - Further options, see {@linkcode List#setByIndex List.setByIndex} for details.
     *
     * @return {List} The List after setting the `i`-th element to `x`, see {@linkcode List#setByIndex List.setByIndex} for details.
     */
    setByIndex(i, x, options = {}) {
        return super.setByIndex(i, this.#sanitize(x), options);
    }

    /**
     * @param {string} name - Name of the List element to set, see {@linkcode List#setByName List.setByName} for details.
     * @param {?number} x - Value of a List element as a number or `null`.
     * @param {Object} [options={}] - Further options, see {@linkcode List#setByName List.setByName} for details.
     *
     * @return {List} The List after setting the `name`d entry to `x`, see {@linkcode List#setByName List.setByName} for details.
     */
    setByName(name, x, options = {}) {
        return super.setByName(name, this.#sanitize(x), options);
    }
}

/**
 * Subclass of a {@linkplain List} that only contains strings or `null`s.
 * If a `null` is present, it should be treated as a missing value.
 * @extends List
 */
export class StringList extends List {
    #sanitize(x) {
        if (x !== null && typeof x !== "string") {
            throw new Error("only strings or nulls can be stored in a StringList");
        }
        return x;
    }

    /**
     * @param {Array|Map|Object} values - Elements of the List.
     * This should only contain strings or `null`s.
     * @param {Object} [options={}] - Further options.
     */
    constructor(values, options = {}) {
        if (arguments.length == 0) {
            super();
        } else {
            super(values, options);
            for (const x of this._values) {
                this.#sanitize(x);
            }
        }
    }

    /**
     * @param {number} i - Index of the List element to set, see {@linkcode List#setByIndex List.setByIndex} for details.
     * @param {?boolean} x - Value of a List element as a string or `null`.
     * @param {Object} [options={}] - Further options, see {@linkcode List#setByIndex List.setByIndex} for details.
     *
     * @return {List} The List after setting the `i`-th element to `x`, see {@linkcode List#setByIndex List.setByIndex} for details.
     */
    setByIndex(i, x, options = {}) {
        return super.setByIndex(i, this.#sanitize(x), options);
    }

    /**
     * @param {string} name - Name of the List element to set, see {@linkcode List#setByName List.setByName} for details.
     * @param {?boolean} x - Value of a List element as a string or `null`.
     * @param {Object} [options={}] - Further options, see {@linkcode List#setByName List.setByName} for details.
     *
     * @return {List} The List after setting the `name`d entry to `x`, see {@linkcode List#setByName List.setByName} for details.
     */
    setByName(name, x, options = {}) {
        return super.setByName(name, this.#sanitize(x), options);
    }
}

/**
 * Subclass of a {@linkplain List} that only contains booleans or `null`s.
 * If a `null` is present, it should be treated as a missing value.
 * @extends List
 */
export class BooleanList extends List {
    #sanitize(x) {
        if (x !== null && typeof x !== "boolean") {
            throw new Error("only booleans or nulls can be stored in a StringList");
        }
        return x;
    }

    /**
     * @param {Array|Map|Object} values - Elements of the List.
     * This should only contain booleans or `null`s.
     * @param {Object} [options={}] - Further options.
     */
    constructor(values, options = {}) {
        if (arguments.length == 0) {
            super();
        } else {
            super(values, options);
            for (const x of this._values) {
                this.#sanitize(x);
            }
        }
    }

    /**
     * @param {number} i - Index of the List element to set, see {@linkcode List#setByIndex List.setByIndex} for details.
     * @param {?boolean} x - Value of a List element as a boolean or `null`.
     * @param {Object} [options={}] - Further options, see {@linkcode List#setByIndex List.setByIndex} for details.
     *
     * @return {List} The List after setting the `i`-th element to `x`, see {@linkcode List#setByIndex List.setByIndex} for details.
     */
    setByIndex(i, x, options = {}) {
        return super.setByIndex(i, this.#sanitize(x), options);
    }

    /**
     * @param {string} name - Name of the List element to set, see {@linkcode List#setByName List.setByName} for details.
     * @param {?boolean} x - Value of a List element as a boolean or `null`.
     * @param {Object} [options={}] - Further options, see {@linkcode List#setByName List.setByName} for details.
     *
     * @return {List} The List after setting the `name`d entry to `x`, see {@linkcode List#setByName List.setByName} for details.
     */
    setByName(name, x, options = {}) {
        return super.setByName(name, this.#sanitize(x), options);
    }
}
