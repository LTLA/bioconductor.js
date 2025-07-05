import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as generics from "./AllGenerics.js";

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
 */
export class List {
    _lookup;
    _values;
    _names;

    /**
     * @param {Array|Map|Object} values - Elements of the List.
     * For Maps or objects, the values (in order of iteration) are used as the List elements.
     * @param {Object} [options={}] - Further options.
     * @param {?Array} [options.names=null] - An array of strings containing the names of the List elements.
     * If provided, this should be of the same length as `values`.
     * If `values` is a Map or object, `names` should have the same keys.
     * If `values` is an array, the names may contain duplicate strings.
     * If `null` and `values` is an array, the List will be unnamed.
     */
    constructor(values, { names = null } = {}) {
        if (arguments.length == 0) {
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
            this._names = names;

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
            this._names = names;
        }

        this._lookup = new Map;
    }

    /**
     * @return {?Array} Array of names of the List elements, or `null` if the List is unnamed.
     */
    names() {
        return this._names;
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

    #check_lookup(name, error) {
        if (this._lookup.has(name)) {
            return this._lookup.get(name);
        }

        for (var i = this._lookup.size; i < this._names.length; i++) {
            const current = this._names[i];
            if (this._lookup.has(current)) {
                continue; // first instance of a duplicated name wins.
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
        let candidate = this.#check_lookup(name, true);
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
     * @return {number} Index of the name in {@linkcode List#names names}.
     * If duplicate names are present, the first occurrence is returned.
     */
    nameToIndex(name) {
        return this.#check_lookup(name, true);
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
            if (target._names !== null) {
                target._names = target._names.slice();
            }
        }

        if (i < 0 || i > this._values.length) {
            throw new Error(" index '" + String(i) + "' out of range for this " + this.constructor.className);
        }

        if (i == target._values.length) {
            target._values.push(x);
            if (name == null) {
                if (target._names != null) {
                    target._names.push("");
                }

            } else {
                if (target._names == null) {
                    target._names = new Array(target._values.length).fill("");
                }
                if (typeof name != "string") {
                    throw new Error("'name' should be a string");
                }
                target._names[target._values.length - 1] = name;
            }

        } else {
            target._values[i] = x;
            if (name !== null) {
                if (target._names === null) {
                    target._names = new Array(this._values.length).fill("");
                }
                target._names[i] = name;
                if (inPlace) {
                    target._lookup = new Map; // lookup is invalidated if the existing names change.
                }
            }
        }

        return target;
    }

    /**
     * @param {number} name - Name of the List element to set.
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
            if (target._names !== null) {
                target._names = target._names.slice();
            }
        }

        if (target._names !== null) {
            let candidate = target.#check_lookup(name, false);
            if (candidate < 0) {
                target._values.push(x);
                target._names.push(name);
            } else {
                target._values[candidate] = x;
            }
        } else {
            target._names = new Array(this._values.length).fill("");
            target._names.push(name);
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

        target._names = names;
        if (inPlace) {
            target._lookup = new Map; // lookup is invalidated if the existing names change.
        }

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
            if (target._names !== null) {
                target._names = target._names.slice();
            }
        }

        this.#check_index(i);
        target._values.splice(i, 1);
        if (target._names !== null) {
            target._names.splice(i, 1);
            if (inPlace) {
                target._lookup = new Map; // lookup is invalidated if existing names change.
            }
        }

        return target;
    }

    /**
     * @param {number} name - Name of the List element to delete.
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
            if (target._names !== null) {
                target._names = target._names.slice();
            }
        }

        if (target._names == null) {
            throw new Error("no available names in this " + this.constructor.className);
        }

        let candidate;
        if (this._lookup.has(name)) {
            candidate = this._lookup.get(name);
        } else {
            // Don't use check_lookup() as we're going to reset the lookup immediately, so it would be needlessly inefficient.
            candidate = this._names.indexOf(name); 
            if (candidate < 0) {
                throw new Error("no matching name for '" + name + "' in this " + this.constructor.className);
            }
        }

        target._values.splice(candidate, 1);
        if (target._names !== null) {
            target._names.splice(candidate, 1);
            if (inPlace) {
                target._lookup = new Map; // lookup is invalidated if existing names change.
            }
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
            target._names = target._names.slice(start, end);
            if (inPlace) {
                target._lookup = new Map;
            }
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
                i = this.#check_lookup(i);
            } else {
                this.#check_index(i);
            }

            new_values.push(this._values[i]);
            if (this._names !== null) {
                new_names.push(this._names[i]);
            }
        }

        let target = cutils.setterTarget(this, inPlace);
        target._values = new_values;
        if (this._names !== null) {
            target._names = new_names;
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

    _bioconductor_SLICE(output, i, { allowView = false }) {
        let sliced = this.sliceIndices(i);
        output._values = sliced._values;
        output._names = sliced._names;
        output._lookup = new Map;
        return output;
    }

    _bioconductor_CLONE(output, { deepCopy = true }) {
        output._values = cutils.cloneField(this._values, deepCopy);
        output._names = cutils.cloneField(this._names, deepCopy);

        // Technically, this is unnecessarily inefficient if the names don't change in the clone.
        // But that would be risking some very dangerous bugs if we forget to reset it after a name change.
        // Better to just reset the lookup so that the clones are independent.
        output._lookup = new Map;
        return;
    }

    _bioconductor_COMBINE(output, objects) {
        let all_values = [];
        let all_names = null;

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

        output._values = all_values;
        output._names = all_names;
        output._lookup = new Map;
        return;
    }
}
