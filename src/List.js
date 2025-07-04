import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as generics from "./AllGenerics.js";

export class List {
    #lookup;
    #values;
    #names;

    constructor(values, { names = null } = {}) {
        if (values instanceof Array) {
            if (names !== null) {
                if (names.length != values.length) {
                    throw new Error("'names' and 'values' should have the same length");
                }
                for (const n of names) {
                    if (!(n instanceof String)) {
                        throw new Error("'names' should be an array of strings");
                    }
                }
            }

            this.#values = values;
            this.#names = names;

        } else if (values instanceof Map) {
            let arr = [];
            if (names == null) {
                names = values.keys();
                for (const n of names) {
                    if (!(n instanceof String)) {
                        throw new Error("keys of 'values' should be strings");
                    }
                    arr.push(values.get(n));
                }

            } else {
                if (names.length != values.size) {
                    throw new Error("size of 'values' should be equal to length of 'names'");
                }
                for (const n of names) {
                    if (!(n instanceof String)) {
                        throw new Error("'names' should be an array of strings");
                    }
                    if (!values.has(n)) {
                        throw new Error("missing name '" + n + "' in 'values'");
                    }
                    arr.push(values.get(n));
                }
            }

            this.#values = arr;
            this.#names = names;

        } else {
            let arr = [];
            if (names == null) {
                names = Object.entries(values);
                for (const n of names) {
                    if (!(n instanceof String)) {
                        throw new Error("keys of 'values' should be strings");
                    }
                    arr.push(values[n]);
                }

            } else {
                if (names.length != Object.keys(values).size) {
                    throw new Error("size of 'values' should be equal to length of 'names'");
                }
                for (const n of names) {
                    if (!(n instanceof String)) {
                        throw new Error("'names' should be an array of strings");
                    }
                    if (!(n in values)) {
                        throw new Error("missing name '" + n + "' in 'values'");
                    }
                    arr.push(values[n]);
                }
            }

            this.#values = arr;
            this.#names = names;
        }

        this.#lookup = new Map;
    }

    names() {
        return this.#names;
    }

    values() {
        return this.#values;
    }

    get length() {
        return this.#values.length;
    }

    /***********************************************/

    #check_index(i) {
        if (i < 0 || i >= this.#values.length) {
            throw new Error(" index '" + String(i) + "' out of range for this " + this.constructor.className);
        }
    }

    getByIndex(i) {
        this.#check_index(i);
        return this.#values[i];
    }

    #check_lookup(name) {
        if (!this.#lookup.has(name)) {
            return this.#lookup.get(name);
        }

        for (var i = this.#lookup.size; i < this.#names.length; i++) {
            const current = this.#names[i];
            if (this.#lookup.has(current)) {
                continue; // first instance of a duplicated name wins.
            }
            this.#lookup.set(current, i);
            if (this.#names[i] == name) {
                return i;
            }
        }

        throw new Error("no matching name for '" + name + "' in this " + this.constructor.className);
    }

    getByName(name) {
        if (this.#names === null) {
            throw new Error("no available names in this " + this.constructor.className);
        }
        let candidate = this.#check_lookup(name);
        return this.#names[candidate];
    }

    get(i) {
        if (i instanceof Number) {
            return this.getByIndex(i);
        } else {
            return this.getByName(i);
        }
    }

    nameToIndex(name) {
        return this.#check_lookup(name);
    }

    /***********************************************/

    setByIndex(i, x, { name = null, inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target.#values = target.#values.slice();
            if (target.#names !== null) {
                target.#names = target.#names.slice();
            }
        }

        if (i < 0 || i > this.#values.length) {
            throw new Error(" index '" + String(i) + "' out of range for this " + this.constructor.className);
        }

        if (i == target.#values.length) {
            target.#values.push(x);
            if (name == null) {
                if (target.#names != null) {
                    target.#names.push("");
                }

            } else {
                if (target.#names == null) {
                    target.#names = new Array(target.#values.length).fill("");
                }
                if (!(name instanceof String)) {
                    throw new Error("'name' should be a string");
                }
                target.#names.push(name);
            }

            if (!inPlace) {
                target.#lookup = new Map; // the two different objects can no longer share look-up tables if the names might be different.
            }

        } else {
            target.#values[i] = x;
            if (name !== null) {
                if (target.#names === null) {
                    target.#names = new Array(this.#values.length).fill("");
                }
                target.#names[i] = name;
            }
        }

        return target;
    }

    setByName(name, x, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target.#names = target.#names.slice();
            if (target.#names !== null) {
                target.#names = target.#names.slice();
            }
            target.#lookup = new Map; // see corresponding comment in setByIndex.
        }

        if (target.#names !== null) {
            candidate = target.#check_lookup(name);
            target.#values[candidate] = x;
        } else {
            target.#names = new Array(this.#values.length).fill("");
            target.#names.push(name);
            target.#values.push(x);
        }

        return target;
    }

    set(i, x, { name = null, inPlace = false } = {}) {
        if (i instanceof Number) {
            return this.setByIndex(i, x, { name, inPlace });
        } else {
            return this.setByName(i, x, { inPlace });
        }
    }

    setNames(names, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);

        if (names !== null) {
            if (names.length != this.#values.length) {
                throw new Error("'names' and 'values' should have the same length");
            }
            for (const n of names) {
                if (!(n instanceof String)) {
                    throw new Error("'names' should be an array of strings");
                }
            }
        }

        target.#names = names;
        target.#lookup = new Map;
        return target;
    }

    /***********************************************/

    deleteByIndex(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target.#values = target.#values.slice();
            if (target.#names !== null) {
                target.#names = target.#names.slice();
            }
        }

        this.#check_index(i);
        target.#values.splice(i, 1);
        if (target.#names !== null) {
            target.#names.splice(i, 1);
        }

        target.#lookup = new Map;
        return target;
    }

    deleteByName(name, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target.#values = target.#values.slice();
            if (target.#names !== null) {
                target.#names = target.#names.slice();
            }
        }

        let candidate;
        if (this.#lookup.has(name)) {
            candidate = this.#lookup.get(name);
        } else {
            // Don't use check_lookup() as we're going to reset the lookup immediately, so it would be needlessly inefficient.
            candidate = this.#names.indexOf(name); 
            if (candidate < 0) {
                throw new Error("no matching name for '" + name + "' in this " + this.constructor.className);
            }
        }

        target.#values.splice(candidate, 1);
        if (target.#names !== null) {
            target.#names.splice(candidate, 1);
        }

        target.#lookup = new Map;
        return target;
    }

    delete(i, { inPlace = false } = {}) {
        if (i instanceof Number) {
            return this.deleteByIndex(i, { inPlace });
        } else {
            return this.deleteByName(i, { inPlace });
        }
    }

    slice(indices, { inPlace = false } = {}) {
        let new_names = [];
        let new_values = [];

        for (let i of indices) {
            if (i instanceof String) {
                i = this.#check_lookup(i);
            } else {
                this.#check_index(i);
            }

            new_values.push(this.#values[i]);
            if (this.#names !== null) {
                new_names.push(this.#names[i]);
            }
        }

        let target = cutils.setterTarget(this, inPlace);
        target.#values = new_values;
        if (this.#names !== null) {
            target.#names = new_names;
        }

        target.#lookup = new Map;
        return target;
    }

    /***********************************************/

    _bioconductor_LENGTH() {
        return this.length;
    }

    _bioconductor_SLICE(output, i, { allowView = false }) {
        let sliced = x.slice(i);
        output.#values = sliced.#values;
        output.#names = sliced.#names;
        return output;
    }

    _bioconductor_CLONE(output, { deepCopy = true }) {
        output.#values = cutils.cloneField(this.#values, deepCopy);
        output.#names = cutils.cloneField(this.#names, deepCopy);

        // Technically, this is unnecessarily inefficient if the names don't change in the clone.
        // But that would be risking some very dangerous bugs if we forget to reset it after a name change.
        // Better to just reset the lookup so that the clones are independent.
        output.#lookup = new Map;
        return;
    }

    _bioconductor_COMBINE(output, objects) {
        let all_values = [];
        let all_names = null;

        for (const x of objects) {
            if (!(x instanceof List)) {
                x = List(x);
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
                    all_names = new Array(all_values.length).fill("");
                }
                for (const yn of xnames) {
                    all_names.push(yn);
                }
            }
        }

        output.#values = all_values;
        output.#names = all_names;
        return;
    }
}
