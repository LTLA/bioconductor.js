import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as generics from "./AllGenerics.js";

export class InternalList {
    constructor(entries, order) {
        if (arguments.length == 0){
            return;
        }

        entries = utils.object2map(entries);

        let expected = Array.from(entries.keys());
        if (order !== null) {
            utils.checkNamesArray(order, "'order'", expected.length, "the length of 'entries'");
            let observed = order.slice().sort();
            expected.sort();

            if (!utils.areArraysEqual(observed, expected)) {
                throw new Error("values of 'order' should be the same as the keys of 'entries'");
            }
        } else {
            order = expected;
        }

        this._entries = entries;
        this._order = order;
    }

    static className = "InternalList";

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    names() {
        return this._order;
    }

    numberOfEntries() {
        return this._order.length;
    }

    #check_entry_index(i) {
        if (i < 0 || i >= this._order.length) {
            throw new Error(" index '" + String(i) + "' out of range for this " + this.constructor.className);
        }
    }

    entry(i) {
        if (typeof i == "string") {
            if (!this._entries.has(i)) {
                throw new Error("no entry '" + i + "' present in this " + this.constructor.className);
            }
            return this._entries.get(i);
        } else {
            this.#check_entry_index(i);
            return this._entries.get(this._order[i]);
        }
    }

    has(name) {
        return this._entries.has(name);
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/
    
    delete(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            // Shallow copies so that we can do our setting.
            target._order = target._order.slice();
            target._entries = new Map(target._entries); 
        }

        if (typeof i == "string") {
            let ii = target._order.indexOf(i);
            if (ii < 0) {
                throw new Error("no entry '" + i + "' present in this " + this.constructor.className);
            }
            target._order.splice(ii, 1); 
            target._entries.delete(i);
        } else {
            this.#check_entry_index(i);
            let n = target._order[i];
            target._order.splice(i, 1);
            target._entries.delete(n);
        }

        return target;
    }

    set(i, value, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            // Shallow copy so that we can do our setting.
            target._entries = new Map(target._entries);
        }

        if (typeof i == "string") {
            if (!target._entries.has(i)) {
                if (!inPlace) {
                    target._order = target._order.slice();
                }
                target._order.push(i);
            }
            target._entries.set(i, value);
        } else {
            this.#check_entry_index(i);
            target._entries.set(target._order[i], value);
        }

        return target;
    }

    setNames(names, { inPlace = false } = {}) {
        utils.checkNamesArray(names, "replacement 'names'", this._order.length, "length of 'names()'");

        let new_entries = new Map;
        for (var i = 0; i < names.length; i++) {
            if (new_entries.has(names[i])) {
                throw new Error("detected duplicate value '" + names[i] + "' in replacement 'names'");
            }
            new_entries.set(names[i], this._entries.get(this._order[i]));
        }

        let target = cutils.setterTarget(this, inPlace);
        target._entries = new_entries;
        target._order = names;
        return target;
    }

    slice(indices, { inPlace = false } = {}) {
        let new_entries = new Map;
        let new_order = [];

        for (var ii of indices) {
            if (typeof ii != "string") {
                this.#check_entry_index(ii);
                ii = this._order[ii];
            }
            if (new_entries.has(ii)) {
                throw new Error("duplicate entries detected in slice request");
            } else if (!this._entries.has(ii)) {
                throw new Error("slice contains missing entry '" + ii + "' ");
            }

            new_entries.set(ii, this._entries.get(ii));
            new_order.push(ii);
        }

        let target = cutils.setterTarget(this, inPlace);
        target._entries = new_entries;
        target._order = new_order;
        return target;
    }

    reorder(indices, { inPlace = false } = {}) {
        // Reorder can be slightly more efficient than slice because we just
        // need to change the ordering vector rather than creating a new Map.
        if (indices.length !== this._order.length) {
            throw utils.formatLengthError("reordered indices", "the number of existing entries");
        }

        let new_order = [];
        for (var ii of indices) {
            if (typeof ii != "string") {
                this.#check_entry_index(ii);
                ii = this._order[ii];
            }
            if (!this._entries.has(ii)) {
                throw new Error("missing entry '" + ii + "' among the reordered indices");
            }
            new_order.push(ii);
        }

        let target = cutils.setterTarget(this, inPlace);
        target._order = new_order;
        return target;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_CLONE(output, { deepCopy = true } = {}) {
        output._entries = (deepCopy ? generics.CLONE(this._entries) : this._entries);
        output._order = (deepCopy ? generics.CLONE(this._order) : this._order);
        return;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    apply(FUN, { inPlace = false } = {}) {
        let new_entries = (inPlace ? this._entries : new Map);
        for (const [k, v] of this._entries) {
            new_entries.set(k, FUN(v));
        }
        return (inPlace ? this : new InternalList(new_entries, this._order));
    }

    static parallelCombine(objects, combiner) {
        let first_order = objects[0]._order;
        for (var i = 1; i < objects.length; i++) {
            if (!utils.areArraysEqual(first_order, objects[i]._order)) {
                throw new Error("detected differences in names between first object and object " + String(i) + " to be combined");
            }
        }

        let new_entries = new Map;
        for (const k of first_order) {
            let found = objects.map(x => x._entries.get(k));
            new_entries.set(k, combiner(found));
        }

        return new InternalList(new_entries, first_order);
    }
}
