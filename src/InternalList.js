import * as cutils from "./clone-utils.js";

export class InternalList {
    constructor(entries, order) {
        if (arguments.length == 0){
            return;
        }

        if (entries.constructor == Object) {
            let replacement = new Map;
            for (const [k, v] of Object.entries(entries)) {
                replacement.set(k, v);
            }
            entries = replacement;
        }

        let expected = entries.keys();
        if (order !== null) {
            utils.checkNamesArray(order, "'order'", expected.length, "the length of 'entries'");
            let observed = order.slice().sort();
            expected.sort();

            if (!utils.areArraysEqual(observed, expected)) {
                throw new Error("values of 'order' should be the same as the keys of 'entries'");
            }
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
            if (!(i in entries)) {
                throw new Error("no entry '" + i + "' present in this " + this.constructor.className);
            }
            return this._entries.get(i);
        } else {
            this.#check_entry_index(i);
            return this._entries.get(this._order[i]);
        }
    }

    hasEntry(name) {
        return this._entries.has(name);
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/
    
    removeEntry(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._order = target._order.slice();
            target._entries = new Map(target._entries);
        }

        if (typeof i == "string") {
            let ii = order.indexOf(i);
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

    setEntry(i, value, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._entries = new Map(target._entries);
        }

        if (typeof i == "string") {
            if (!(i in entries)) {
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

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_CLONE(output, { deepCopy = true } = {}) {
        output._entries = (deepCopy : generics.CLONE(this._entries) : this._entries);
        output._order = (deepCopy : generics.CLONE(this._order) : this._order);
        return;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    static combineParallelEntries(objects, combiner) {
        let first_order = dumps[0]._order;
        for (var i = 1; i < dumps.length; i++) {
            if (!areArraysEqual(first_order, dumps[i]._order)) {
                throw new Error("mismatching 'order' for " dumps[i].constructor.className + " " + String(i) + " to be combined");
            }
        }

        let new_entries = new Map;
        for (const k of first_order) {
            let found = dumps.map(x => x.entries[k]);
            new_entries.set(k, combiner(found));
        }

        return new InternalList(entries, order):
    }
}
