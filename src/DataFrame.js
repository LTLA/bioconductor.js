import * as generics from "./AllGenerics.js";

export class DataFrame {
    constructor(x, { numberOfRows = null, rowNames = null, columnOrder = null }) {
        this._nrow = numberOfRows;
        this._row_names = rowNames;
        this._columns = x;
        this._column_order = columnOrder;

        let vals = Object.values(x);
        if (vals.length) {
            for (var i = 0; i < vals.length; i++) {
                let n = generics._length(vals[i]);
                if (this._nrow == null) {
                    this._nrow = n;
                } else if (n != this._nrow) {
                    throw new Error("expected all arrays in 'x' to have equal length");
                }
            }
        }

        if (columnOrder == null) {
            this._column_order = Object.keys(x);
        } else if (columnOrder.length != vals.length) {
            throw new Error("'columnOrder' should have the same length as 'x'");
        } else {
            let uniq = Array.from(new Set(columnOrder));
            uniq.sort();
            let expected = Object.keys(x);
            expected.sort();
            if (!areArraysEqual(uniq, expected)) {
                throw new Error("values of 'columnOrder' should be the same as the keys of 'x'");
            }
        }

        if (rowNames != null) {
            if (this._nrow == null) {
                this._nrow = rowNames.length;
            } else if (this._nrow != rowNames.length) {
                throw new Error("length of 'rowNames' is inconsistent with the number of rows of 'x'");
            }
        }

        if (this._nrow == null) {
            this._nrow = 0;
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    rowNames() {
        return this._row_names;
    }

    columnNames() {
        return this._column_order;
    }
    
    numberOfRows() {
        return this._nrow;
    }

    numberOfColumns() {
        return this._column_order.length;
    }

    column(i) {
        if (typeof i == "string") {
            return this._columns[i];
        } else {
            return this._columns[this._column_order[i]];
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    setColumn(i, value) {
        if (generics._length(value) != this._nrow) {
            throw new Error("expected 'value' to have the same length as the number of rows in 'x'");
        }

        if (typeof i == "string") {
            if (!(i in this._columns)) {
                this._column_order.push(i);
            }
            this._columns[i] = value;
            return;
        }

        if (i < 0 || i >= this.nrow()) {
            throw new Error("index 'i' is out of range of the number of columns");
        }
        this._columns[this._column_order[i]] = value;
        return;
    }

    setColumnNames(names) {
        if (names.length != this._column_order.length) {
            throw new Error("length of replacement 'names' must be equal to the number of columns");
        }

        let new_columns = {};
        for (var i = 0; i < names.length; i++) {
            if (names[i] in new_columns) {
                throw new Error("detected duplicates in replacement 'names'");
            }
            new_columns[names[i]] = this._columns[this._column_order[i]];
        }

        this._columns = new_columns;
        this._column_order = names.slice(); // make a copy.
        return;
    }

    sliceColumns(i) {
        // TODO.
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_length() {
        return this.numberOfRows();
    }

    _bioconductor_slice(i, { allowInPlace = false, allowView = false }) {
        let new_columns = {};
        for (const [k, v] of Object.entries(this._columns)) {
            new_columns[k] = generics._slice(v, i, { allowInPlace, allowView });
        }

        let new_row_names = (this._row_names == null ? null : generics._splice(this._row_names, i));

        let new_nrow;
        if (i.constructor == Object) {
            new_nrow = i.end - i.start;
            if (new_nrow < 0) {
                new_nrow *= -1;
            }
        } else {
            new_nrow = i.length;
        }

        if (allowInPlace) {
            this._columns = new_columns;
            this._row_names = new_row_names;
            this._nrow = new_nrow;
            return this;
        } else {
            let output = Object.create(this.prototype);
            output._columns = new_columns;
            output._row_names = new_row_names;
            output._column_order = this._column_order;
            output._nrow = new_nrow;
            return output;
        }
    }

    _bioconductor_combine(y, { allowAppend = false }) {
        let new_columns = {};
        for (const x of this._column_order) {
            let yarr = [this._columns[x]];
            for (const yi of y) {
                if (!(x in yi._columns)) {
                    throw new Error("missing column '" + x + "' across DataFrames to be combined");
                }
                yarr.push(yi._columns[x]);
            }
            
            let combined = generics._combine(yarr);
            new_columns[x] = combined;
        }

        let new_nrow = this._nrow;
        let has_row_names = this._row_names !== null;
        for (const yi of y) {
            if (yi.numberOfColumns() != this.numberOfColumns()) {
                throw new Error("mismatching number of columns across DataFrames to be combined");
            }

            if (!has_row_names && yi._row_names !== null) {
                has_row_names = true;
            }

            new_nrow += yi._nrow;
        }

        let new_row_names = null;
        if (has_row_names) {
            new_row_names = new Array(new_nrow);

            let counter = 0;
            let set_row_names = obj => {
                if (obj._row_names == null) {
                    new_row_names.fill(null, counter, counter + obj._nrow);
                } else {
                    new_row_names.set(object._row_names, counter);
                }
                counter += obj._nrow;
            }

            set_row_names(this);
            for (const yi of y) {
                set_row_names(yi);
            }
        }

        if (allowAppend) {
            this._columns = new_columns;
            this._nrow = new_nrow;
            this._row_names = new_row_names;
            return this;
        } else {
            let output = Object.create(this.prototype);
            output._columns = new_columns;
            output._row_names = new_row_names;
            output._column_order = this._column_order;
            output._nrow = new_nrow;
            return output;
        }
    }
};
