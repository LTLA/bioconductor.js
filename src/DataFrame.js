import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";

export class DataFrame {
    constructor(columns, { numberOfRows = null, rowNames = null, columnOrder = null } = {}) {
        this._numberOfRows = numberOfRows;
        this._rowNames = rowNames;
        this._columns = columns;
        this._columnOrder = columnOrder;

        let vals = Object.values(columns);
        if (vals.length) {
            for (var i = 0; i < vals.length; i++) {
                let n = generics.LENGTH(vals[i]);
                if (this._numberOfRows == null) {
                    this._numberOfRows = n;
                } else if (n != this._numberOfRows) {
                    throw new Error("expected all arrays in 'x' to have equal length");
                }
            }
        }

        if (columnOrder == null) {
            this._columnOrder = Object.keys(columns);
        } else if (columnOrder.length != vals.length) {
            throw new Error("'columnOrder' should have the same length as 'x'");
        } else {
            let uniq = Array.from(new Set(columnOrder));
            uniq.sort();
            let expected = Object.keys(columns);
            expected.sort();
            if (!utils.areArraysEqual(uniq, expected)) {
                throw new Error("values of 'columnOrder' should be the same as the keys of 'x'");
            }
        }

        if (rowNames != null) {
            if (this._numberOfRows == null) {
                this._numberOfRows = rowNames.length;
            } else if (this._numberOfRows != rowNames.length) {
                throw new Error("length of 'rowNames' is inconsistent with the number of rows of 'x'");
            }
        }

        if (this._numberOfRows == null) {
            this._numberOfRows = 0;
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _check_index(i) {
        if (i < 0 || i >= this._columnOrder.length) {
            throw new Error("column index '" + String(i) + "' out of range for this DataFrame");
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    rowNames() {
        return this._rowNames;
    }

    columnNames() {
        return this._columnOrder;
    }

    hasColumn(name) {
        return name in this._columns;
    }

    numberOfRows() {
        return this._numberOfRows;
    }

    numberOfColumns() {
        return this._columnOrder.length;
    }

    column(i) {
        if (typeof i == "string") {
            if (!(i in this._columns)) {
                throw new Error("no column '" + i + "' present in this DataFrame");
            }
            return this._columns[i];
        } 

        this._check_index(i);
        return this._columns[this._columnOrder[i]];
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    $removeColumn(i) {
        if (typeof i == "string") {
            let ii = this._columnOrder.indexOf(i);
            if (ii < 0) {
                throw new Error("no column '" + i + "' in this DataFrame");
            }
            this._columnOrder.splice(ii, 1);
            delete this._columns[i];
        } else {
            this._check_index(i);
            let n = this._columnOrder[i];
            this._columnOrder.splice(i, 1);
            delete this._columns[n];
        }
    }

    $setColumn(i, value) {
        if (generics.LENGTH(value) != this._numberOfRows) {
            throw new Error("expected 'value' to have the same length as the number of rows in 'x'");
        }

        if (typeof i == "string") {
            if (!(i in this._columns)) {
                this._columnOrder.push(i);
            }
            this._columns[i] = value;
            return;
        }

        this._check_index(i);
        this._columns[this._columnOrder[i]] = value;
        return;
    }

    $setColumnNames(names) {
        if (names.length != this._columnOrder.length) {
            throw new Error("length of replacement 'names' must be equal to the number of columns");
        }

        let new_columns = {};
        for (var i = 0; i < names.length; i++) {
            if (names[i] in new_columns) {
                throw new Error("detected duplicates in replacement 'names'");
            }
            new_columns[names[i]] = this._columns[this._columnOrder[i]];
        }

        this._columns = new_columns;
        this._columnOrder = names;
        return;
    }

    $setRowNames(names) {
        if (names != null) {
            if (names.length != this._numberOfRows) {
                throw new Error("length of replacement 'names' must be equal to the number of rows");
            }
        }
        this._rowNames = names;
        return;
    }

    $sliceColumns(i) {
        let new_columns = {};
        let new_columnOrder = [];

        for (var ii of i) {
            if (typeof ii != "string") {
                this._check_index(ii);
                ii = this._columnOrder[ii];
            }
            if (ii in new_columns) {
                throw new Error("duplicate columns detected in slice request");
            }

            new_columns[ii] = this._columns[ii];
            new_columnOrder.push(ii);
        }

        this._columns = new_columns;
        this._columnOrder = new_columnOrder;
        return;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_LENGTH() {
        return this.numberOfRows();
    }

    _bioconductor_SLICE(i, { allowInPlace = false, allowView = false }) {
        let options = { allowInPlace, allowView };
        let new_columns = {};
        for (const [k, v] of Object.entries(this._columns)) {
            new_columns[k] = generics.SLICE(v, i, options);
        }

        let new_rowNames = (this._rowNames == null ? null : generics.SLICE(this._rowNames, i, options));

        let new_numberOfRows;
        if (i.constructor == Object) {
            new_numberOfRows = i.end - i.start;
        } else {
            new_numberOfRows = i.length;
        }

        if (allowInPlace) {
            this._columns = new_columns;
            this._rowNames = new_rowNames;
            this._numberOfRows = new_numberOfRows;
            return this;
        } else {
            let output = Object.create(this.constructor.prototype);
            output._columns = new_columns;
            output._rowNames = new_rowNames;
            output._columnOrder = this._columnOrder;
            output._numberOfRows = new_numberOfRows;
            return output;
        }
    }

    _bioconductor_COMBINE(y, { allowAppend = false }) {
        let options = { allowAppend };
        let new_columns = {};
        for (const x of this._columnOrder) {
            let yarr = [];
            for (const yi of y) {
                if (!(x in yi._columns)) {
                    throw new Error("missing column '" + x + "' across DataFrames to be combined");
                }
                yarr.push(yi._columns[x]);
            }
            new_columns[x] = generics.COMBINE(this._columns[x], yarr, options);
        }

        let new_numberOfRows = this._numberOfRows;
        let has_rowNames = this._rowNames !== null;
        for (const yi of y) {
            if (yi.numberOfColumns() != this.numberOfColumns()) {
                throw new Error("mismatching number of columns across DataFrames to be combined");
            }

            if (!has_rowNames && yi._rowNames !== null) {
                has_rowNames = true;
            }

            new_numberOfRows += yi._numberOfRows;
        }

        let new_rowNames = null;
        if (has_rowNames) {
            new_rowNames = new Array(new_numberOfRows);

            let counter = 0;
            let set_rowNames = obj => {
                if (obj._rowNames == null) {
                    new_rowNames.fill(null, counter, counter + obj._numberOfRows);
                    counter += obj._numberOfRows;
                } else {
                    obj._rowNames.forEach(x => {
                        new_rowNames[counter] = x;
                        counter++;
                    });
                }
            }

            set_rowNames(this);
            for (const yi of y) {
                set_rowNames(yi);
            }
        }

        if (allowAppend) {
            this._columns = new_columns;
            this._numberOfRows = new_numberOfRows;
            this._rowNames = new_rowNames;
            return this;
        } else {
            let output = Object.create(this.constructor.prototype); // avoid validity checks.
            output._columns = new_columns;
            output._rowNames = new_rowNames;
            output._columnOrder = this._columnOrder;
            output._numberOfRows = new_numberOfRows;
            return output;
        }
    }

    _bioconductor_CLONE({ deepCopy = true }) {
        let new_columnOrder = this._columnOrder.slice();
        let new_rowNames = (this._rowNames == null ? null : this._rowNames.slice());

        let new_columns = { ...(this._columns) };
        if (deepCopy) {
            for (const [k, v] of Object.entries(new_columns)) {
                new_columns[k] = generics.CLONE(v, { deepCopy });
            }
        }

        let output = Object.create(this.constructor.prototype); // avoid validity checks.
        output._columns = new_columns;
        output._rowNames = new_rowNames;
        output._columnOrder = new_columnOrder;
        output._numberOfRows = this._numberOfRows;
        return output;
    }
};
