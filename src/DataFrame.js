import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";

/**
 * A DataFrame is a collection of equilength vector-like objects as "columns".
 * The number of rows in the DataFrame is equal to the length of the columns, where the i-th row consists of the i-th element from each column.
 *
 * This class supports optional row names, which are either `null` or an array of strings of length equal to the number of rows.
 *
 * This class supports empty instances with a non-zero number of rows, which may be useful for piece-wise construction.
 *
 * The vector-like object for each column is expected to have methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 *
 * The DataFrame itself defines methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 */
export class DataFrame {
    /**
     * @param {Object} columns - Object where keys are the column names and the values are equilength vector-like objects.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?number} [options.numberOfRows=null] - Non-negative value specifying the number of rows in the DataFrame.
     * If `null`, this is automatically determined from the length of the vectors in `columns`, or from the length of `rowNames`.
     * If non-`null`, this should not conflict with the inferred lengths from `columns` or `rowNames`.
     * @param {?Array} [options.rowNames=null] - Array of strings containing the names for each row.
     * If non-`null`, this should have the same length as the vectors inside `columns`, if any exist.
     * If `null`, no row names are used.
     * @param {?Array} [options.columnOrder=null] - Array of strings specifying the ordering of the columns.
     * If non-`null`, this should have the same values as the keys of `columns`.
     * If `null`, an arbitrary ordering is obtained from `columns`.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(columns, { numberOfRows = null, rowNames = null, columnOrder = null, metadata = {} } = {}) {
        this._numberOfRows = numberOfRows;
        this._rowNames = rowNames;
        this._columns = columns;
        this._columnOrder = columnOrder;
        this._metadata = metadata;

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
            DataFrame._check_names(columnOrder, "column names");
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
            DataFrame._check_names(rowNames, "row names");
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

    static _check_names(names, msg) {
        for (const x of names) {
            if (typeof x !== "string") {
                throw new Error("array of " + msg + " should only contain strings");
            }
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {?Array} Array of strings containing row names, or `null` if no row names are available.
     */
    rowNames() {
        return this._rowNames;
    }

    /**
     * @return {Array} Array of strings containing the column names in the specified order.
     */
    columnNames() {
        return this._columnOrder;
    }

    /**
     * @param {string} name - Name of a column.
     * @return {boolean} Whether the column exists in this DataFrame.
     */
    hasColumn(name) {
        return name in this._columns;
    }

    /**
     * @return {number} Number of rows in this DataFrame.
     */
    numberOfRows() {
        return this._numberOfRows;
    }

    /**
     * @return {number} Number of columns in this DataFrame.
     */
    numberOfColumns() {
        return this._columnOrder.length;
    }

    /**
     * @param {string|number} i - Column to retrieve, either by name or index.
     * @return {*} The contents of column `i` as a vector-like object.
     */
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

    /**
     * @return {Object} Object containing arbitrary metadata.
     */
    metadata() {
        return this._metadata;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {string|number} i - Column to remove, either by name or index.
     * @return {DataFrame} Reference to this DataFrame after the column is removed.
     */
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
        return this;
    }

    /**
     * @param {string|number} i - Column to add, either by name or index.
     * Numeric `i` should be non-negative and less than the number of columns.
     * @return {DataFrame} Reference to this DataFrame with modified columns.
     * - If `i` is a number, the column at the specified index is replaced.
     * - If `i` is a string, any column with the same name is replaced.
     *   If no such column exists, a new column is appended to the DataFrame.
     */
    $setColumn(i, value) {
        if (generics.LENGTH(value) != this._numberOfRows) {
            throw new Error("expected 'value' to have the same length as the number of rows in 'x'");
        }

        if (typeof i == "string") {
            if (!(i in this._columns)) {
                this._columnOrder.push(i);
            }
            this._columns[i] = value;
        } else {
            this._check_index(i);
            this._columns[this._columnOrder[i]] = value;
        }
        return this;
    }

    /**
     * @param {Array} names - Array of unique strings containing the new name for each column.
     * This should have the same length as {@linkcode DataFrame#columnNames DataFrame.columnNames}.
     * @return {DataFrame} Reference to this DataFrame with modified column names.
     */
    $setColumnNames(names) {
        DataFrame._check_names(names, "column names");
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
        return this;
    }

    /**
     * @param {?Array} names - Array of unique strings containing the new name for each row.
     * This should have the same length as {@linkcode DataFrame#numberOfRows DataFrame.numberOfRows}.
     *
     * Alternatively, this may be `null` to remove any existing column names.
     * @return {DataFrame} Reference to this DataFrame with modified row names.
     */
    $setRowNames(names) {
        if (names != null) {
            if (names.length != this._numberOfRows) {
                throw new Error("length of replacement 'names' must be equal to the number of rows");
            }
            DataFrame._check_names(names, "row names");
        }
        this._rowNames = names;
        return this;
    }

    /**
     * @param {Array} i - Array of strings or indices specifying the columns to retain in the slice.
     * This should refer to unique column names.
     *
     * @return {DataFrame} Reference to this DataFrame after slicing to the specified columns.
     */
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
        return this;
    }

    /**
     * @param {Object} meta - Object containing the metadata.
     *
     * @return {DataFrame} Reference to this DataFrame after replacing the metadata.
     */
    $setMetadata(meta) {
        this._metadata = meta;
        return this;
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
            output._metadata = this._metadata;
            return output;
        }
    }

    _bioconductor_COMBINE(objects, { allowAppend = false }) {
        let options = { allowAppend };
        let new_columns = {};
        for (const x of this._columnOrder) {
            let yarr = [];
            for (const yi of objects) {
                if (!(x in yi._columns)) {
                    throw new Error("missing column '" + x + "' across DataFrames to be combined");
                }
                yarr.push(yi._columns[x]);
            }
            new_columns[x] = generics.COMBINE(yarr, options);
        }

        let new_numberOfRows = 0;
        let has_rowNames = false;
        for (const yi of objects) {
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
            for (const obj of objects) {
                if (obj._rowNames == null) {
                    new_rowNames.fill("", counter, counter + obj._numberOfRows);
                    counter += obj._numberOfRows;
                } else {
                    obj._rowNames.forEach(x => {
                        new_rowNames[counter] = x;
                        counter++;
                    });
                }
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
            output._metadata = this._metadata;
            return output;
        }
    }

    _bioconductor_CLONE({ deepCopy = true }) {
        let new_columnOrder = this._columnOrder.slice();
        let new_rowNames = (this._rowNames == null ? null : this._rowNames.slice());
        let new_columns = generics.CLONE(this._columns, { deepCopy });
        let new_meta = generics.CLONE(this._metadata, { deepCopy });

        let output = Object.create(this.constructor.prototype); // avoid validity checks.
        output._columns = new_columns;
        output._rowNames = new_rowNames;
        output._columnOrder = new_columnOrder;
        output._numberOfRows = this._numberOfRows;
        output._metadata = new_meta; 
        return output;
    }
};

/**
 * Flexibly combine multiple DataFrames by row by filling in missing columns with an array of `null`s.
 * This is equivalent to calling {@linkcode COMBINE} on an array of DataFrames that may have mismatching columns.
 *
 * @param {Array} objects - Array of {@linkplain DataFrame}s to be combined.
 *
 * @return {DataFrame} The combined DataFrame, where the number of rows is equal to sum of rows across `objects`,
 * and the columns is equal to the union of columns across `objects`.
 */
export function flexibleCombineRows(objects) {
    let ckeys = new Set();
    let corder = [];
    for (const current of objects) {
        let cnames = current.columnNames();
        for (const a of cnames) {
            if (!ckeys.has(a)) {
                ckeys.add(a);
                corder.push(a);
            }
        }
    }

    let copies = [];
    for (const current of objects) {
        let dummy = new Array(current.numberOfRows());
        dummy.fill(null);
        let copy = generics.CLONE(current, { deepCopy: false });

        for (const a of corder) {
            if (!current.hasColumn(a)) {
                copy.$setColumn(a, dummy);
            }
        }

        copies.push(copy);
    }

    return generics.COMBINE(copies);
}
