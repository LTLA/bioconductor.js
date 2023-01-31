import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as ann from "./Annotated.js";

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
 *
 * @augments Annotated
 */
export class DataFrame extends ann.Annotated {
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
        if (arguments.length == 0) {
            super();
            return;
        }

        super(metadata);
        this._numberOfRows = numberOfRows;
        this._rowNames = rowNames;

        let vals = Object.values(columns);
        for (var i = 0; i < vals.length; i++) {
            let n = generics.LENGTH(vals[i]);
            if (this._numberOfRows == null) {
                this._numberOfRows = n;
            } else if (n != this._numberOfRows) {
                throw new Error("expected all arrays in 'columns' to have equal length");
            }
        }
        this._columns = {
            entries: columns,
            order: utils.checkEntryOrder(columns, columnOrder, "column")
        };

        if (rowNames != null) {
            if (this._numberOfRows == null) {
                this._numberOfRows = rowNames.length;
            }
            utils.checkNamesArray(rowNames, "'rowNames'", this._numberOfRows, "'numberOfRows' or the length of arrays in 'columns'");
        }

        if (this._numberOfRows == null) {
            this._numberOfRows = 0;
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
        return this._columns.order;
    }

    /**
     * @param {string} name - Name of a column.
     * @return {boolean} Whether the column exists in this DataFrame.
     */
    hasColumn(name) {
        return name in this._columns.entries;
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
        return this._columns.order.length;
    }

    /**
     * @param {string|number} i - Column to retrieve, either by name or index.
     * @return {*} The contents of column `i` as a vector-like object.
     */
    column(i) {
        return utils.retrieveSingleEntry(this._columns.entries, this._columns.order, i, "column", "DataFrame");
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {string|number} i - Column to remove, either by name or index.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DataFrame instance in place.
     * If `false`, a new instance is returned.
     * 
     * @return {DataFrame} The DataFrame after removing the specified column.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    removeColumn(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (inPlace) {
            target._columns = cutils.shallowCloneEntries(target._columns);
        }

        utils.removeSingleEntry(target._columns.entries, target._columns.order, i, "column", "SummarizedExperiment");
        return target;
    }

    /**
     * @param {string|number} i - Column to remove, either by name or index.
     * @return {DataFrame} A reference to this DataFrame after removing the specified column.
     */
    $removeColumn(i) {
        return this.removeColumn(i, { inPlace: true });
    }

    /**
     * @param {string|number} i - Identity of the column to add, either by name or index.
     * - If `i` is a number, the column at the specified index is replaced.
     *   `i` should be non-negative and less than the number of columns.
     * - If `i` is a string, any column with the same name is replaced.
     *   If no such column exists, a new column is appended to the DataFrame.
     * @param {*} value - Array-like column to set/add as the column.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DataFrame instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {DataFrame} The DataFrame after adding/replacing the specified column.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setColumn(i, value, { inPlace = false } = {}) {
        if (generics.LENGTH(value) != this._numberOfRows) {
            throw new Error("expected 'value' to have the same length as the number of rows in 'x'");
        }

        let target = cutils.setterTarget(this, inPlace);
        if (inPlace) {
            target._columns = cutils.shallowCloneEntries(target._columns);
        }

        utils.setSingleEntry(target._columns.entries, target._columns.order, i, value, "column", "SummarizedExperiment");
        return target;
    }

    /**
     * @param {string|number} i - Identity of the column to add, either by name or index.
     * - If `i` is a number, the column at the specified index is replaced.
     *   `i` should be non-negative and less than the number of columns.
     * - If `i` is a string, any column with the same name is replaced.
     *   If no such column exists, a new column is appended to the DataFrame.
     * @param {*} value - Array-like column to set/add as the column.
     *
     * @return {DataFrame} A reference to this DataFrame after adding/replacing the specified column.
     */
    $setColumn(i, value) {
        return this.setColumn(i, value, { inPlace: true });
    }

    /**
     * @param {Array} names - Array of unique strings containing the new name for each column.
     * This should have the same length as {@linkcode DataFrame#columnNames DataFrame.columnNames}.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DataFrame instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {DataFrame} The DataFrame with modified column names.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setColumnNames(names, { inPlace = false } = {}) {
        utils.checkNamesArray(names, "replacement 'names'", this._columns.order.length, "'numberOfColumns()'");

        let new_columns = {};
        for (var i = 0; i < names.length; i++) {
            if (names[i] in new_columns) {
                throw new Error("detected duplicates in replacement 'names'");
            }
            new_columns[names[i]] = this._columns.entries[this._columns.order[i]];
        }

        let target = cutils.setterTarget(this, inPlace);
        target._columns = { entries: new_columns, order: names };
        return target;
    }

    /**
     * @param {Array} names - Array of unique strings containing the new name for each column.
     * This should have the same length as {@linkcode DataFrame#columnNames DataFrame.columnNames}.
     * @return {DataFrame} A reference to this DataFrame with modified column names.
     */
    $setColumnNames(names) {
        return this.setColumnNames(names, { inPlace: true });
    }

    /**
     * @param {?Array} names - Array of unique strings containing the new name for each row.
     * This should have the same length as {@linkcode DataFrame#numberOfRows DataFrame.numberOfRows}.
     *
     * Alternatively, this may be `null` to remove any existing column names.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DataFrame instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {DataFrame} The DataFrame with modified row names.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setRowNames(names, { inPlace = false } = {}) {
        if (names != null) {
            utils.checkNamesArray(names, "replacement 'names'", this._numberOfRows, "'numberOfRows()'");
        }

        let target = cutils.setterTarget(this, inPlace);
        target._rowNames = names;
        return target;
    }

    /**
     * @param {?Array} names - Array of unique strings containing the new name for each row.
     * This should have the same length as {@linkcode DataFrame#numberOfRows DataFrame.numberOfRows}.
     *
     * Alternatively, this may be `null` to remove any existing column names.
     * @return {DataFrame} A reference to this DataFrame with modified row names.
     */
    $setRowNames(names) {
        return this.setRowNames(names, { inPlace: true });
    }

    /**
     * @param {Array} i - Array of strings or indices specifying the columns to retain in the slice.
     * This should refer to unique column names.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DataFrame instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {DataFrame} Reference to this DataFrame after slicing to the specified columns.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    sliceColumns(i, { inPlace = false } = {}) {
        let new_columns = {};
        let new_order = [];

        for (var ii of i) {
            if (typeof ii != "string") {
                utils.check_entry_index(this._columns.order, ii, "column", "DataFrame");
                ii = this._columns.order[ii];
            }
            if (ii in new_columns) {
                throw new Error("duplicate columns detected in slice request");
            }

            new_columns[ii] = this._columns.entries[ii];
            new_order.push(ii);
        }

        let target = cutils.setterTarget(this, inPlace);
        target._columns = { entries: new_columns, order: new_order };
        return target;
    }

    /**
     * @param {Array} i - Array of strings or indices specifying the columns to retain in the slice.
     * This should refer to unique column names.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DataFrame instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {DataFrame} Reference to this DataFrame after slicing to the specified columns.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    $sliceColumns(i) {
        return this.sliceColumns(i, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_LENGTH() {
        return this.numberOfRows();
    }

    _bioconductor_SLICE(output, i, { allowView = false }) {
        let options = { allowView };

        let new_columns = {};
        for (const [k, v] of Object.entries(this._columns.entries)) {
            new_columns[k] = generics.SLICE(v, i, options);
        }

        let new_rowNames = (this._rowNames == null ? null : generics.SLICE(this._rowNames, i, options));

        let new_numberOfRows;
        if (i.constructor == Object) {
            new_numberOfRows = i.end - i.start;
        } else {
            new_numberOfRows = i.length;
        }

        output._rowNames = new_rowNames;
        output._columns = { entries: new_columns, order: this._columns.order };
        output._numberOfRows = new_numberOfRows;
        output._metadata = this._metadata;
        return; 
    }

    _bioconductor_COMBINE(output, objects) {
        let new_columns = utils.combineEntries(objects.map(x => x._columns), generics.COMBINE, "columnNames", "DataFrame");

        let all_n = [];
        let all_l = [];
        for (const yi of objects) {
            all_n.push(yi.rowNames());
            all_l.push(yi.numberOfRows());
        }

        let new_numberOfRows = utils.sum(all_l);
        let new_rowNames = utils.combineNames(all_n, all_l, new_numberOfRows);

        output._rowNames = new_rowNames;
        output._columns = new_columns;
        output._numberOfRows = new_numberOfRows;
        output._metadata = this._metadata;
        return;
    }

    _bioconductor_CLONE(output, { deepCopy = true }) {
        super._bioconductor_CLONE(output, { deepCopy });
        output._columns = cutils.cloneField(this._columns, deepCopy);
        output._rowNames = (this._rowNames == null ? null : cutils.cloneField(this._rowNames, deepCopy));
        output._numberOfRows = this._numberOfRows;
        return;
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

        copy._columns = { entries: copy._columns.entries, order: corder };
        copies.push(copy);
    }

    return generics.COMBINE(copies);
}
