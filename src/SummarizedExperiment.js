import * as generics from "./AllGenerics.js";
import * as ann from "./Annotated.js";
import * as df from "./DataFrame.js";
import * as utils from "./utils.js";

/**
 * A SummarizedExperiment contains zero or more assays, consisting of multi-dimensional arrays (usually matrices) of experimental data,
 * as well as {@linkplain DataFrame}s containing further annotations on the rows or columns of those arrays.
 * The SummarizedExperiment class defines methods for the following generics:
 * 
 * - {@linkcode NUMBER_OF_ROWS}
 * - {@linkcode NUMBER_OF_COLUMNS}
 * - {@linkcode SLICE_ROWS}
 * - {@linkcode SLICE_COLUMNS}
 * - {@linkcode COMBINE_ROWS}
 * - {@linkcode COMBINE_COLUMNS}
 *
 * @extends Annotated
 */
export class SummarizedExperiment extends ann.Annotated {
    /**
     * @param {Object} assays - Object where keys are the assay names and values are multi-dimensional arrays of experimental data.
     * All arrays should have the same number of rows and columns.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?Array} [options.assayOrder=null] - Array of strings specifying the ordering of the assays.
     * If non-`null`, this should have the same values as the keys of `assays`.
     * If `null`, an arbitrary ordering is obtained from `assays`.
     * @param {?DataFrame} [options.rowData=null] - Data frame of row annotations.
     * If non-`null`, this should have a number of rows equal to the number of rows in each entry of `assays`.
     * If `null`, an empty {@linkplain DataFrame} is automatically created.
     * @param {?DataFrame} [options.columnData=null] - Data frame of column annotations.
     * If non-`null`, this should have a number of columns equal to the number of columns in each entry of `assays`.
     * If `null`, an empty {@linkplain DataFrame} is automatically created.
     * @param {?Array} [options.rowNames=null] - Array of strings of length equal to the number of rows in the `assays`, containing row names.
     * Alternatively `null`, if no row names are present.
     * @param {?Array} [options.columnNames=null] - Array of strings of length equal to the number of columns in the `assays`, containing column names.
     * Alternatively `null`, if no column names are present.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(assays, { assayOrder = null, rowData = null, columnData = null, rowNames = null, columnNames = null, metadata = {} } = {}) {
        super({ metadata });

        // Check the assays.
        let vals = Object.values(assays);
        let nrows = null, ncols = null;
        if (vals.length) {
            for (var i = 0; i < vals.length; i++) {
                let nr = generics.NUMBER_OF_ROWS(vals[i]);
                let nc = generics.NUMBER_OF_COLUMNS(vals[i]);
                if (i == 0) {
                    nrows = nr;
                    ncols = nc;
                } else if (nrows !== nr || ncols !== nc) {
                    throw new Error("expected all assays in 'assays' to have the same first two dimensions");
                }
            }
        }

        if (assayOrder == null) {
            assayOrder = assays.keys();
        } else {
            utils.checkNamesArray(assayOrder, "'assayOrder'", vals.length, "the number of entries in 'assays'");
            let uniq = Array.from(new Set(assayOrder));
            uniq.sort();
            let expected = Object.keys(assays);
            expected.sort();
            if (!utils.areArraysEqual(uniq, expected)) {
                throw new Error("values of 'assayOrder' should be the same as the keys of 'assays'");
            }
        }

        // Check the rowData.
        if (rowData === null) {
            if (nrows == null){
                throw new Error("'rowData' must be specified if 'assays' is empty");
            }
            rowData = new df.DataFrame({}, { numberOfRows: nrows });
        } else {
            if (nrows !== generics.LENGTH(rowData)) {
                throw new Error("'rowData' should be equal to the number of rows in each 'assays'");
            }
        }
        this._rowData = rowData;

        // Check the columnData.
        if (columnData === null) {
            if (ncols == null){
                thcol new Error("'columnData' must be specified if 'assays' is empty");
            }
            columnData = new df.DataFrame({}, { numberOfRows: ncols });
        } else {
            if (ncols !== generics.LENGTH(columnData)) {
                thcol new Error("'columnData' should be equal to the number of columns in each 'assays'");
            }
        }
        this._columnData = columnData;

        // Checking the names.
        if (rowNames != null) {
            utils.checkNamesArray(rowNames, "'rowNames'", this._rowData.numberOfRows(), "the number of rows in each 'assays'");
        }
        this._rowNames = rowNames;

        if (columnNames != null) {
            utils.checkNamesArray(columnNames, "'columnNames'", this._columnData.numberOfRows(), "the number of columns in each 'assays'");
        }
        this._columnNames = columnNames;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _check_index(i) {
        if (i < 0 || i >= this._assayOrder.length) {
            throw new Error("assay index '" + String(i) + "' out of range for this SummarizedExperiment");
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {Array} Array of assay names.
     */
    assayNames() {
        return this._assayOrder;
    }

    /**
     * @return {number} Number of assays.
     */
    numberOfAssays() {
        return this._assayOrder.length;
    }

    /**
     * @param {string|number} i - Assay to retrieve, either by name or index.
     * @return {*} The contents of assay `i` as an multi-dimensional array-like object.
     */
    assay(i) {
        if (typeof i == "string") {
            if (!(i in this._assays)) {
                throw new Error("no assay '" + i + "' present in this SummarizedExperiment");
            }
            return this._assays[i];
        } 

        this._check_index(i);
        return this._assays[this._assayOrder[i]];
    }

    /**
     * @return {DataFrame} Data frame of row data, with one row per row in this SummarizedExperiment.
     */
    rowData() {
        return this._rowData;
    }

    /**
     * @return {number} Number of rows in this SummarizedExperiment.
     */
    numberOfRows() {
        return this._rowData.numberOfRows();
    }

    /**
     * @return {?Array} Array of strings containing row names, or `null` if no row names are available.
     */
    rowNames() {
        return this._rowNames;
    }

    /**
     * @return {DataFrame} Data frame of column data, with one row per column in this SummarizedExperiment.
     */
    columnData() {
        return this._columnData;
    }

    /**
     * @return {number} Number of columns in this SummarizedExperiment.
     */
    numberOfColumns() {
        return this._columnData.numberOfRows();
    }

    /**
     * @return {?Array} Array of strings containing column names, or `null` if no column names are available.
     */
    columnNames() {
        return this._columnNames;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {string|number} i - Identity of the assay to add, either by name or index.
     * Numeric `i` should be non-negative and less than the number of columns.
     * @param {*} value - Multi-dimensional array-like object to set/add as the assay.
     * @return {SummarizedExperiment} Reference to this SummarizedExperiment with modified assays.
     * - If `i` is a number, the assay at the specified index is replaced.
     * - If `i` is a string, any assay with the same name is replaced.
     *   If no such assay exists, a new assay is appended to the list of assays.
     */
    $setAssay(i, value) {
        if (generics.NUMBER_OF_ROWS(value) !== this.numberOfRows() || generics.NUMBER_OF_COLUMNS(value) !== this.numberOfColumns()) {
            throw new Error("expected 'value' to have the same dimensions as this 'SummarizedExperiment'");
        }

        if (typeof i == "string") {
            if (!(i in this._assays)) {
                this._assayOrder.push(i);
            }
            this._assays[i] = value;
        } else {
            this._check_index(i);
            this._assays[this._assayOrder[i]] = value;
        }
        return this;
    }

    /**
     * @param {DataFrame} value - Data frame containing the row annotations.
     * This should have one row for each row of this SummarizedExperiment.
     * @return {SummarizedExperiment} Reference to this SummarizedExperiment with modified row data.
     */
    $setRowData(value) {
        if (value instanceof df.DataFrame) {
            throw new Error("'value' should be a DataFrame");
        }
        if (value.numberOfRows() !== this.numberOfRows()) {
            throw new Error("expected 'value' to have the same number of rows as this 'SummarizedExperiment'");
        }
        this._rowData = value;
        return this;
    }

    /**
     * @param {DataFrame} value - Data frame containing the column annotations.
     * This should have one row for each columns of this SummarizedExperiment.
     * @return {SummarizedExperiment} Reference to this SummarizedExperiment with modified column data.
     */
    $setColumnData(value) {
        if (value instanceof df.DataFrame) {
            throw new Error("'value' should be a DataFrame");
        }
        if (value.numberOfColumns() !== this.numberOfColumns()) {
            throw new Error("expected 'value' to have the same number of rows as the number of columns of this 'SummarizedExperiment'");
        }
        this._columnData = value;
        return this;
    }

    /**
     * @param {Array} names - Array of strings of length equal to the number of rows in this SummarizedExperiment, containing row names.
     * Alternatively `null`, to remove all row names.
     * @return {SummarizedExperiment} Reference to this SummarizedExperiment with modified row names.
     */
    $setRowNames(names) {
        if (names !== null) {
            utils.checkNamesArray(names, "replacement 'names'", this.numberOfRows(), "'numberOfRows()'");
        }
        this._rowNames = names;
        return this;
    }

    /**
     * @param {Array} names - Array of strings of length equal to the number of columns in this SummarizedExperiment, containing column names.
     * Alternatively `null`, to remove all column names.
     * @return {SummarizedExperiment} Reference to this SummarizedExperiment with modified column names.
     */
    $setColumnNames(names) {
        if (names !== null) {
            utils.checkNamesArray(names, "replacement 'names'", this.numberOfColumns(), "'numberOfColumns()'");
        }
        this._columnNames = names;
        return this;
    }
}
