import * as generics from "./AllGenerics.js";
import * as ann from "./Annotated.js";
import * as df from "./DataFrame.js";
import * as utils from "./utils.js";

export class SummarizedExperiment extends ann.Annotated {
    constructor(assays, { assayOrder = null, rowData = null, columnData = null, metadata = null, rowNames = null, columnNames = null } = {}) {
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

    $setAssay(i, value) {
        if (generics.LENGTH(value) != this._numberOfRows) {
            throw new Error("expected 'value' to have the same length as the number of rows in 'x'");
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
}
