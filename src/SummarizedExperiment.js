import * as generics from "./AllGenerics.js";
import * as ann from "./Annotated.js";
import * as df from "./DataFrame.js";
import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as il from "./InternalList.js";

/**
 * A SummarizedExperiment contains zero or more assays, consisting of multi-dimensional arrays (usually matrices) of experimental data,
 * as well as {@linkplain DataFrame}s containing further annotations on the rows or columns of those arrays.
 * The SummarizedExperiment class defines methods for the following generics:
 * 
 * - {@linkcode NUMBER_OF_ROWS}
 * - {@linkcode NUMBER_OF_COLUMNS}
 * - {@linkcode SLICE_2D}
 * - {@linkcode COMBINE_ROWS}
 * - {@linkcode COMBINE_COLUMNS}
 * - {@linkcode CLONE}
 *
 * Assays are expected to provide methods for the following generics:
 *
 * - {@linkcode NUMBER_OF_ROWS}
 * - {@linkcode NUMBER_OF_COLUMNS}
 * - {@linkcode SLICE_2D}
 * - {@linkcode COMBINE_ROWS}
 * - {@linkcode COMBINE_COLUMNS}
 * - {@linkcode CLONE}
 *
 * Constructors of SummarizedExperiment subclasses should be callable with no arguments, possibly creating an empty object with no properties.
 * This will be used by the `_bioconductor_CLONE`, `_bioconductor_COMBINE_ROWS`, `_bioconductor_COMBINE_COLUMNS` and `_bioconductor_SLICE_2D` methods to return an instance of the subclass.
 *
 * @extends Annotated
 */
export class SummarizedExperiment extends ann.Annotated {
    /**
     * @param {Object|Map} assays - Object or Map where keys are the assay names and values are multi-dimensional arrays of experimental data.
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
     * @param {Object|Array|Map|List} [options.metadata={}] - Arbitrary metadata, see the {@link Annotated} constructor. 
     */
    constructor(assays, { assayOrder = null, rowData = null, columnData = null, rowNames = null, columnNames = null, metadata = {} } = {}) {
        if (arguments.length == 0) {
            super();
            return;
        }

        super(metadata);

        // Check the assays.
        try {
            this._assays = new il.InternalList(assays, assayOrder);
        } catch (e) {
            throw new Error("failed to initialize assay list for this SummarizedExperiment; " + e.message, { cause: e });
        }

        let nrows = null;
        let ncols = null;
        for (const k of this._assays.names()) {
            let current = this._assays.entry(k);
            let nr = generics.NUMBER_OF_ROWS(current);
            let nc = generics.NUMBER_OF_COLUMNS(current);
            if (nrows == null) {
                nrows = nr;
                ncols = nc;
            } else if (nrows !== nr || ncols !== nc) {
                throw new Error("expected all assays in 'assays' to have the same number of rows and columns");
            }
        }

        // Check the rowData.
        if (rowData === null) {
            if (nrows == null){
                throw new Error("'rowData' must be specified if 'assays' is empty");
            }
            rowData = new df.DataFrame({}, { numberOfRows: nrows });
        } else {
            if (nrows !== null && nrows !== generics.LENGTH(rowData)) {
                throw new Error("'rowData' should be equal to the number of rows in each 'assays'");
            }
        }
        this._rowData = rowData;

        // Check the columnData.
        if (columnData === null) {
            if (ncols == null){
                throw new Error("'columnData' must be specified if 'assays' is empty");
            }
            columnData = new df.DataFrame({}, { numberOfRows: ncols });
        } else {
            if (ncols !== null && ncols !== generics.LENGTH(columnData)) {
                throw new Error("'columnData' should be equal to the number of columns in each 'assays'");
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

    static className = "SummarizedExperiment";

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {Array} Array of assay names.
     */
    assayNames() {
        return this._assays.names();
    }

    /**
     * @return {number} Number of assays.
     */
    numberOfAssays() {
        return this._assays.numberOfEntries();
    }

    /**
     * @param {string|number} i - Assay to retrieve, either by name or index.
     * @return {*} The contents of assay `i` as an multi-dimensional array-like object.
     */
    assay(i) {
        let output;
        try {
            output = this._assays.entry(i);
        } catch (e) {
            throw new Error("failed to retrieve the specified assay from this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return output;
    }

    /**
     * @return {Map} Map where keys are the assay names and values are the assays.
     */
    assays() {
        return this._assays.entries();
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
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SummarizedExperiment} The SummarizedExperiment after removing the specified assay.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    removeAssay(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._assays = target._assays.delete(i, { inPlace });
        } catch (e) {
            throw new Error("failed to remove assay " + (typeof i == "string" ? "'" + i + "'" : String(i)) + " from this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $removeAssay(i) {
        return this.removeAssay(i, { inPlace: true });
    }

    /**
     * @param {string|number} i - Identity of the assay to add, either by name or index.
     * - If `i` is a number, the assay at the specified index is replaced.
     *   `i` should be non-negative and less than the number of assays.
     * - If `i` is a string, any assay with the same name is replaced.
     *   If no such assay exists, a new assay is appended to the list of assays.
     * @param {*} value - Multi-dimensional array-like object to set/add as the assay.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SummarizedExperiment} A SummarizedExperiment with modified assays.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setAssay(i, value, { inPlace = false } = {}) {
        if (generics.NUMBER_OF_ROWS(value) !== this.numberOfRows() || generics.NUMBER_OF_COLUMNS(value) !== this.numberOfColumns()) {
            throw new Error("expected 'value' to have the same dimensions as this 'SummarizedExperiment'");
        }
        let target = cutils.setterTarget(this, inPlace);
        target._assays = target._assays.set(i, value, { inPlace });
        return target;
    }

    /**
     * @param {Object|Map} value - Object containing zero, one or more multi-dimensional array-like objects in the values.
     * All arrays should have the same number of rows and columns.
     * Keys are assay names, each of which should be present in `order`.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     * @param {Array|boolean} [options.newOrder=false] - Whether to replace the order of assays with the order of keys in `value`.
     * If `false`, the existing order in {@linkcode SummarizedExperiment#assayNames assayNames} is used.
     * If an array is provided, this is used as the order.
     * If `null`, this has the same effect as `true`.
     *
     * @return {SummarizedExperiment} A SummarizedExperiment with modified assays.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setAssays(value, { inPlace = false, newOrder = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);

        if (newOrder === false) {
            newOrder = target._assays.names();
        } else if (newOrder == true) {
            newOrder = null;
        }
        try {
            target._assays = new il.InternalList(value, newOrder);
        } catch (e) {
            throw new Error("failed to replace assay list for this SummarizedExperiment; " + e.message, { cause: e });
        }

        let se_nr = target.numberOfRows();
        let se_nc = target.numberOfColumns();
        for (const k of target._assays.names()) {
            let current = target._assays.entry(k);
            let nr = generics.NUMBER_OF_ROWS(current);
            let nc = generics.NUMBER_OF_COLUMNS(current);
            if (nr !== se_nr || nc !== se_nc) {
                throw new Error("mismatch in the dimensions of assay '" + k + "' compared to the SummarizedExperiment");
            }
        }

        return target;
    }

    $setAssays(value, { newOrder = false } = {}) {
        return this.setAssays(value, { inPlace: true, newOrder });
    }

    $setAssay(i, value) {
        return this.setAssay(i, value, { inPlace: true });
    }

    /**
     * @param {Array} names - Array of strings containing the assay names.
     * This should be of the same length as the number of assays and contain unique values.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SummarizedExperiment} The SummarizedExperiment with modified assay names.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setAssayNames(names, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._assays = target._assays.setNames(names, { inPlace });
        } catch (e) {
            throw new Error("failed to set the assay names for this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $setAssayNames(names) {
        return this.setAssayNames(names, { inPlace: true });
    }

    /**
     * @param {Array} i - Array of strings or indices specifying the assays to retain in the slice.
     * This should refer to unique assay names.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SummarizedExperiment} The SummarizedExperiment with sliced assays.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    sliceAssays(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._assays = this._assays.slice(i, { inPlace });
        } catch (e) {
            throw new Error("failed to slice the assays for this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $sliceAssays(i) {
        return this.sliceAssays(i, { inPlace: true });
    }

    /**
     * @param {DataFrame} value - Data frame containing the row annotations.
     * This should have one row for each row of this SummarizedExperiment.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SummarizedExperiment} The SummarizedExperiment with modified row data.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setRowData(value, { inPlace = false } = {}) {
        if (!(value instanceof df.DataFrame)) {
            throw new Error("'value' should be a DataFrame");
        }

        if (value.numberOfRows() !== this.numberOfRows()) {
            throw new Error("expected 'value' to have the same number of rows as this 'SummarizedExperiment'");
        }

        let target = cutils.setterTarget(this, inPlace);
        target._rowData = value;
        return target;
    }

    $setRowData(value) {
        return this.setRowData(value, { inPlace: true });
    }

    /**
     * @param {DataFrame} value - Data frame containing the column annotations.
     * This should have one row for each columns of this SummarizedExperiment.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SummarizedExperiment} The SummarizedExperiment with modified column data.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setColumnData(value, { inPlace = false } = {}) {
        if (!(value instanceof df.DataFrame)) {
            throw new Error("'value' should be a DataFrame");
        }

        if (value.numberOfRows() !== this.numberOfColumns()) {
            throw new Error("expected 'value' to have the same number of rows as the number of columns of this 'SummarizedExperiment'");
        }

        let target = cutils.setterTarget(this, inPlace);
        target._columnData = value;
        return target;
    }

    $setColumnData(value) {
        return this.setColumnData(value, { inPlace: true });
    }

    /**
     * @param {Array} names - Array of strings of length equal to the number of rows in this SummarizedExperiment, containing row names.
     * Alternatively `null`, to remove all row names.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SummarizedExperiment} The SummarizedExperiment with modified row names.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setRowNames(names, { inPlace = false } = {}) {
        if (names !== null) {
            utils.checkNamesArray(names, "replacement 'names'", this.numberOfRows(), "'numberOfRows()'");
        }

        let target = cutils.setterTarget(this, inPlace);
        target._rowNames = names;
        return target;
    }

    $setRowNames(names) {
        return this.setRowNames(names, { inPlace: true });
    }

    /**
     * @param {Array} names - Array of strings of length equal to the number of columns in this SummarizedExperiment, containing column names.
     * Alternatively `null`, to remove all column names.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SummarizedExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SummarizedExperiment} The SummarizedExperiment with modified column names.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setColumnNames(names, { inPlace = false } = {}) {
        if (names !== null) {
            utils.checkNamesArray(names, "replacement 'names'", this.numberOfColumns(), "'numberOfColumns()'");
        }

        let target = cutils.setterTarget(this, inPlace);
        target._columnNames = names;
        return target;
    }

    $setColumnNames(names) {
        return this.setColumnNames(names, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_NUMBER_OF_ROWS() {
        return this.numberOfRows();
    }

    _bioconductor_NUMBER_OF_COLUMNS() {
        return this.numberOfColumns();
    }

    _bioconductor_SLICE_2D(rows, columns, { allowView = false }) {
        let output = new this.constructor;
        output._assays = this._assays.apply(v => generics.SLICE_2D(v, rows, columns, { allowView }));

        if (rows !== null) {
            output._rowData = generics.SLICE(this._rowData, rows, { allowView });
            output._rowNames = (this._rowNames == null ? null : generics.SLICE(this._rowNames, rows, { allowView }));
        } else {
            output._rowData = this._rowData;
            output._rowNames = this._rowNames;
        }

        if (columns !== null) {
            output._columnData = generics.SLICE(this._columnData, columns, { allowView });
            output._columnNames = (this._columnNames == null ? null : generics.SLICE(this._columnNames, columns, { allowView }));
        } else {
            output._columnData = this._columnData;
            output._columnNames = this._columnNames;
        }

        output._metadata = this._metadata;
        return output;
    }

    _bioconductor_COMBINE_ROWS(objects) {
        let all_assays = [this._assays];
        let all_dfs = [this._rowData];
        let all_n = [this._rowNames];
        let all_l = [this.numberOfRows()];
        for (const x of objects) {
            all_assays.push(x._assays);
            all_dfs.push(x._rowData);
            all_n.push(x._rowNames);
            all_l.push(x.numberOfRows());
        }

        let output = new this.constructor;
        output._assays = il.InternalList.parallelCombine(all_assays, generics.COMBINE_ROWS);
        output._rowData = generics.COMBINE(all_dfs);
        output._rowNames = utils.combineNames(all_n, all_l);

        output._columnData = this._columnData;
        output._columnNames = this._columnNames;
        output._metadata = this._metadata;
        return output;
    }

    _bioconductor_COMBINE_COLUMNS(objects) {
        let all_assays = [this._assays];
        let all_dfs = [this._columnData];
        let all_n = [this._columnNames];
        let all_l = [this.numberOfColumns()];
        for (const x of objects) {
            all_assays.push(x._assays);
            all_dfs.push(x._columnData);
            all_n.push(x._columnNames);
            all_l.push(x.numberOfColumns());
        }

        let output = new this.constructor;
        output._assays = il.InternalList.parallelCombine(all_assays, generics.COMBINE_COLUMNS);
        output._columnData = generics.COMBINE(all_dfs);
        output._columnNames = utils.combineNames(all_n, all_l);

        output._rowData = this._rowData;
        output._rowNames = this._rowNames;
        output._metadata = this._metadata;
        return output;
    }

    _bioconductor_CLONE({ deepCopy = true }) {
        let output = super._bioconductor_CLONE({ deepCopy });

        output._assays = cutils.cloneField(this._assays, deepCopy);
        output._rowData = cutils.cloneField(this._rowData, deepCopy);
        output._rowNames = cutils.cloneField(this._rowNames, deepCopy);

        output._columnData = cutils.cloneField(this._columnData, deepCopy);
        output._columnNames = cutils.cloneField(this._columnNames, deepCopy);
        return output;
    }
}
