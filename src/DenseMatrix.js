import * as cutils from "./clone-utils.js";

/**
 * Dense matrix of numbers.
 * Not really a Bioconductor-exclusive data structure, but we need this at a minimum for the {@linkplain SummarizedExperiment} to be useful.
 *
 * - {@linkcode NUMBER_OF_ROWS}
 * - {@linkcode NUMBER_OF_COLUMNS}
 * - {@linkcode SLICE_2D}
 * - {@linkcode COMBINE_ROWS}
 * - {@linkcode COMBINE_COLUMNS}
 * - {@linkcode CLONE}
 *
 * Constructors of DataFrame subclasses should be callable with no arguments, possibly creating an empty object with no properties.
 * This will be used by the `_bioconductor_CLONE`, `_bioconductor_COMBINE_ROWS`, `_bioconductor_COMBINE_COLUMNS` and `_bioconductor_SLICE_2D` methods to return an instance of the subclass.
 */
export class DenseMatrix {
    /**
     * @param {number} numberOfRows - Number of rows, duh.
     * @param {number} numberOfColumns - Number of columns.
     * @param {TypedArray} values - 1-dimensional array of the matrix contents.
     * This should have length equal to the product of `numberOfRows` and `numberOfColumns`.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.columnMajor=true] - Whether `values` represents a column-major layout.
     */
    constructor(numberOfRows, numberOfColumns, values, { columnMajor = true } = {}) {
        if (arguments.length == 0) {
            return;
        }

        this._numberOfRows = numberOfRows;
        this._numberOfColumns = numberOfColumns;
        this._values = values;
        this._columnMajor = columnMajor;
        if (numberOfRows * numberOfColumns != values.length) {
            throw new Error("length of 'values' should be equal to the product of 'dimensions'");
        }
    }

    static name = "DenseMatrix";

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {number} Number of rows.
     */
    numberOfRows() {
        return this._numberOfRows;
    }

    /**
     * @return {number} Number of columns.
     */
    numberOfColumns() {
        return this._numberOfColumns;
    }

    /**
     * @return {boolean} Whether the matrix is column-major.
     */
    isColumnMajor() {
        return this._columnMajor;
    }

    /**
     * @return {TypedArray} Matrix contents as a 1-dimensional array.
     */
    values() {
        return this._values;
    }

    #extractor(i, nprimary, nsecondary, allowView, primaryMajor) {
        if (!primaryMajor) {
            let output = new this._values.constructor(nsecondary);
            let offset = i;
            for (var s = 0; s < nsecondary; s++) {
                output[s] = this._values[offset];
                offset += nprimary;
            }
            return output;

        } else {
            let start = i * nsecondary;
            let end = start + nsecondary;
            if (allowView) {
                return this._values.subarray(start, end);
            } else {
                return this._values.slice(start, end);
            }
        }
    }

    /**
     * Retrieve the contents of a particular row.
     *
     * @param {number} i - Index of the row of interest.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.allowView=false] - Whether to allow a view to be returned, if possible.
     *
     * @return {TypedArray} Contents of the row `i`.
     * This may be a view on the array returned by {@linkcode DenseMatrix#values values}, if permitted by the layout.
     */
    row(i, { allowView = false } = {}) {
        return this.#extractor(i, this._numberOfRows, this._numberOfColumns, allowView, !this._columnMajor);
    }

    /**
     * Retrieve the contents of a particular column.
     *
     * @param {number} i - Index of the column of interest.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.allowView=false] - Whether to allow a view to be returned, if possible.
     *
     * @return {TypedArray} Contents of the column `i`.
     * This may be a view on the array returned by {@linkcode DenseMatrix#values values}, if permitted by the layout.
     */
    column(i, { allowView = false } = {}) {
        return this.#extractor(i, this._numberOfColumns, this._numberOfRows, allowView, this._columnMajor);
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {TypedArray} values - 1-dimensional array of matrix contents,
     * of the same length as the array returned by {@linkcode DenseMatrix#values values}.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DenseMatrix instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {DenseMatrix} The DenseMatrix after modifying the matrix contents.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setValues(values, { inPlace = false } = {}) {
        if (values.length !== this._values.length) {
            throw new Error("replacement 'values' should have length equal to 'values()'");
        }

        let target = cutils.setterTarget(this, inPlace);
        target._values = values;
        return target;
    }

    $setValues(values) {
        return this.setValues(values, { inPlace: true });
    }

    #inserter(i, nprimary, nsecondary, primaryMajor, replacement) {
        if (!primaryMajor) {
            let output = new this._values.constructor(nsecondary);
            let offset = i;
            for (var s = 0; s < nsecondary; s++) {
                this._values[offset] = replacement[s];
                offset += nprimary;
            }
        } else {
            let start = i * nsecondary;
            this._values.set(replacement, start);
        }
    }

    /**
     * @param {number} i - Row index to set.
     * @param {TypedArray} values - Row contents, of length equal to the number of columns in this DenseMatrix.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DenseMatrix instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {DenseMatrix} The DenseMatrix after modifying the matrix contents.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setRow(i, values, { inPlace = false } = {}) {
        if (values.length !== this._numberOfColumns) {
            throw new Error("replacement row should have length equal to 'numberOfColumns()'");
        }

        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._values = target._values.slice();
        }

        target.#inserter(i, target._numberOfRows, target._numberOfColumns, !target._columnMajor, values);
        return target;
    }

    $setRow(i, value) {
        return this.setRow(i, value, { inPlace: true });
    }

    /**
     * @param {number} i - Column index to set.
     * @param {TypedArray} values - Column contents, of length equal to the number of rows in this DenseMatrix.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this DenseMatrix instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {DenseMatrix} The DenseMatrix after modifying the matrix contents.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setColumn(i, values, { inPlace = false } = {}) {
        if (values.length !== this._numberOfRows) {
            throw new Error("replacement column should have length equal to 'numberOfRows()'");
        }

        let target = cutils.setterTarget(this, inPlace);
        if (!inPlace) {
            target._values = target._values.slice();
        }

        target.#inserter(i, target._numberOfColumns, target._numberOfRows, target._columnMajor, values);
        return target;
    }

    $setColumn(i, value) {
        return this.setColumn(i, value, { inPlace: true });
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

    _bioconductor_SLICE_2D(rows, columns, {}) {
        let output = new this.constructor;

        let full_rows = (rows === null);
        let is_row_range = (!full_rows && rows.constructor == Object);
        let new_rows = full_rows ? this._numberOfRows : (is_row_range ? rows.end - rows.start : rows.length);
        output._numberOfRows = new_rows;

        let full_columns = (columns === null);
        let is_column_range = (!full_columns && columns.constructor == Object);
        let new_columns = full_columns ? this._numberOfColumns : (is_column_range ? columns.end - columns.start : columns.length);
        output._numberOfColumns = new_columns;

        let new_values = new this._values.constructor(new_rows * new_columns);
        output._values = new_values;

        if (this._columnMajor) {
            this.#primarySlicer(columns, full_columns, is_column_range, this._numberOfColumns, rows, full_rows, is_row_range, this._numberOfRows, new_rows, new_values);
        } else {
            this.#primarySlicer(rows, full_rows, is_row_range, this._numberOfRows, columns, full_columns, is_column_range, this._numberOfColumns, new_columns, new_values);
        }
        output._columnMajor = this._columnMajor;
        return output;
    }

    #primarySlicer(primarySlice, fullPrimary, isPrimaryRange, primaryDim, secondarySlice, fullSecondary, isSecondaryRange, inSecondaryDim, outSecondaryDim, outputValues) {
        if (fullPrimary) {
            for (var p = 0; p < primaryDim; p++) {
                this.#secondarySlicer(secondarySlice, fullSecondary, isSecondaryRange, inSecondaryDim, outSecondaryDim, outputValues, p, p);
            }
        } else if (isPrimaryRange) {
            for (var p = primarySlice.start; p < primarySlice.end; p++) {
                this.#secondarySlicer(secondarySlice, fullSecondary, isSecondaryRange, inSecondaryDim, outSecondaryDim, outputValues, p, p - primarySlice.start);
            }
        } else {
            for (var pi = 0; pi < primarySlice.length; pi++) {
                this.#secondarySlicer(secondarySlice, fullSecondary, isSecondaryRange, inSecondaryDim, outSecondaryDim, outputValues, primarySlice[pi], pi);
            }
        }
    }

    #secondarySlicer(secondarySlice, fullSecondary, isSecondaryRange, inSecondaryDim, outSecondaryDim, outputValues, inPrimary, outPrimary) {
        let in_offset = inPrimary * inSecondaryDim;
        let out_offset = outPrimary * outSecondaryDim;

        if (fullSecondary) {
            let view = this._values.subarray(in_offset, in_offset + inSecondaryDim);
            outputValues.set(view, out_offset);
        } else if (isSecondaryRange) {
            for (var s = secondarySlice.start; s < secondarySlice.end; s++) {
                outputValues[out_offset + s - secondarySlice.start] = this._values[in_offset + s];
            }
        } else {
            for (var si = 0; si < secondarySlice.length; si++) {
                outputValues[out_offset + si] = this._values[in_offset + secondarySlice[si]];
            }
        }
    }

    _combiner(objects, primaryFun, secondaryFun, isPrimaryMajor, secondaryName) {
        let num_primary = primaryFun(this);
        let num_secondary = secondaryFun(this);
        for (const x of objects) {
            if (secondaryFun(x) !== num_secondary) {
                throw new Error("all objects must have the same number of " + secondaryName);
            }
            num_primary += primaryFun(x);
        }

        let primary_major = isPrimaryMajor(this);
        let values = new this._values.constructor(num_primary * num_secondary);

        if (primary_major) {
            let used_primary = 0;
            for (var i = 0; i <= objects.length; i++) {
                let current = (i == 0 ? this : objects[i - 1]);
                let cur_primary = primaryFun(current);
                let out_offset = used_primary * num_secondary;

                if (isPrimaryMajor(current)) {
                    values.set(current._values, out_offset);
                } else {
                    for (var s = 0; s < num_secondary; s++) {
                        let in_offset = s * cur_primary;
                        let out_offset2 = out_offset + s;
                        for (var p = 0; p < cur_primary; p++) {
                            values[out_offset2 + p * num_secondary] = current._values[in_offset + p];
                        }
                    }
                }

                used_primary += cur_primary;
            }
        } else {
            let used_primary = 0;
            for (var i = 0; i <= objects.length; i++) {
                let current = (i == 0 ? this : objects[i - 1]);
                let cur_primary = primaryFun(current);

                if (!isPrimaryMajor(current)) {
                    for (var s = 0; s < num_secondary; s++) {
                        let view_offset = s * cur_primary;
                        let view = current._values.subarray(view_offset, view_offset + cur_primary);
                        values.set(view, used_primary + s * num_primary);
                    }
                } else {
                    for (var p = 0; p < cur_primary; p++) {
                        let in_offset = p * num_secondary;
                        let out_offset = used_primary + p;
                        for (var s = 0; s < num_secondary; s++) {
                            values[out_offset + s * num_primary] = current._values[in_offset + s];
                        }
                    }
                }

                used_primary += cur_primary;
            }
        }

        return { num_primary, num_secondary, values, primary_major };
    }

    _bioconductor_COMBINE_ROWS(objects) {
        let combined = this._combiner(objects,
            x => x._numberOfRows,
            x => x._numberOfColumns,
            x => !(x._columnMajor),
            "columns"
        );

        let output = new this.constructor;
        output._numberOfRows = combined.num_primary;
        output._numberOfColumns = combined.num_secondary;
        output._values = combined.values;
        output._columnMajor = !(combined.primary_major);
        return output;
    }

    _bioconductor_COMBINE_COLUMNS(objects) {
        let combined = this._combiner(objects,
            x => x._numberOfColumns,
            x => x._numberOfRows,
            x => x._columnMajor,
            "rows"
        );

        let output = new this.constructor;
        output._numberOfColumns = combined.num_primary;
        output._numberOfRows = combined.num_secondary;
        output._values = combined.values;
        output._columnMajor = combined.primary_major;
        return output;
    }

    _bioconductor_CLONE({ deepCopy = true } = {}) {
        let output = new this.constructor;
        output._values = (deepCopy ? this._values.slice() : this._values);
        output._numberOfRows = this._numberOfRows;
        output._numberOfColumns = this._numberOfColumns;
        output._columnMajor = this._columnMajor;
        return output;
    }
}
