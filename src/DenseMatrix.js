export class DenseMatrix {
    constructor(numberOfRows, numberOfColumns, values, { columnMajor = true } = {}) {
        this._numberOfRows = numberOfRows;
        this._numberOfColumns = numberOfColumns;
        this._values = values;
        this._columnMajor = columnMajor;
        if (numberOfRows * numberOfColumns != values.length) {
            throw new Error("length of 'values' should be equal to the product of 'dimensions'");
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    numberOfRows() {
        return this._numberOfRows;
    }

    numberOfColumns() {
        return this._numberOfColumns;
    }

    values() {
        return this._values;
    }

    #extractor(i, nprimary, nsecondary, allowView, primaryFastest) {
        if (primaryFastest) {
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

    row(i, { allowView = false } = {}) {
        return this.#extractor(i, this._numberOfRows, this._numberOfColumns, allowView, !this._columnMajor);
    }

    column(i, { allowView = false } = {}) {
        return this.#extractor(i, this._numberOfColumns, this._numberOfRows, allowView, this._columnMajor);
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    $setValues(values) {
        if (values.length !== this._values.length) {
            throw new Error("replacement 'values' should have length equal to 'values()'");
        }
        this._values = values;
        return this;
    }

    #inserter(i, nprimary, nsecondary, primaryFastest, replacement) {
        if (primaryFastest) {
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

    $setRow(i, values) {
        if (values.length !== this._numberOfColumns) {
            throw new Error("replacement row should have length equal to 'numberOfColumns()'");
        }
        this.#inserter(i, this._numberOfRows, this._numberOfColumns, !this._columnMajor, values);
        return this;
    }

    $setColumn(i, values) {
        if (values.length !== this._numberOfRows) {
            throw new Error("replacement column should have length equal to 'numberOfRows()'");
        }
        this.#inserter(i, this._numberOfColumns, this._numberOfRows, this._columnMajor, values);
        return this;
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

    _bioconductor_SLICE_2D(output, rows, columns, {}) {
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
    }

    #primarySlicer(primarySlice, fullPrimary, isPrimaryRange, primaryDim, secondarySlice, fullSecondary, isSecondaryRange, inSecondaryDim, outSecondaryDim, outputValues) {
        if (fullPrimary) {
            for (var p = 0; p < primaryDim; p++) {
                this.#secondarySlicer(secondarySlice, fullSecondary, isSecondaryRange, inSecondaryDim, outSecondaryDim, outputValues, p, p);
            }
        } else if (isPrimaryRange) {
            for (var p = primarySlice.start; p < primarySlice.end; p++) {
                this.#secondarySlicer(secondarySlice, fullSecondary, isSecondaryRange, inSecondaryDim, outSecondaryDim, outputValues, p - primarySlice.start, p);
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
                outputValues[outputOffset + s - secondarySlice.start] = this._values[in_offset + s];
            }
        } else {
            for (var si = 0; si < secondarySlice.length; s++) {
                outputValues[outputOffset + si] = this._values[in_offset + secondarySlice[si]];
            }
        }
    }

    #combiner(objects, primaryFun, secondaryFun, isPrimaryMajor, secondaryName) {
        let num_primary = primaryFun(objects[0]);
        let num_secondary = secondaryFun(objects[0]);
        for (var i = 1; i < objects.length; i++) {
            if (secondaryFun(objects[i]) !== num_secondary) {
                throw new Error("all objects must have the same number of " + secondaryName);
            }
            num_primary += primaryFun(objects[i]);
        }

        let primary_major = isPrimaryMajor(objects[0]);
        let values = new objects[0]._values.constructor(num_primary * num_secondary);

        if (primary_major) {
            let used_primary = 0;
            for (var i = 0; i < objects.length; i++) {
                let current = objects[i];
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
            for (var i = 0; i < objects.length; i++) {
                let current = objects[i];
                let cur_primary = primaryFun(current);

                if (!isPrimaryMajor(current)) {
                    for (var s = 0; s < num_secondary; s++) {
                        let view_offset = s * cur_primary;
                        let view = current._values.view(view_offset, view_offset + cur_primary);
                        final_values.set(view, used_primary + s * num_primary);
                    }
                } else {
                    for (var p = 0; p < cur_primary; p++) {
                        let in_offset = p * num_secondary;
                        let out_offset = used_primary + p;
                        for (var s = 0; s < num_secondary; s++) {
                            final_values[out_offset + s * num_primary] = current._values[in_offset + s];
                        }
                    }
                }

                used_primary += cur_primary;
            }
        }

        return { num_primary, num_secondary, values, primary_major };
    }

    _bioconductor_COMBINE_ROWS(output, objects) {
        let combined = this.#combiner(objects,
            x => x._numberOfRows,
            x => x._numberOfColumns,
            x => !(x._columnMajor),
            "columns"
        );

        output._numberOfRows = combined.num_primary;
        output._numberOfColumns = combined.num_secondary;
        output._values = combined.values;
        output._columnMajor = !(combined.primary_major);
        return;
    }

    _bioconductor_COMBINE_COLUMNS(output, objects) {
        let combined = this.#combiner(objects,
            x => x._numberOfColumns,
            x => x._numberOfRows,
            x => x._columnMajor,
            "rows"
        );

        output._numberOfColumns = combined.num_primary;
        output._numberOfRows = combined.num_secondary;
        output._values = combined.values;
        output._columnMajor = combined.primary_major;
        return;
    }

    _bioconductor_CLONE(output, { deepCopy = true } = {}) {
        output._values = (deepCopy ? this._values.slice() : this._values);
        output._numberOfRows = this._numberOfRows;
        output._numberOfColumns = this._numberOfColumns;
        output._columnMajor = this._columnMajor;
        return;
    }
}
