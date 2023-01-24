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
            if (allowView && ArrayBuffer.isView(this._values)) {
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
            if (ArrayBuffer.isView(this._values)) {
                this._values.set(replacement, start);
            } else {
                for (var s = 0; s < nsecondary; s++) {
                    this._values[s + start] = replacement[s];
                }
            }
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

    _bioconductor_COMBINE_ROWS(output, objects) {
        let NC = objects[0]._numberOfColumns;
        let NR = objects[0]._numberOfRows;
        for (var i = 1; i < objects.length; i++) {
            if (objects[i]._numberOfColumns !== NC) {
                throw new Error("all objects must have the same number of columns");
            }
            NR += objects[i]._numberOfRows;
        }

        output._numberOfRows = NR;
        output._numberOfColumns = NC;
        output._columnMajor = objects[0]._columnMajor;

        let final_values = new objects[0]._values.constructor(NR * NC);
        output._values = final_values;

        if (output._columnMajor) {
            let used_rows = 0;
            for (var i = 0; i < objects.length; i++) {
                let current = objects[i];
                let currows = current._numberOfRows;

                if (current._columnMajor) {
                    for (var c = 0; c < NC; c++) {
                        let view_offset = c * currows;
                        let view = current._values.view(view_offset, view_offset + currows);
                        final_values.set(view, used_rows + c * NR);
                    }
                } else {
                    for (var r = 0; r < currows; r++) {
                        let offset = used_rows + r;
                        for (var c = 0; c < NC; c++) {
                            final_values[c * NR + offset] = current._values[r * NC + c];
                        }
                    }
                }
                 
                used_rows += currows;
            }

        } else {
            let used_rows = 0;
            for (var i = 0; i < objects.length; i++) {
                let current = objects[i];
                let currows = current._numberOfRows;

                if (!current._columnMajor) {
                    final_values.set(current._values, used_rows * NC);
                } else {
                    for (var c = 0; c < NC; c++) {
                        let offset = used_rows * NC + c;
                        for (var r = 0; r < currows; r++) {
                            final_values[offset + r * NC] = current._values[c * currows + r];
                        }
                    }
                }

                used_rows += currows;
            }
        }
    }
}
