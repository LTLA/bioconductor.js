import * as bioc from "../src/index.js";

function spawn_random_vector(n) {
    let v = new Float64Array(n);
    v.forEach((x, i) => { v[i] = Math.random(); });
    return v;
}

test("Constructing a DenseMatrix works", () => {
    let NR = 9;
    let NC = 11;
    let v = new spawn_random_vector(NR * NC);

    // Column-major.
    let mat = new bioc.DenseMatrix(NR, NC, v);
    expect(mat.numberOfRows()).toBe(NR);
    expect(mat.numberOfColumns()).toBe(NC);
    expect(mat.values()).toEqual(v);
    expect(mat.isColumnMajor()).toBe(true);

    // Row-major
    let rmat = new bioc.DenseMatrix(NR, NC, v, { columnMajor: false });
    expect(rmat.numberOfRows()).toBe(NR);
    expect(rmat.numberOfColumns()).toBe(NC);
    expect(rmat.values()).toEqual(v);
    expect(rmat.isColumnMajor()).toBe(false);

    // Error.
    expect(() => new bioc.DenseMatrix(NR - 1, NC, v)).toThrow("product");
})

test("row and column extractors work as expected", () => {
    let NR = 7;
    let NC = 6;
    let v = new spawn_random_vector(NR * NC);

    // Column-major.
    let mat = new bioc.DenseMatrix(NR, NC, v);
    for (var c = 0; c < NC; c++) {
        expect(mat.column(c)).toEqual(v.slice(c * NR, c * NR + NR));
        expect(mat.column(c, { allowView: true })).toEqual(v.slice(c * NR, c * NR + NR));
    }

    for (var r = 0; r < NR; r++) {
        let row = mat.row(r);
        let expected = new Float64Array(NC);
        for (var c = 0; c < NC; c++) {
            expected[c] = v[r + c * NR];
        }
        expect(row).toEqual(expected);
    }

    // Row-major.
    let tmat = new bioc.DenseMatrix(NC, NR, v, { columnMajor: false });
    for (var c = 0; c < NC; c++) {
        expect(mat.column(c)).toEqual(tmat.row(c));
        expect(mat.column(c)).toEqual(tmat.row(c, { allowView: true }));
    }
    for (var r = 0; r < NR; r++) {
        expect(mat.row(r)).toEqual(tmat.column(r));
    }
})

test("values setters work as expected", () => {
    let NR = 8;
    let NC = 13;
    let v = new spawn_random_vector(NR * NC);
    let mat = new bioc.DenseMatrix(NR, NC, v);

    let v2 = new spawn_random_vector(NR * NC);
    mat.$setValues(v2);
    expect(mat.values()).toEqual(v2);

    expect(() => mat.$setValues(new Float64Array(0))).toThrow("replacement");
})

test("row and column setters work as expected", () => {
    let NR = 8;
    let NC = 13;
    let v = new spawn_random_vector(NR * NC);
    let v2 = new spawn_random_vector(NR * NC);

    // Column-major.
    let mat2 = new bioc.DenseMatrix(NR, NC, v2);
    {
        let mat = new bioc.DenseMatrix(NR, NC, v.slice());
        for (var c = 0; c < NC; c++) {
            mat.$setColumn(c, mat2.column(c));
        }
        expect(mat.values()).toEqual(v2);
    }

    {
        let mat = new bioc.DenseMatrix(NR, NC, v.slice());
        for (var r = 0; r < NR; r++) {
            mat.$setRow(r, mat2.row(r));
        }
        expect(mat.values()).toEqual(v2);
    }

    // Row-major
    {
        let mat = new bioc.DenseMatrix(NC, NR, v.slice(), { columnMajor: false });
        for (var c = 0; c < NC; c++) {
            mat.$setRow(c, mat2.column(c));
        }
        expect(mat.values()).toEqual(v2);
    }

    {
        let mat = new bioc.DenseMatrix(NC, NR, v.slice(), { columnMajor: false });
        for (var r = 0; r < NR; r++) {
            mat.$setColumn(r, mat2.row(r));
        }
        expect(mat.values()).toEqual(v2);
    }

    // Errors
    {
        let mat = new bioc.DenseMatrix(NC, NR, v.slice());
        expect(() => mat.$setRow(0, new Float64Array(1))).toThrow("numberOfColumns")
        expect(() => mat.$setColumn(0, new Float64Array(1))).toThrow("numberOfRows")
    }
})

test("NUMBER_OF generics work as expected", () => {
    let NR = 8;
    let NC = 13;
    let v = new spawn_random_vector(NR * NC);

    let mat = new bioc.DenseMatrix(NR, NC, v);
    expect(bioc.NUMBER_OF_ROWS(mat)).toBe(NR);
    expect(bioc.NUMBER_OF_COLUMNS(mat)).toBe(NC);
})

test("SLICE_2D generic works as expected (none)", () => {
    let NR = 10;
    let NC = 15;
    let v = new spawn_random_vector(NR * NC);

    for (const major of [ true, false ]) {
        let mat = new bioc.DenseMatrix(NR, NC, v, { columnMajor: major });
        let sliced = bioc.SLICE_2D(mat, null, null);
        expect(sliced.numberOfRows()).toEqual(NR);
        expect(sliced.numberOfColumns()).toEqual(NC);
        expect(sliced.values()).toEqual(v);
        expect(sliced.isColumnMajor()).toBe(major);
    }
})

test("SLICE_2D generic works as expected (rows only)", () => {
    let NR = 5;
    let NC = 11;
    let v = new spawn_random_vector(NR * NC);

    for (const major of [ true, false ]) {
        let mat = new bioc.DenseMatrix(NR, NC, v, { columnMajor: major });

        {
            let sub = { start: 1, end: 4 };
            let sliced = bioc.SLICE_2D(mat, sub, null);
            expect(sliced.numberOfRows()).toEqual(3);
            expect(sliced.numberOfColumns()).toEqual(NC);

            if (major) {
                var expected = [];
                for (var i = 0; i < NC; i++) {
                    expected.push(bioc.SLICE(mat.column(i), sub));
                }
                expect(sliced.values()).toEqual(bioc.COMBINE(expected));
            } else {
                expect(sliced.values()).toEqual(bioc.COMBINE([ mat.row(1), mat.row(2), mat.row(3) ]));
            }
        }

        {
            let indices = [ 4, 3, 3, 2, 0 ];
            let sliced = bioc.SLICE_2D(mat, indices, null);
            expect(sliced.numberOfRows()).toEqual(indices.length);
            expect(sliced.numberOfColumns()).toEqual(NC);

            if (major) {
                var expected = [];
                for (var i = 0; i < NC; i++) {
                    expected.push(bioc.SLICE(mat.column(i), indices));
                }
                expect(sliced.values()).toEqual(bioc.COMBINE(expected));
            } else {
                expect(sliced.values()).toEqual(bioc.COMBINE(indices.map(i => mat.row(i))));
            }
        }
    }
})

test("SLICE_2D generic works as expected (columns only)", () => {
    let NR = 8;
    let NC = 7;
    let v = new spawn_random_vector(NR * NC);

    for (const major of [ true, false ]) {
        let mat = new bioc.DenseMatrix(NR, NC, v, { columnMajor: major });

        {
            let sub = { start: 2, end: 7 };
            let sliced = bioc.SLICE_2D(mat, null, sub);
            expect(sliced.numberOfRows()).toEqual(NR);
            expect(sliced.numberOfColumns()).toEqual(5);

            var expected = [];
            if (major) {
                for (var i = 2; i < 7; i++) {
                    expected.push(mat.column(i));
                }
            } else {
                var expected = [];
                for (var i = 0; i < NR; i++) {
                    expected.push(bioc.SLICE(mat.row(i), sub));
                }
            }
            expect(sliced.values()).toEqual(bioc.COMBINE(expected));
        }

        {
            let indices = [ 0, 1, 6, 2, 3, 0, 5, 5 ];
            let sliced = bioc.SLICE_2D(mat, null, indices);
            expect(sliced.numberOfRows()).toEqual(NR);
            expect(sliced.numberOfColumns()).toEqual(indices.length);

            if (major) {
                expect(sliced.values()).toEqual(bioc.COMBINE(indices.map(i => mat.column(i))));
            } else {
                var expected = [];
                for (var i = 0; i < NR; i++) {
                    expected.push(bioc.SLICE(mat.row(i), indices));
                }
                expect(sliced.values()).toEqual(bioc.COMBINE(expected));
            }
        }
    }
})

test("SLICE_2D generic works as expected (both rows and columns)", () => {
    let NR = 10;
    let NC = 15;
    let v = new spawn_random_vector(NR * NC);

    for (const major of [ true, false ]) {
        let mat = new bioc.DenseMatrix(NR, NC, v, { columnMajor: major });

        {
            let rsub = { start: 2, end: 7 };
            let csub = { start: 0, end: 10 };
            let sliced = bioc.SLICE_2D(mat, rsub, csub);
            expect(sliced.numberOfRows()).toEqual(5);
            expect(sliced.numberOfColumns()).toEqual(10);

            var expected = [];
            if (major) {
                for (var i = 0; i < 10; i++) {
                    expected.push(bioc.SLICE(mat.column(i), rsub));
                }
            } else {
                var expected = [];
                for (var i = 2; i < 7; i++) {
                    expected.push(bioc.SLICE(mat.row(i), csub));
                }
            }
            expect(sliced.values()).toEqual(bioc.COMBINE(expected));
        }

        {
            let rindices = [ 0, 2, 4, 6, 8 ];
            let cindices = [ 1, 3, 5, 7, 9, 11, 13 ];
            let sliced = bioc.SLICE_2D(mat, rindices, cindices);
            expect(sliced.numberOfRows()).toEqual(rindices.length);
            expect(sliced.numberOfColumns()).toEqual(cindices.length);

            var expected = [];
            if (major) {
                for (const c of cindices) {
                    expected.push(bioc.SLICE(mat.column(c), rindices));
                }
            } else {
                for (const r of rindices) {
                    expected.push(bioc.SLICE(mat.row(r), cindices));
                }
            }
            expect(sliced.values()).toEqual(bioc.COMBINE(expected));
        }
    }
})

test("COMBINE_ROWS generic works as expected", () => {
    let NC = 15;

    for (const major1 of [ true, false ]) {
        let NR1 = 10;
        let v1 = new spawn_random_vector(NR1 * NC);
        let mat1 = new bioc.DenseMatrix(NR1, NC, v1, { columnMajor: major1 });

        for (const major2 of [ true, false ]) {
            let NR2 = 20;
            let v2 = new spawn_random_vector(NR2 * NC);
            let mat2 = new bioc.DenseMatrix(NR2, NC, v2, { columnMajor: major2 });

            let combined = bioc.COMBINE_ROWS([mat1, mat2]);
            expect(combined.numberOfRows()).toBe(NR1 + NR2);
            expect(combined.numberOfColumns()).toBe(NC);
            expect(combined.isColumnMajor()).toBe(major1);

            let expected = [];
            if (major1) {
                for (var c = 0; c < NC; c++) {
                    expected.push(mat1.column(c));
                    expected.push(mat2.column(c));
                }
            } else {
                for (var r = 0; r < NR1; r++) {
                    expected.push(mat1.row(r));
                }
                for (var r = 0; r < NR2; r++) {
                    expected.push(mat2.row(r));
                }
            }
            expect(combined.values()).toEqual(bioc.COMBINE(expected));
        }
    }

    // Errors.
    {
        let NR = 8;
        let NC = 7;
        let v = new spawn_random_vector(NR * NC);
        let mat1 = new bioc.DenseMatrix(NR, NC, v);
        let mat2 = new bioc.DenseMatrix(4, 14, v);
        expect(() => bioc.COMBINE_ROWS([mat1, mat2])).toThrow("same number of columns");
    }
})

test("COMBINE_COLUMNS generic works as expected", () => {
    let NR = 21;

    for (const major1 of [ true, false ]) {
        let NC1 = 13;
        let v1 = new spawn_random_vector(NR * NC1);
        let mat1 = new bioc.DenseMatrix(NR, NC1, v1, { columnMajor: major1 });

        for (const major2 of [ true, false ]) {
            let NC2 = 11;
            let v2 = new spawn_random_vector(NR * NC2);
            let mat2 = new bioc.DenseMatrix(NR, NC2, v2, { columnMajor: major2 });

            let combined = bioc.COMBINE_COLUMNS([mat1, mat2]);
            expect(combined.numberOfRows()).toBe(NR);
            expect(combined.numberOfColumns()).toBe(NC1 + NC2);
            expect(combined.isColumnMajor()).toBe(major1);

            let expected = [];
            if (major1) {
                for (var c = 0; c < NC1; c++) {
                    expected.push(mat1.column(c));
                }
                for (var c = 0; c < NC2; c++) {
                    expected.push(mat2.column(c));
                }
            } else {
                for (var r = 0; r < NR; r++) {
                    expected.push(mat1.row(r));
                    expected.push(mat2.row(r));
                }
            }
            expect(combined.values()).toEqual(bioc.COMBINE(expected));
        }
    }

    // Errors.
    {
        let NR = 8;
        let NC = 7;
        let v = new spawn_random_vector(NR * NC);
        let mat1 = new bioc.DenseMatrix(NR, NC, v);
        let mat2 = new bioc.DenseMatrix(4, 14, v);
        expect(() => bioc.COMBINE_COLUMNS([mat1, mat2])).toThrow("same number of rows");
    }
})

test("CLONE generic works as expected", () => {
    let NR = 7;
    let NC = 11;
    let v = new spawn_random_vector(NR * NC);
    let mat = new bioc.DenseMatrix(NR, NC, v);

    let cloned = bioc.CLONE(mat);
    expect(cloned.numberOfRows()).toBe(NR);
    expect(cloned.numberOfColumns()).toBe(NC);
    expect(cloned.isColumnMajor()).toBe(true);
    expect(cloned.values()).toEqual(v);

    // Deepcopies by default.
    cloned.values()[0] = -100;
    expect(cloned.values()[0]).toBe(-100);
    expect(mat.values()[0]).not.toBe(-100);

    // Non-deep copy.
    let cloned2 = bioc.CLONE(mat, { deepCopy: false });
    cloned2.values()[0] = -100;
    expect(cloned2.values()[0]).toBe(-100);
    expect(mat.values()[0]).toBe(-100);
})
