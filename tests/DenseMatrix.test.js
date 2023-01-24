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
})

test("NUMBER_OF generics work as expected", () => {
    let NR = 8;
    let NC = 13;
    let v = new spawn_random_vector(NR * NC);

    let mat = new bioc.DenseMatrix(NR, NC, v);
    expect(bioc.NUMBER_OF_ROWS(mat)).toBe(NR);
    expect(bioc.NUMBER_OF_COLUMNS(mat)).toBe(NC);
})

test("SLICE_2D generic works as expected (rows only)", () => {
    let NR = 5;
    let NC = 11;
    let v = new spawn_random_vector(NR * NC);

    for (const major of [ true, false ]) {
        let mat = new bioc.DenseMatrix(NR, NC, v, { columnMajor: major });

        {
            let sliced = bioc.SLICE_2D(mat, null, null);
            expect(sliced.numberOfRows()).toEqual(NR);
            expect(sliced.numberOfColumns()).toEqual(NC);
            expect(sliced.values()).toEqual(v);
            expect(sliced.isColumnMajor()).toBe(major);
        }

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

