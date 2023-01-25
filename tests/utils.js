import * as bioc from "../src/index.js";

export function spawn_random_vector(N) {
    let v = new Float64Array(N);
    v.forEach((x, i) => { v[i] = Math.random(); });
    return v;
}

export function spawn_random_matrix(NR, NC) {
    let v = spawn_random_vector(NR * NC);
    return new bioc.DenseMatrix(NR, NC, v);
}

export function spawn_random_GRanges(N) {
    let s = new Int32Array(N);
    s.forEach((x, i) => { s[i] = Math.random() * 1000; });

    let w = new Int32Array(N);
    w.forEach((x, i) => { w[i] = Math.random() * 100; });

    let n = new Array(N);
    for (var i = 0; i < N; ++i) {
        n[i] = "chr" + String(Math.random() * 3 + 1);
    }

    return new bioc.GRanges(n, new bioc.IRanges(s, w));
}
