# Bioconductor objects in Javascript

🚧🚧🚧🚧 **Under construction!** 🚧🚧🚧🚧

This package aims to provide Javascript implementations of [Bioconductor](https://github.com/Bioconductor) data structures for use in web applications.
Much like the original R code, we focus on the use of common generics to provide composability, allowing users to construct complex objects that "just work".

## Quick start

Here, we perform some generic operations on a `DataFrame` object, equivalent to Bioconductor's `S4Vectors::DFrame` class.

```js
// Import using ES6 notation
import * as bioc from "bioconductor";

// Construct a DataFrame
let results = new bioc.DataFrame(
    { 
        logFC: new Float64Array([-1, -2, 1.3, 2.1]),
        pvalue: new Float64Array([0.01, 0.02, 0.001, 1e-8])
    },
    {
        rowNames: [ "p53", "SNAP25", "MALAT1", "INS" ]
    }
);

// Run generics
bioc.LENGTH(results);
bioc.SLICE(results, [ 2, 3, 1 ]); 
bioc.CLONE(results);

let more_results = new bioc.DataFrame(
    { 
        logFC: new Float64Array([0, 0.1, -0.1]),
        pvalue: new Float64Array([1e-5, 1e-4, 0.5])
    },
    {
        rowNames: [ "GFP", "mCherry", "tdTomato" ]
    }
);

bioc.COMBINE([results, more_results]);
```

See the [reference documentation](https://ltla.github.io/bioconductor.js) for more details.

# Representing (genomic) ranges

We can construct equivalents of Bioconductor's `IRanges` and `GRanges` objects, representing integer and genomic ranges respectively.

```js
let ir = new bioc.IRanges(/* start = */ [1,2,3], /* width = */ [ 10, 20, 30 ]);
let gr = new bioc.GRanges([ "chrA", "chrB", "chrC" ], ir, { strand: [ 1, 0, -1 ] });

// Generics still work on these range objects:
bioc.LENGTH(gr);
bioc.SLICE(gr, [ 2, 1, 0 ]);
bioc.CLONE(gr);
```

We can find overlaps between two sets of ranges, akin to Bioconductor's `findOverlaps()` function:

```js
let index = gr.buildOverlapIndex();
let gr2 = new bioc.GRanges([ "chrA", "chrC", "chrA" ], new bioc.IRanges([5, 3, 2], [9, 9, 9]));
let overlaps = index.overlap(gr2);
```

We can store per-range metadata in the `elementMetadata` field of each object, just like Bioconductor's `mcols()`.

```js
let meta = gr.elementMetadata();
meta.$setColumn("symbol", [ "Nanog", "Snap25", "Malat1" ]);
gr.$setElementMetadata(meta);
gr.elementMetadata().columnNames();
```

# Handling experimental data

The `SummarizedExperiment` object is a data structure for storing experimental data in a matrix-like object, 
along with further annotations on the rows (usually features) and samples (usually columns).
To illustrate, let's mock up a small count matrix, ostensibly from an RNA-seq experiment:

```js
// Making a column-major dense matrix of random data.
let ngenes = 100;
let nsamples = 20;
let expression = new Int32Array(ngenes * nsamples);
expression.forEach((x, i) => expression[i] = Math.random() * 10);
let mat = new bioc.DenseMatrix(ngenes, nsamples, expression);

// Mocking up row names, column annotations.
let rownames = [];
for (var g = 0; g < ngenes; g++) {
    rownames.push("Gene_" + String(g));
}

let treatment = new Array(nsamples);
treatment.fill("control", 0, 10);
treatment.fill("treated", 10, nsamples);
let sample_meta = new bioc.DataFrame({ group: treatment });
```

We can now store all of this information in a `SummarizedExperiment`:

```js
let se = new bioc.SummarizedExperiment({ counts: mat }, 
    { rowNames: rownames, columnData: sample_meta });
```

This can be manipulated by generics for two-dimensional objects:

```js
bioc.NUMBER_OF_ROWS(se);
bioc.SLICE_2D(se, { start: 0, end: 50 }, [0, 2, 4, 8, 10, 12, 14, 16, 18]);
bioc.COMBINE_COLUMNS([se, se]);
```

# Using generics

Our generics allow us to operate on different objects in a consistent manner.
For example, a `DataFrame` allows us to store any object as a column as long as it defines methods for the `LENGTH`, `SLICE`, `CLONE` and `COMBINE` generics.
This allows us to construct complex objects like a `DataFrame` nested inside another `DataFrame`.

```js
let genomic_results = new bioc.DataFrame(
    { 
        logFC: new Float64Array([-1, -2, 1.3, 2.1]),
        pvalue: new Float64Array([0.01, 0.02, 0.001, 1e-8]),
        location: new bioc.DataFrame({
            "chromosome": [ "chrA", "chrB", "chrC", "chrD" ],
            "start": [ 1, 2, 3, 4 ],
            "width": [ 10, 20, 30, 40 ],
            "strand": new Uint8Array([-1, 1, 1, -1 ])
        })
    },
    {
        rowNames: [ "p53", "SNAP25", "MALAT1", "INS" ]
    }
);

let subset = bioc.SLICE(genomic_results, { start: 2, end: 4 });
bioc.LENGTH(subset); 
subset.column("location");
```

Indeed, we can even store an `IRanges` as a column of our `DataFrame`, and all generics on the `DataFrame` will propagate to the column.

```js
let old_location = genomic_results.column("location");
let new_location = new bioc.GRanges(old_location.column("chromosome"),
    new bioc.IRanges(old_location.column("start"), old_location.column("width")),
    { strand: old_location.column("strand") });
genomic_results.$setColumn("location", new_location);

subset = bioc.SLICE(genomic_results, { start: 2, end: 4 });
subset.column("location");
```

## Supported classes and generics

For classes:

|**Javascript**|**R/Bioconductor equivalent**|
|---|---|
| [`DataFrame`](https://ltla.github.io/bioconductor.js/DataFrame.html) | `S4Vectors::DFrame` |
| [`IRanges`](https://ltla.github.io/bioconductor.js/IRanges.html) | `IRanges::IRanges` |
| [`GRanges`](https://ltla.github.io/bioconductor.js/GRanges.html) | `GenomicRanges::GRanges` |
| [`GRanges`](https://ltla.github.io/bioconductor.js/SummarizedExperiment.html) | `SummarizedExperiment::SummarizedExperiment` |

For generics:

|**Javascript**|**R/Bioconductor equivalent**|
|---|---|
| [`LENGTH`](https://ltla.github.io/bioconductor.js/LENGTH.html) | `base::NROW` |
| [`SLICE`](https://ltla.github.io/bioconductor.js/SLICE.html) | `S4Vectors::extractROWS` |
| [`COMBINE`](https://ltla.github.io/bioconductor.js/COMBINE.html) | `S4Vectors::bindROWS` |
| [`CLONE`](https://ltla.github.io/bioconductor.js/CLONE.html) | - |
| [`NUMBER_OF_ROWS`](https://ltla.github.io/bioconductor.js/NUMBER_OF_ROWS.html) | `base::NROW` |
| [`NUMBER_OF_COLUMNS`](https://ltla.github.io/bioconductor.js/NUMBER_OF_COLUMNS.html) | `base::NCOL` |
| [`SLICE_2D`](https://ltla.github.io/bioconductor.js/SLICE_2D.html) | `[` |
| [`COMBINE_ROWS`](https://ltla.github.io/bioconductor.js/COMBINE_ROWS.html) | `S4Vectors::bindROWS` |
| [`COMBINE_COLUMNS`](https://ltla.github.io/bioconductor.js/COMBINE_COLUMNS.html) | `S4Vectors::bindCOLS` |

## Design considerations

We mimic R's S4 generics using methods in Javascript classes.
For example, each supported class defines a `_bioconductor_LENGTH` method to quantify its concept of "length".
The `LENGTH` function will then call this method to obtain a length value for any instance of any supported class.
We prefix this method with `_bioconductor_` to avoid collisions with other properties;
this allows safe monkey patching of third-party classes if they are sufficiently vector-like.

Admittedly, the `LENGTH` function is not really necessary, as users could just call `_bioconductor_LENGTH` directly.
However, the latter is long and unpleasant to type, so we might as well wrap it in something that's easier to remember.
It would also require monkey patching of built-in classes like Arrays and TypedArrays, which is somewhat concerning as it risks interfering with the behavior of other packages.
By defining our own `LENGTH` function, we can safely handle the built-in classes as special cases without modifying their prototypes.

Another difference from R is that Javascript functions can modify objects in place. 
To avoid confusion, all class methods that modify their instances are prefixed with the `$` symbol.
The same applies to all functions that modify their input arguments.
Functions or methods without `$` should not modify their inputs when used with default parameters,
though modifications are allowed if the user provides explicit settings to do so.
