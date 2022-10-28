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

bioc.COMBINE(results, [more_results]);
```

The generics allow us to manipulate complex objects like a nested `DataFrame`:

```js
let results = new bioc.DataFrame(
    { 
        logFC: new Float64Array([-1, -2, 1.3, 2.1]),
        pvalue: new Float64Array([0.01, 0.02, 0.001, 1e-8]),
        location: new bioc.DataFrame(
            "chromosome": [ "chrA", "chrB", "chrC", "chrD" ],
            "start": [ 1, 2, 3, 4 ],
            "end": [ 10, 20, 30, 40 ],
            "strand": new Uint8Array([-1, 1, 1, -1 ])
        )
    },
    {
        rowNames: [ "p53", "SNAP25", "MALAT1", "INS" ]
    }
);

let subset = bioc.SLICE(results, { start: 2, end: 4 });
bioc.LENGTH(subset); 
subset.column("location");
```

See the [reference documentation](https://ltla.github.io/bioconductor.js) for more details.

## Supported classes and generics

For classes:

|**Javascript**|**R/Bioconductor equivalent**|
|---|---|
| `DataFrame` | `S4Vectors::DFrame` |

For generics:

|**Javascript**|**R/Bioconductor equivalent**|
|---|---|
| `LENGTH` | `base::NROW` |
| `SLICE` | `S4Vectors::extractROWS` |
| `COMBINE` | `S4Vectors::bindROWS` |
| `CLONE` | - |
