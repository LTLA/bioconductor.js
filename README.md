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
| [`DataFrame`](https://ltla.github.io/bioconductor.js/DataFrame.html) | `S4Vectors::DFrame` |

For generics:

|**Javascript**|**R/Bioconductor equivalent**|
|---|---|
| [`LENGTH`](https://ltla.github.io/bioconductor.js/LENGTH.html) | `base::NROW` |
| [`SLICE`](https://ltla.github.io/bioconductor.js/SLICE.html) | `S4Vectors::extractROWS` |
| [`COMBINE`](https://ltla.github.io/bioconductor.js/COMBINE.html) | `S4Vectors::bindROWS` |
| [`CLONE`](https://ltla.github.io/bioconductor.js/CLONE.html) | - |

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
