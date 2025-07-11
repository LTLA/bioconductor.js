import * as vec from "./Vector.js";
import * as gr from "./GRanges.js";
import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as generics from "./AllGenerics.js";

/**
 * A GroupedGRanges object is a collection of groups of genomic ranges, inspired by the `GRangesList` class from the Bioconductor ecosystem.
 * Each group consists of a {@linkplain GRanges} object of arbitrary length, which is most often used to represent a multi-exon gene.
 * The GroupedGRanges can be considered a vector of groups, and defines methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 *
 * Our implementation re-uses Bioconductor's strategy of storing the groups in a single concatenated GRanges.
 * This improves efficiency for large numbers of small GRanges, especially in placeholder objects where all the GRanges are zero-length.
 *
 * Constructors of GroupedGRanges subclasses should be callable with no arguments, possibly creating an empty object with no properties.
 * This will be used by the `_bioconductor_CLONE`, `_bioconductor_SLICE` and `_bioconductor_COMBINE` methods to return an instance of the subclass.
 * 
 * @extends Vector
 */
export class GroupedGRanges extends vec.Vector {
    static #computeStarts(lengths) {
        let starts = new Int32Array(lengths.length);
        let last = 0;
        for (var i = 0; i < lengths.length; i++) {
            starts[i] = last;
            last += lengths[i];
        }
        return { starts: starts, total: last };
    }

    #staged_setGroup = null;

    /**
     * @param {Array|GRanges} ranges - An array of {@linkplain GRanges} objects, where each element represents a group of genomic ranges.
     * All objects should have compatible columns in their {@linkplain Vector#elementMetadata elementMetadata}.
     * 
     * Alternatively, a single GRanges containing a concatenation of ranges from all groups.
     * In this case, `rangeLengths` must be supplied.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?(TypedArray|Array)} [options.rangeLengths=null] - Length of the ranges within each group.
     * This should be coercible to an Int32Array, contain non-negative values, and have a sum equal to the length of `ranges`.
     * Only used if `ranges` is a single {@linkplain GRanges} object, where each group's ranges are assumed to form contiguous intervals along `ranges`.
     * @param {?Array} [options.names=null] - Array of strings of length equal to `start`, containing names for each genomic range.
     * Alternatively `null`, in which case the ranges are assumed to be unnamed.
     * @param {?DataFrame} [options.elementMetadata=null] - A {@linkplain DataFrame} with number of rows equal to the length of `start`, containing arbitrary per-range annotations.
     * Alternatively `null`, in which case a zero-column DataFrame is automatically constructed.
     * @param {Object|Array|Map|List} [options.metadata={}] - Arbitrary metadata, see the {@link Annotated} constructor. 
     */
    constructor(ranges, { rangeLengths = null, names = null, elementMetadata = null, metadata = {} } = {}) {
        if (arguments.length == 0) {
            super();
            return;
        }

        if (ranges.constructor == Array) {
            super(ranges.length, { names, elementMetadata, metadata });
            rangeLengths = new Int32Array(ranges.length);
            for (var i = 0; i < rangeLengths.length; i++) {
                if (!(ranges[i] instanceof gr.GRanges)) {
                    throw new Error("'ranges' must either be a 'GRanges' or an array of 'GRanges'");
                }
                rangeLengths[i] = generics.LENGTH(ranges[i]);
            }
            ranges = generics.COMBINE(ranges);

        } else {
            if (!(ranges instanceof gr.GRanges)) {
                throw new Error("'ranges' must either be a 'GRanges' or an array of 'GRanges'");
            }
            if (rangeLengths == null) {
                throw new Error("'rangeLengths' must be specified when 'ranges' is a 'GRanges'");
            }
            super(rangeLengths.length, { names, elementMetadata, metadata });
            rangeLengths = utils.convertToInt32Array(rangeLengths);
            utils.checkNonNegative(rangeLengths);
        }

        this._ranges = ranges;
        this._rangeLengths = rangeLengths;

        let accumulated = GroupedGRanges.#computeStarts(rangeLengths);
        this._rangeStarts = accumulated.starts;

        if (accumulated.total !== generics.LENGTH(ranges)) {
            throw new Error("sum of 'rangeLengths' must be equal to the length of 'ranges'");
        }
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {GRanges} The concatenated set of ranges across all groups. 
     */
    ranges() {
        this.#flush_staged_setGroup();
        return this._ranges;
    }

    /**
     * @return {Int32Array} The start indices for each group's ranges along the concatenated set of ranges returned by {@linkcode GroupedGRanges#ranges ranges}.
     */
    rangeStarts() {
        this.#flush_staged_setGroup();
        return this._rangeStarts;
    }

    /**
     * @return {Int32Array} The length of each group's ranges along the concatenated set of ranges returned by {@linkcode GroupedGRanges#ranges ranges}.
     */
    rangeLengths() {
        this.#flush_staged_setGroup();
        return this._rangeLengths;
    }

    /**
     * @param {number} i - Index of the group of interest.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.allowView=false] - Whether a view can be created in any internal slicing operations.
     *
     * @return {GRanges} The genomic ranges for group `i`.
     */
    group(i, { allowView = false } = {}) {
        this.#flush_staged_setGroup();
        let s = this._rangeStarts[i];
        return generics.SLICE(this._ranges, { start: s, end: s + this._rangeLengths[i] }, { allowView });
    }

    /**
     * @return {number} Number of groups in this object.
     */
    numberOfGroups() {
        return this._rangeStarts.length;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {GRanges} ranges - Genomic ranges of length equal to the concatenated set of ranges returned by {@linkcode GroupedGRanges#ranges ranges}.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this GroupedGRanges instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {GroupedGRanges} The GroupedGRanges object after modifying the internal ranges.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setRanges(ranges, { inPlace = false } = {}) {
        if (!(ranges instanceof gr.GRanges)) {
            throw new Error("'ranges' must be a 'GRanges'");
        }

        this.#flush_staged_setGroup();
        if (generics.LENGTH(ranges) !== generics.LENGTH(this._ranges)) {
            throw utils.formatLengthError("'ranges'", "number of ranges");
        }

        let target = cutils.setterTarget(this, inPlace);
        target._ranges = ranges;
        return target;
    }

    $setRanges(ranges) {
        return this.setRanges(ranges, { inPlace: true });
    }

    #flush_staged_setGroup() {
        let staged = this.#staged_setGroup;
        if (staged === null) {
            return;
        }

        staged.sort((a, b) => {
            let diff = a[0] - b[0];
            return (diff === 0 ? a[1] - b[1] : diff);
        });

        let counter = 0;
        let accumulated = 0;
        let last_start = 0;
        let more_ranges = [];

        let ngroups = this.numberOfGroups();
        for (var g = 0; g < ngroups; g++) {
            if (counter < staged.length && g == staged[counter][0]) { 
                let current_start = this._rangeStarts[g];
                if (last_start < current_start) {
                    more_ranges.push(generics.SLICE(this._ranges, { start: last_start, end: current_start }));
                }
                last_start = current_start + this._rangeLengths[g];

                let replacement;
                do {
                    replacement = staged[counter][2];
                    counter++;
                } while (counter < staged.length && g == staged[counter][0]);

                more_ranges.push(replacement);
                this._rangeLengths[g] = generics.LENGTH(replacement);
            }

            this._rangeStarts[g] = accumulated;
            accumulated += this._rangeLengths[g];
        }

        let nranges = generics.LENGTH(this._ranges);
        if (last_start < nranges) {
            more_ranges.push(generics.SLICE(this._ranges, { start: last_start, end: nranges }));
        }

        try {
            this._ranges = generics.COMBINE(more_ranges);
        } catch (e) {
            throw new Error("failed to combine staged '$setGroup' operations; " + e.message);
        }

        this.#staged_setGroup = null;
        return;
    }

    /**
     * Multiple consecutive calls to `$setGroup` are not executed immediately.
     * Rather, the operations are staged and executed in batch once the modified GroupedGRanges is used in other methods.
     * This enables efficient setting of individual groups inside a single concatenated {@linkplain GRanges}. 
     *
     * @param {number} i - Index of the group of interest.
     * @param {GRanges} ranges - Genomic ranges for group `i`.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this GroupedGRanges instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {GroupedGRanges} The GroupedGRanges object after setting group `i`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setGroup(i, ranges, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        if (target.#staged_setGroup === null) {
            target.#staged_setGroup = [];
        } else if (!inPlace) {
            target.#staged_setGroup = target.#staged_setGroup.slice();
        }

        if (!inPlace) {
            target._rangeStarts = target._rangeStarts.slice();
            target._rangeLengths = target._rangeLengths.slice();
        }

        let nops = target.#staged_setGroup.length;
        target.#staged_setGroup.push([i, nops, ranges]);
        return target;
    }

    $setGroup(i, ranges) {
        return this.setGroup(i, ranges, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {Object} [options={}] - Optional parameters.
     * @param {?(Array|Set)} [options.restrictToSeqnames=null] - Array or Set containing the sequence names to use in the index.
     * If `null`, all available sequence names are used.
     * @param {?(Array|Set)} [options.restrictToStrand=null] - Array or Set containing the strands to use in the index.
     * If `null`, all available strands are used.
     *
     * @return {GroupedGRangesOverlapIndex} A pre-built index for computing overlaps with other {@linkplain GRanges} instances.
     */
    buildOverlapIndex({ restrictToSeqnames = null, restrictToStrand = null } = {}) {
        this.#flush_staged_setGroup();
        return new GroupedGRangesOverlapIndex(
            this._ranges.buildOverlapIndex({ restrictToSeqnames, restrictToStrand }),
            generics.LENGTH(this._ranges),
            this._rangeStarts,
            this._rangeLengths
        );
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_LENGTH() {
        return this._rangeStarts.length;
    }

    _bioconductor_SLICE(i, { allowView = false } = {}) {
        let output = super._bioconductor_SLICE(i, { allowView });
        this.#flush_staged_setGroup();

        output._rangeLengths = generics.SLICE(this._rangeLengths, i, { allowView });
        let accumulated = GroupedGRanges.#computeStarts(output._rangeLengths);
        output._rangeStarts = accumulated.starts;

        if (i.constructor == Object) {
            // Handle this specially for optimizing allowView = true.
            let s = this._rangeStarts[i.start];
            output._ranges = generics.SLICE(this._ranges, { start: s, end: s + accumulated.total }, { allowView });
        } else {
            let keep = new Int32Array(accumulated.total);

            let counter = 0;
            for (const j of i) {
                let start = this._rangeStarts[j];
                let end = start + this._rangeLengths[j];
                for (var k = start; k < end; k++) {
                    keep[counter] = k;
                    counter++;
                }
            }

            output._ranges = generics.SLICE(this._ranges, keep, { allowView });
        }

        return output;
    }

    _bioconductor_COMBINE(objects) {
        // We need to flush the staged operations in each object.
        this.#flush_staged_setGroup();
        for (const x of objects) {
            x.#flush_staged_setGroup();
        }

        let all_rl = [this._rangeLengths];
        let all_ranges = [this._ranges];
        for (const x of objects) {
            all_rl.push(x._rangeLengths);
            all_ranges.push(x._ranges);
        }

        let output = super._bioconductor_COMBINE(objects);
        output._rangeLengths = generics.COMBINE(all_rl);
        let accumulated = GroupedGRanges.#computeStarts(output._rangeLengths);
        output._rangeStarts = accumulated.starts;
        output._ranges = generics.COMBINE(all_ranges);

        return output;
    }

    _bioconductor_CLONE({ deepCopy = true }) {
        let output = super._bioconductor_CLONE({ deepCopy });

        output.#staged_setGroup = cutils.cloneField(this.#staged_setGroup, deepCopy);
        output._rangeLengths = cutils.cloneField(this._rangeLengths, deepCopy);
        output._rangeStarts = cutils.cloneField(this._rangeStarts, deepCopy);
        output._ranges = cutils.cloneField(this._ranges, deepCopy);

        return output;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {number} [numberOfGroups=0] - Numbe of empty groups to create.
     * @return {GroupedGRanges} A GroupedGRanges object of length equal to `numberOfGroups`,
     * where each group is of zero length.
     */
    static empty(numberOfGroups) {
        let runs = new Int32Array(numberOfGroups);
        runs.fill(0);
        return new GroupedGRanges(gr.GRanges.empty(), { rangeLengths: runs });
    }
}

/**
 * Pre-built index for overlapping {@linkplain GroupedGRanges} objects.
 * This is typically constructed using the {@linkcode GroupedGRanges#buildOverlapIndex GroupedGRanges.buildOverlapIndex} method for a "reference" object,
 * and can be applied to different query GroupedGRanges or {@linkplain GRanges} to identify overlaps with the reference.
 *
 * @hideconstructor
 */
export class GroupedGRangesOverlapIndex {
    constructor(index, fullLength, rangeStarts, rangeLengths) {
        this._index = index;
        this._rangeStarts = rangeStarts;
        this._rangeLengths = rangeLengths;

        let rev_map = new Int32Array(fullLength);
        for (var i = 0; i < rangeStarts.length; i++) {
            let start = rangeStarts[i];
            let end = start + rangeLengths[i];
            for (var s = start; s < end; s++) {
                rev_map[s] = i;
            }
        }
        this._reverseMapping = rev_map;
    }

    /**
     * @param {GroupedGRanges|GRanges} query - The query object, containing ranges to be overlapped with those in the reference GroupedGRanges (that was used to construct this GroupedGRangesOverlapIndex object).
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.ignoreStrand=true] - Whether to ignore differences in strandedness between the ranges in `query` and the reference object.
     *
     * @return {Array} An array of length equal to the number of ranges or groups in `query`,
     * where each element is an array containing the indices of the overlapping ranges in the reference {@linkplain GRanges} object.
     */
    overlap(query, { ignoreStrand = true } = {}) {
        let output = new Array(this._rangeStarts.length);
        let rev_map = this._reverseMapping;

        if (query instanceof GroupedGRanges) {
            let overlaps = this._index.overlap(query._ranges);
            for (var i = 0; i < query._rangeStarts.length; i++) {
                let start = query._rangeStarts[i];
                let end = start + query._rangeLengths[i];

                let results = new Set;
                for (var s = start; s < end; s++) {
                    overlaps[s].forEach(x => results.add(rev_map[x]));
                }
                output[i] = Array.from(results);
            }

        } else {
            let overlaps = this._index.overlap(query);
            for (var i = 0; i < overlaps.length; i++) {
                let results = new Set;
                overlaps[i].forEach(x => results.add(rev_map[x]));
                output[i] = Array.from(results);
            }
        }

        return output;
    }

}
