import * as generics from "./AllGenerics.js";
import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as ir from "./IRanges.js";
import * as vec from "./Vector.js";
import * as olap from "./overlap-utils.js";

/**
 * A GRanges object is a collection of genomic ranges, inspired by the class of the same name from the Bioconductor ecosystem.
 * Each range consists of a sequence name, a start position on that sequence, and a width.
 * Each range may also be associated with arbitrary range-level metadata in a {@linkplain DataFrame}.
 * The GRanges defines methods for the following generics:
 *
 * - {@linkcode LENGTH}
 * - {@linkcode SLICE}
 * - {@linkcode COMBINE}
 * - {@linkcode CLONE}
 *
 * @extends Vector
 */
export class GRanges extends vec.Vector {
    static #convertToInt8Array(x) {
        if (x instanceof Int8Array) {
            return x;
        } else {
            return new Int8Array(x);
        }
    }

    static #checkStrandedness(strand) {
        for (const y of strand) {
            if (y < -1 || y > 1) {
                throw new Error("'strand' must be -1, 0 or 1");
            }
        }
    }

    /**
     * @param {Array} seqnames - Array of strings containing the sequence names for each genomic range.
     * @param {IRanges} ranges - Position and width of the range on its specified sequence.
     * This should have the same length as `seqnames`.
     * @param {Object} [options={}] - Optional parameters.
     * @param {?(Array|TypedArray)} [options.strand=null] - Array containing the strandedness of each genomic range.
     * This should be 0 (any strand), 1 (forward strand) or -1 (reverse strand).
     * If `null`, this is assumed to be 0 for all genomic ranges.
     * @param {?Array} [options.names=null] - Array of strings of length equal to `start`, containing names for each genomic range.
     * Alternatively `null`, in which case the ranges are assumed to be unnamed.
     * @param {?DataFrame} [options.elementMetadata=null] - A {@linkplain DataFrame} with number of rows equal to the length of `start`, containing arbitrary per-range annotations.
     * Alternatively `null`, in which case a zero-column DataFrame is automatically constructed.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor(seqnames, ranges, { strand = null, names = null, elementMetadata = null, metadata = {} } = {}) {
        if (arguments.length == 0) {
            super();
            return;
        }

        super(seqnames.length, { names, elementMetadata, metadata });

        utils.checkStringArray(seqnames, "seqnames");
        this._seqnames = seqnames;

        let n = seqnames.length;
        if (n !== generics.LENGTH(ranges)) {
            throw utils.formatLengthError("'ranges'", "'seqnames'");
        }
        this._ranges = ranges;

        if (strand !== null) {
            if (n !== strand.length) {
                throw utils.formatLengthError("'strand'", "'seqnames'");
            }
            strand = GRanges.#convertToInt8Array(strand);
            GRanges.#checkStrandedness(strand);
        } else {
            strand = new Int8Array(n);
            strand.fill(0);
        }
        this._strand = strand;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {Int32Array} Array of integers containing the start position for each genomic range.
     */
    start() {
        return this._ranges.start();
    }

    /**
     * @return {Int32Array} Array of integers containing the end position (specifically, one-past-the-end) for each genomic range.
     */
    end() {
        return this._ranges.end();
    }

    /**
     * @return {Int32Array} Array of integers containing the width of each genomic range.
     */
    width() {
        return this._ranges.width();
    }

    /**
     * @return {Array} Array of strings containing the sequence name for each genomic range.
     */
    seqnames() {
        return this._seqnames;
    }

    /**
     * @return {IRanges} Start positions and widths for all ranges on their specified sequence names.
     */
    ranges() {
        return this._ranges;
    }

    /**
     * @return {Int8Array} Array containing the strandedness for each genomic range - 0 (any strand), 1 (forward strand) or -1 (reverse strand).
     */
    strand() {
        return this._strand;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {Array} seqnames - Array of strings containing the sequence names for each genomic range.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this GRanges instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {GRanges} The GRanges object after setting the sequence names to `seqnames`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setSeqnames(seqnames, { inPlace = false } = {}) {
        utils.checkNamesArray(seqnames, "replacement 'seqnames'", generics.LENGTH(this), "'LENGTH(<GRanges>)'");
        let target = cutils.setterTarget(this, inPlace); 
        target._seqnames = seqnames;
        return target;
    }

    /**
     * @param {Array} seqnames - Array of strings containing the sequence names for each genomic range.
     * @return {GRanges} A reference to this GRanges object after setting the sequence names to `seqnames`.
     */
    $setSeqnames(seqnames) {
        return this.setSeqnames(seqnames, { inPlace: true });
    }

    /**
     * @param {IRanges} ranges - Start positions and widths for each genomic range.
     * This should have length equal to the number of ranges. 
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this GRanges instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {GRanges} The GRanges object after setting the ranges to `ranges`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setRanges(ranges, { inPlace = false } = {}) {
        if (!(ranges instanceof ir.IRanges)) {
            throw new Error("'ranges' should be an IRanges object");
        }

        if (generics.LENGTH(ranges) !== generics.LENGTH(this._ranges)) {
            throw utils.formatLengthError("replacement 'ranges'", "'LENGTH(<GRanges>)'");
        }

        let target = cutils.setterTarget(this, inPlace); 
        target._ranges = ranges;
        return target;
    }

    /**
     * @param {IRanges} ranges - Start positions and widths for each genomic range.
     * This should have length equal to the number of ranges. 
     * @return {GRanges} A reference to this GRanges object after setting the ranges to `ranges`.
     */
    $setRanges(ranges) {
        return this.setRanges(ranges, { inPlace: true });
    }

    /**
     * @param {Array|TypedArray} strand - Array of strands for each genomic range.
     * This should have length equal to the number of ranges. 
     * Entries may be 0 (any strand), 1 (forward strand) or -1 (reverse strand).
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this GRanges instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {GRanges} The GRanges object after setting the strands to `strand`.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setStrand(strand, { inPlace = false }) {
        if (this._strand.length !== strand.length) {
            throw utils.formatLengthError("'strand'", "'seqnames'");
        }
        strand = GRanges.#convertToInt8Array(strand);
        GRanges.#checkStrandedness(strand);

        let target = cutils.setterTarget(this, inPlace); 
        target._strand = strand;
        return target;
    }

    /**
     * @param {Array|TypedArray} strand - Array of strands for each genomic range.
     * This should have length equal to the number of ranges. 
     * Entries may be 0 (any strand), 1 (forward strand) or -1 (reverse strand).
     *
     * @return {GRanges} A reference to this GRanges object after setting the strands to `strand`.
     */
    $setStrand(strand) {
        return this.setStrand(strand, { inPlace: true });
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
     * @return {GRangesOverlapIndex} A pre-built index for computing overlaps with other {@linkplain GRanges} instances.
     */
    buildOverlapIndex({ restrictToSeqnames = null, restrictToStrand = null } = {}) {
        let indices = utils.createSequence(generics.LENGTH(this));
        let by_seqname = generics.SPLIT(indices, this._seqnames);
        let starts = this.start();
        let ends = this.end();

        if (restrictToSeqnames !== null && restrictToSeqnames instanceof Array) {
            restrictToSeqnames = new Set(restrictToSeqnames);
        }
        if (restrictToStrand !== null && restrictToStrand instanceof Array) {
            restrictToStrand = new Set(restrictToStrand);
        }

        for (const name of Object.keys(by_seqname)) {
            if (restrictToSeqnames !== null && !restrictToSeqnames.has(name)) {
                delete by_seqname[name];
                continue;
            }
            let seqname_indices = by_seqname[name];
            let seqname_strand = generics.SLICE(this._strand, seqname_indices);
            let by_strand = generics.SPLIT(seqname_indices, seqname_strand);

            for (const str of Object.keys(by_strand)) {
                if (restrictToStrand !== null && !restrictToStrand.has(Number(str))) {
                    delete by_strand[str];
                    continue;
                }
                let str_indices = by_strand[str];
                by_strand[str] = olap.buildIntervalTree(starts, ends, { slice: str_indices });
            }
            by_seqname[name] = by_strand;
        }

        return new GRangesOverlapIndex(by_seqname);
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_LENGTH() {
        return this._seqnames.length;
    }

    _bioconductor_SLICE(output, i, { allowView = false }) {
        super._bioconductor_SLICE(output, i, { allowView });
        output._seqnames = generics.SLICE(this._seqnames, i, { allowView });
        output._ranges = generics.SLICE(this._ranges, i, { allowView });
        output._strand = generics.SLICE(this._strand, i, { allowView });
        return;
    }

    _bioconductor_COMBINE(output, objects) {
        super._bioconductor_COMBINE(output, objects);

        let all_sn = [];
        let all_rr = [];
        let all_st = [];
        for (const x of objects) {
            all_sn.push(x._seqnames);
            all_rr.push(x._ranges);
            all_st.push(x._strand);
        }

        output._seqnames = generics.COMBINE(all_sn);
        output._ranges = generics.COMBINE(all_rr);
        output._strand = generics.COMBINE(all_st);
        return;
    }

    _bioconductor_CLONE(output, { deepCopy = true }) {
        super._bioconductor_CLONE(output, { deepCopy });
        output._seqnames = cutils.cloneField(this._seqnames, deepCopy);
        output._ranges = cutils.cloneField(this._ranges, deepCopy);
        output._strand = cutils.cloneField(this._strand, deepCopy);
        return;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {GRanges} A zero-length GRanges object.
     */
    static empty() {
        return new GRanges([], ir.IRanges.empty());
    }
}

/**
 * Pre-built index for overlapping {@linkplain GRanges} objects.
 * This is typically constructed using the {@linkcode GRanges#buildOverlapIndex GRanges.buildOverlapIndex} method for a "reference" object,
 * and can be applied to different query GRanges to identify overlaps with the reference.
 *
 * @hideconstructor
 */
export class GRangesOverlapIndex {
    constructor(index) {
        this._index = index;
    }

    /**
     * @param {GRanges} query - The query object, containing ranges to be overlapped with those in the reference GRanges (that was used to construct this GRangesOverlapIndex object).
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.ignoreStrand=true] - Whether to ignore differences in strandedness between the ranges in `query` and the reference object.
     *
     * @return {Array} An array of length equal to the number of ranges in `query`,
     * where each element is an array containing the indices of the overlapping ranges in the reference {@linkplain GRanges} object.
     */
    overlap(query, { ignoreStrand = true } = {}) {
        let n = generics.LENGTH(query);
        let results = new Array(n);
        let starts = query.start();
        let ends = query.end();

        for (var i = 0; i < n; i++) {
            results[i] = [];
            let my_results = results[i];

            let name = query._seqnames[i];
            if (!(name in this._index)) {
                continue;
            }
            let seq_index = this._index[name];

            let strand = query._strand[i];
            let allowed_strands;
            if (ignoreStrand || strand == 0) {
                allowed_strands = Object.keys(seq_index);
            } else {
                let sstr = String(strand);
                if (!(sstr in seq_index)) {
                    continue;
                }
                allowed_strands = [ sstr ];
            }

            let start = starts[i];
            let end = ends[i];
            for (const str of allowed_strands) {
                let str_results = olap.queryIntervalTree(start, end, seq_index[str]);
                str_results.forEach(x => my_results.push(x));
            }
        }

        return results;
    }
}
