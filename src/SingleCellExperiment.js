import * as generics from "./AllGenerics.js";
import * as rse from "./RangedSummarizedExperiment.js";
import * as se from "./SummarizedExperiment.js";
import * as utils from "./utils.js";
import * as cutils from "./clone-utils.js";
import * as il from "./InternalList.js";

/**
 * A SingleCellExperiment is a {@linkplain RangedSummarizedExperiment} subclass that contains additional fields for storing reduced dimensions and alternative experiments.
 * It supports the same set of generics as the {@linkplain SummarizedExperiment}.
 *
 * Each reduced dimension instance should have number of rows equal to the number of columns of the SingleCellExperiment.
 * Each instance is expected to provide methods for the following generics:
 *
 * - {@linkcode NUMBER_OF_ROWS}
 * - {@linkcode SLICE_2D}
 * - {@linkcode COMBINE_ROWS}
 * - {@linkcode CLONE}
 *
 * Each alternative experiment should be a {@linkplain SummarizedExperiment} with number of columns equal to that of the SingleCellExperiment.
 *
 * @extends RangedSummarizedExperiment
 */
export class SingleCellExperiment extends rse.RangedSummarizedExperiment {
    /**
     * @param {Object} assays - Object where keys are the assay names and values are multi-dimensional arrays of experimental data.
     * @param {Object} [options={}] - Optional parameters, including those used in the {@linkplain RangedSummarizedExperiment} constructor.
     * @param {?(GRanges|GroupedGRanges)} [options.rowRanges=null] - Genomic ranges corresponding to each row, see the {@linkplain RangedSummarizedExperiment} constructor.
     * @param {Object|Map} [options.reducedDimensions={}] - Object containing named reduced dimensions.
     * Each value should be a 2-dimensional object with number of rows equal to the number of columns of the assays.
     * @param {?Array} [options.reducedDimensionOrder=null] - Array containing the order of the reduced dimensions.
     * This should have the same values as the keys of `reducedDimensions`, and defaults to those keys if `null`.
     * @param {Object|Map} [options.alternativeExperiments={}] - Object containing named alternative experiments.
     * Each value should be a 2-dimensional object with number of columns equal to that of the assays.
     * @param {?Array} [options.alternativeExperimentOrder=null] - Array containing the order of the alternative experiments.
     * This should have the same values as the keys of `alternativeExperiments`, and defaults to those keys if `null`.
     * @param {?string} [options.mainExperimentName=null] - Main experiment name, possibly `null` if the main experiment is unnamed.
     */
    constructor(assays, options={}) {
        if (arguments.length == 0) {
            super();
            return;
        }

        let { reducedDimensions = {}, reducedDimensionOrder = null, alternativeExperiments = {}, alternativeExperimentOrder = null, rowRanges = null, mainExperimentName = null } = options;
        super(assays, rowRanges, options);
        let ncols = this.numberOfColumns();

        try {
            this._reducedDimensions = new il.InternalList(reducedDimensions, reducedDimensionOrder);
        } catch (e) {
            throw new Error("failed to initialize reduced dimension list for this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        for (const k of this._reducedDimensions.names()) {
            let v = this._reducedDimensions.entry(k);
            if (generics.NUMBER_OF_ROWS(v) !== ncols) {
                throw new Error("number of rows for reduced dimension '" + k + "' is not equal to number of columns for this " + this.constructor.className);
            }
        }

        try {
            this._alternativeExperiments = new il.InternalList(alternativeExperiments, alternativeExperimentOrder);
        } catch (e) {
            throw new Error("failed to initialize alternative experiment list for this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        for (const k of this._alternativeExperiments.names()) {
            let v = this._alternativeExperiments.entry(k);
            if (!(v instanceof se.SummarizedExperiment)) {
                throw new Error("alternative experiment '" + k + "' is not a SummarizedExperiment");
            }
            if (v.numberOfColumns(v) !== ncols) {
                throw new Error("number of columns for alternative experiment '" + k + "' is not equal to number of columns for this " + this.constructor.className);
            }
        }

        this._mainExperimentName = mainExperimentName;
        return;
    }

    static className = "SingleCellExperiment";

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {Array} Array of strings containing the (ordered) names of the reduced dimensions.
     */
    reducedDimensionNames() {
        return this._reducedDimensions.names();
    }

    /**
     * @param {string|number} i - Reduced dimension to retrieve, either by name or index.
     * @return {*} The contents of reduced dimension `i` as an multi-dimensional array-like object.
     */
    reducedDimension(i) {
        let output;
        try {
            output = this._reducedDimensions.entry(i);
        } catch (e) {
            throw new Error("failed to retrieve the specified reduced dimension from this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return output;
    }

    /**
     * @return {Map} Map where keys are the names and values are the reduced dimensions.
     */
    reducedDimensions() {
        return this._reducedDimensions.entries();
    }

    /**
     * @return {Array} Array of strings containing the (ordered) names of the alternative experiments.
     */
    alternativeExperimentNames() {
        return this._alternativeExperiments.names();
    }

    /**
     * @param {string|number} i - Alternative experiment to retrieve, either by name or index.
     * @return {SummarizedExperiment} The specified alternative experiment `i`. 
     */
    alternativeExperiment(i) {
        let output;
        try {
            output = this._alternativeExperiments.entry(i);
        } catch (e) {
            throw new Error("failed to retrieve the specified alternative experiment from this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return output;
    }

    /**
     * @return {Map} Map where keys are the names and values are the alternative experiments.
     */
    alternativeExperiments() {
        return this._alternativeExperiments.entries();
    }

    /**
     * @return {?string} The name of the main experiment, possibly `null` if this is unnamed.
     */
    mainExperimentName() {
        return this._mainExperimentName;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {string|number} i - Identity of the reduced dimension to remove, either by name or index.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SingleCellExperiment} The SingleCellExperiment after removing the specified assay.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    removeReducedDimension(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._reducedDimensions = target._reducedDimensions.delete(i, { inPlace });
        } catch (e) {
            throw new Error("failed to remove the specified reduced dimension from this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $removeReducedDimension(i) {
        return this.removeReducedDimension(i, { inPlace: true });
    }

    /**
     * @param {string|number} i - Identity of the reduced dimension to add, either by name or index.
     * - If `i` is a number, the reduced dimension at the specified index is replaced.
     *   `i` should be non-negative and less than the number of reduced dimensions.
     * - If `i` is a string, any reduced dimension with the same name is replaced.
     *   If no such reduced dimension exists, a new reduced dimension is appended to the list of reduced dimensions.
     * @param {*} value - Multi-dimensional array-like object to set/add as the reduced dimension.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SingleCellExperiment} The SingleCellExperiment with modified reduced dimensions.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setReducedDimension(i, value, { inPlace = false } = {}) {
        if (generics.NUMBER_OF_ROWS(value) != this.numberOfColumns()) {
            throw new Error("number of rows of 'value' should be the same as the number of columns of this SingleCellExperiment");
        }
        let target = cutils.setterTarget(this, inPlace);
        target._reducedDimensions = target._reducedDimensions.set(i, value, { inPlace });
        return target;
    }

    $setReducedDimension(i, value) {
        return this.setReducedDimension(i, value, { inPlace: true });
    }

    /**
     * @param {Array} names - Array of strings containing the reduced dimension names.
     * This should be of the same length as the number of reduced dimensions and contain unique values.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SingleCellExperiment} The SingleCellExperiment with modified reduced dimension names.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setReducedDimensionNames(names, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._reducedDimensions = target._reducedDimensions.setNames(names, { inPlace });
        } catch (e) {
            throw new Error("failed to set the reduced dimension names for this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $setReducedDimensionNames(names) {
        return this.setReducedDimensionNames(names, { inPlace: true });
    }

    /**
     * @param {Object|Map} value - Object containing zero, one or more multi-dimensional array-like objects in the values.
     * Each value should be a 2-dimensional object with number of rows equal to the number of columns in this SingleCellExperiment.
     * Keys are reduced dimension names, each of which should be present in `order`.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     * @param {Array|boolean} [options.newOrder=false] - Whether to replace the order of reduced dimensions with the order of keys in `value`.
     * If `false`, the existing order in {@linkcode SingleCellExperiment#reducedDimensionNames reducedDimensionNames} is used.
     * If an array is provided, this is used as the order.
     * If `null`, this has the same effect as `true`.
     *
     * @return {SingleCellExperiment} A SingleCellExperiment with the new reduced dimensions.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setReducedDimensions(value, { inPlace = false, newOrder = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);

        if (newOrder === false) {
            newOrder = target._reducedDimensions.names();
        } else if (newOrder == true) {
            newOrder = null;
        }
        try {
            target._reducedDimensions = new il.InternalList(value, newOrder);
        } catch (e) {
            throw new Error("failed to replace reduced dimension list for this SingleCellExperiment; " + e.message, { cause: e });
        }

        let sce_nc = target.numberOfColumns();
        for (const k of target._reducedDimensions.names()) {
            let current = target._reducedDimensions.entry(k);
            let nr = generics.NUMBER_OF_ROWS(current);
            if (nr !== sce_nc) {
                throw new Error("mismatch in the number of rows for reduced dimension '" + k + "' compared to the number of columns in the SingleCellExperiment");
            }
        }

        return target;
    }

    $setReducedDimensions(value, { newOrder = false } = {}) {
        return this.setReducedDimensions(value, { inPlace: true, newOrder });
    }

    /**
     * @param {Array} i - Array of strings or indices specifying the reduced dimensions to retain in the slice.
     * This should refer to unique reduced dimension names.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SingleCellExperiment} The SingleCellExperiment with sliced reduced dimensions.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    sliceReducedDimensions(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._reducedDimensions = this._reducedDimensions.slice(i, { inPlace });
        } catch (e) {
            throw new Error("failed to slice the reduced dimensions for this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $sliceReducedDimensions(i) {
        return this.sliceReducedDimensions(i, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {string|number} i - Identity of the reduced dimension to remove, either by name or index.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SingleCellExperiment} The SingleCellExperiment after removing the specified assay.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    removeAlternativeExperiment(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._alternativeExperiments = target._alternativeExperiments.delete(i, { inPlace });
        } catch (e) {
            throw new Error("failed to remove the specified alternative experiment from this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $removeAlternativeExperiment(i) {
        return this.removeAlternativeExperiment(i, { inPlace: true });;
    }

    /**
     * @param {string|number} i - Identity of the alternative experiment to add, either by name or index.
     * - If `i` is a number, the alternative experiment at the specified index is replaced.
     *   `i` should be non-negative and less than the number of alternative experiments.
     * - If `i` is a string, any alternative experiment with the same name is replaced.
     *   If no such alternative experiment exists, a new alternative experiment is appended to the list of alternative experiments.
     * @param {SummarizedExperiment} value - A SummarizedExperiment to set/add as the alternative experiment.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SingleCellExperiment} The SingleCellExperiment with modified alternative experiments.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setAlternativeExperiment(i, value, { inPlace = false } = {}) {
        if (!(value instanceof se.SummarizedExperiment) || generics.NUMBER_OF_COLUMNS(value) != this.numberOfColumns()) {
            throw new Error("'value' should be a SummarizedExperiment with the same number of columns as this SingleCellExperiment");
        }
        let target = cutils.setterTarget(this, inPlace);
        target._alternativeExperiments = target._alternativeExperiments.set(i, value, { inPlace });
        return target;
    }

    $setAlternativeExperiment(i, value) {
        return this.setAlternativeExperiment(i, value, { inPlace: true });
    }

    /**
     * @param {Array} names - Array of strings containing the alternative experiment names.
     * This should be of the same length as the number of alternative experiments and contain unique values.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SingleCellExperiment} The SingleCellExperiment with modified alternative experiment names.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setAlternativeExperimentNames(names, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._alternativeExperiments = target._alternativeExperiments.setNames(names, { inPlace });
        } catch (e) {
            throw new Error("failed to set the alternative experiment names for this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $setAlternativeExperimentNames(names) {
        return this.setAlternativeExperimentNames(names, { inPlace: true });
    }

    /**
     * @param {Object|Map} value - Object containing zero, one or more {@link SummarizedExperiment} objects in the values.
     * Each value should have the same number of columns as this SingleCellExperiment.
     * Keys are alternative experiment names, each of which should be present in `order`.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     * @param {Array|boolean} [options.newOrder=false] - Whether to replace the order of alternative experiments with the order of keys in `value`.
     * If `false`, the existing order in {@linkcode SingleCellExperiment#alternativeExperimentNames alternativeExperimentNames} is used.
     * If an array is provided, this is used as the order.
     * If `null`, this has the same effect as `true`.
     *
     * @return {SingleCellExperiment} A SingleCellExperiment with the new alternative experiments.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setAlternativeExperiments(value, { inPlace = false, newOrder = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);

        if (newOrder === false) {
            newOrder = target._alternativeExperiments.names();
        } else if (newOrder == true) {
            newOrder = null;
        }
        try {
            target._alternativeExperiments = new il.InternalList(value, newOrder);
        } catch (e) {
            throw new Error("failed to replace alternative experiment list for this SingleCellExperiment; " + e.message, { cause: e });
        }

        let sce_nc = target.numberOfColumns();
        for (const k of target._alternativeExperiments.names()) {
            let current = target._alternativeExperiments.entry(k);
            let nr = generics.NUMBER_OF_COLUMNS(current);
            if (nr !== sce_nc) {
                throw new Error("mismatch in the number of columns for alternative experiment '" + k + "' compared to the SingleCellExperiment");
            }
        }

        return target;
    }

    $setAlternativeExperiments(value, { newOrder = false } = {}) {
        return this.setAlternativeExperiments(value, { inPlace: true, newOrder });
    }

    /**
     * @param {Array} i - Array of strings or indices specifying the alternative experiments to retain in the slice.
     * This should refer to unique alternative experiment names.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {SingleCellExperiment} The SingleCellExperiment with sliced alternative experiments.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    sliceAlternativeExperiments(i, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        try {
            target._alternativeExperiments = this._alternativeExperiments.slice(i, { inPlace });
        } catch (e) {
            throw new Error("failed to slice the alternative experiments for this " + this.constructor.className + "; " + e.message, { cause: e });
        }
        return target;
    }

    $sliceAlternativeExperiments(i) {
        return this.sliceAlternativeExperiments(i, { inPlace: true });
    }

    /**
     * @return {?string} name - The name of the main experiment, possibly `null` if this is unnamed.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this SingleCellExperiment instance in place.
     * If `false`, a new instance is returned.
     * @return {SingleCellExperiment} A SingleCellExperiment with a new main experiment name.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setMainExperimentName(name, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        target._mainExperimentName = name;
        return target;
    }

    $setMainExperimentName(name) {
        return this.setMainExperimentName(name, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_SLICE_2D(rows, columns, { allowView = false }) {
        let output = super._bioconductor_SLICE_2D(rows, columns, { allowView });

        if (columns !== null) {
            output._reducedDimensions = this._reducedDimensions.apply(v => generics.SLICE_2D(v, columns, null, { allowView }));
            output._alternativeExperiments = this._alternativeExperiments.apply(v => generics.SLICE_2D(v, null, columns, { allowView }));
        } else {
            output._reducedDimensions = this._reducedDimensions;
            output._alternativeExperiments = this._alternativeExperiments;
        }

        output._mainExperimentName = this._mainExperimentName;
        return output;
    }

    _bioconductor_COMBINE_ROWS(objects) {
        let output = super._bioconductor_COMBINE_ROWS(objects);
        output._reducedDimensions = this._reducedDimensions;
        output._alternativeExperiments = this._alternativeExperiments;
        output._mainExperimentName = this._mainExperimentName;
        return output;
    }

    _bioconductor_COMBINE_COLUMNS(objects) {
        let output = super._bioconductor_COMBINE_COLUMNS(objects);

        let all_rd = [this._reducedDimensions];
        let all_ae = [this._alternativeExperiments];
        for (const x of objects) {
            all_rd.push(x._reducedDimensions);
            all_ae.push(x._alternativeExperiments);
        }

        try {
            output._reducedDimensions = il.InternalList.parallelCombine(all_rd, generics.COMBINE_ROWS);
        } catch (e) {
            throw new Error("failed to combine reduced dimensions for " + this.constructor.className + " objects; " + e.message, { cause: e });
        }

        try {
            output._alternativeExperiments = il.InternalList.parallelCombine(all_ae, generics.COMBINE_COLUMNS);
        } catch (e) {
            throw new Error("failed to combine alternative experiments for " + this.constructor.className + " objects; " + e.message, { cause: e });
        }

        output._mainExperimentName = this._mainExperimentName;
        return output;
    }

    _bioconductor_CLONE({ deepCopy }) {
        let output = super._bioconductor_CLONE({ deepCopy });
        output._reducedDimensions = cutils.cloneField(this._reducedDimensions, deepCopy);
        output._alternativeExperiments = cutils.cloneField(this._alternativeExperiments, deepCopy);
        output._mainExperimentName = this._mainExperimentName;
        return output;
    }
}
