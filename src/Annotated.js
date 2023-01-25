import * as generics from "./AllGenerics.js";

/**
 * The Annotated class provides a store for arbitrary object-wide metadata.
 * It is intended as a base class for other structures and should not be constructed directly.
 */
export class Annotated {
    /**
     * @param {Object} [options={}] - Optional parameters.
     * @param {Object} [options.metadata={}] - Object containing arbitrary metadata as key-value pairs.
     */
    constructor({ metadata = {} } = {}) {
        this._metadata = metadata;
    }

    /**
     * @return {Object} Object containing arbitrary metadata.
     */
    metadata() {
        return this._metadata;
    }

    /**
     * @param {Object} value - Object containing the metadata.
     *
     * @return {Annotated} Reference to this Annotated object after replacing the metadata.
     */
    $setMetadata(value) {
        this._metadata = value;
        return this;
    }

    _bioconductor_CLONE(output, { deepCopy = true }) {
        output._metadata = generics.CLONE(this._metadata, { deepCopy });
        return;
    }
}
