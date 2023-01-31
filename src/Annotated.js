import * as generics from "./AllGenerics.js";
import * as cutils from "./clone-utils.js";
import * as utils from "./utils.js";

/**
 * The Annotated class provides a store for arbitrary object-wide metadata.
 * It is intended as a base class for other structures and should not be constructed directly.
 */
export class Annotated {
    /**
     * @param {Object|Map} metadata - Object or Map containing arbitrary metadata as key-value pairs.
     */
    constructor(metadata) {
        if (arguments.length == 0) {
            return;
        }

        this._metadata = utils.object2map(metadata);
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {Map} Map containing arbitrary metadata.
     */
    metadata() {
        return this._metadata;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {Object|Map} value - Object containing the metadata.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this Annotated instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {Annotated} The Annotated object after replacing the metadata.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setMetadata(value, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);
        target._metadata = utils.object2map(value);
        return target;
    }

    /**
     * @param {Object} value - Object containing the metadata.
     * @return {Annotated} A reference to this Annotated object.
     */
    $setMetadata(value) {
        return this.setMetadata(value, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_CLONE(output, { deepCopy = true }) {
        output._metadata = cutils.cloneField(this._metadata, deepCopy);
        return;
    }
}
