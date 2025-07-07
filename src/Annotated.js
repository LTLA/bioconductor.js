import * as generics from "./AllGenerics.js";
import * as cutils from "./clone-utils.js";
import * as utils from "./utils.js";
import * as list from "./List.js";

/**
 * The Annotated class provides a store for arbitrary object-wide metadata.
 * It is intended as a base class for other structures and should not be constructed directly.
 *
 * Constructors of Annotated subclasses should be callable with no arguments, possibly creating an empty object with no properties.
 * This will be used by the `_bioconductor_CLONE` method to return an instance of the subclass.
 */
export class Annotated {
    /**
     * @param {Object|Map|Array|List} metadata - Arbitrary list of metadata.
     * An object/Map is converted to a named {@link List}, while an array is converted to an unnamed List.
     */
    constructor(metadata) {
        if (arguments.length == 0) {
            return;
        }

        if (!(metadata instanceof list.List)) {
            metadata = new list.List(metadata);
        }
        this._metadata = metadata;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @return {List} List of arbitrary metadata.
     */
    metadata() {
        return this._metadata;
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    /**
     * @param {Object|Map|Array|List} value - Object containing the metadata.
     * @param {Object} [options={}] - Optional parameters.
     * @param {boolean} [options.inPlace=false] - Whether to mutate this Annotated instance in place.
     * If `false`, a new instance is returned.
     *
     * @return {Annotated} The Annotated object after replacing the metadata.
     * If `inPlace = true`, this is a reference to the current instance, otherwise a new instance is created and returned.
     */
    setMetadata(value, { inPlace = false } = {}) {
        let target = cutils.setterTarget(this, inPlace);

        if (!(value instanceof list.List)) {
            value = new list.List(value);
        }
        target._metadata = value;

        return target;
    }

    $setMetadata(value) {
        return this.setMetadata(value, { inPlace: true });
    }

    /**************************************************************************
     **************************************************************************
     **************************************************************************/

    _bioconductor_CLONE({ deepCopy = true }) {
        let output = new this.constructor;
        output._metadata = cutils.cloneField(this._metadata, deepCopy);
        return output;
    }
}
