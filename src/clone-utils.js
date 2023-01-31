import * as generics from "./AllGenerics.js";

export function setterTarget(object, inPlace) {
    return (inPlace ? object : generics.CLONE(object, { deepCopy: false }));
}

export function cloneField(value, deepCopy) {
    return (deepCopy ? generics.CLONE(value) : value);
}
