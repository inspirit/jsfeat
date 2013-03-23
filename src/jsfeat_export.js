/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(lib) {
    "use strict";

    if (typeof module === "undefined" || typeof module.exports === "undefined") {
        // in a browser, define its namespaces in global
        window.jsfeat = lib;
    } else {
        // in commonjs, or when AMD wrapping has been applied, define its namespaces as exports
        module.exports = lib;
    }
})(jsfeat);
