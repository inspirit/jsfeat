/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(lib) {
    "use strict";

    if (!module || !module.exports) {
        // in a browser, define its namespaces in global
        window.jsfeat = lib;
    } else {
        // in commonjs, or when AMD wrapping has been applied, define its namespaces as exports
        module.exports = lib;
    }
})(jsfeat);
