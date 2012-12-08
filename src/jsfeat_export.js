/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */
 
(function(lib) {
    "use strict";

    var shim = {};
    if (typeof(exports) === 'undefined') {
        // in a browser, define its namespaces in global
        shim.exports = window;
    } else {    
        // in commonjs, define its namespaces in exports
        shim.exports = exports;
    }

    if(typeof(shim.exports) !== 'undefined') {
        shim.exports.jsfeat = lib;
    }
})(jsfeat);