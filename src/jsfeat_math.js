/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(global) {
    "use strict";
    //

    var math = (function() {

        return {
            // looks like latest V8 engine is super smart :)
            // so it turns none of hacks i published here can beat Math class
            // probably next time...
        };

    })();

    global.math = math;

})(jsfeat);