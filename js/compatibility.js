/**
* this code is from all around the web :)
* if u want to put some credits u are welcome!
*/
var compatibility = (function() {
        var lastTime = 0,
        isLittleEndian = true,

        URL = window.URL || window.webkitURL,

        requestAnimationFrame = function(callback, element) {
            var requestAnimationFrame =
                window.requestAnimationFrame        || 
                window.webkitRequestAnimationFrame  || 
                window.mozRequestAnimationFrame     || 
                window.oRequestAnimationFrame       ||
                window.msRequestAnimationFrame      ||
                function(callback, element) {
                    var currTime = new Date().getTime();
                    var timeToCall = Math.max(0, 16 - (currTime - lastTime));
                    var id = window.setTimeout(function() {
                        callback(currTime + timeToCall);
                    }, timeToCall);
                    lastTime = currTime + timeToCall;
                    return id;
                };

            return requestAnimationFrame.call(window, callback, element);
        },

        cancelAnimationFrame = function(id) {
            var cancelAnimationFrame = window.cancelAnimationFrame ||
                                        function(id) {
                                            clearTimeout(id);
                                        };
            return cancelAnimationFrame.call(window, id);
        },

        getUserMedia = function(options, success, error) {
            var getUserMedia =
                window.navigator.getUserMedia ||
                window.navigator.mozGetUserMedia ||
                window.navigator.webkitGetUserMedia ||
                window.navigator.msGetUserMedia ||
                function(options, success, error) {
                    error();
                };

            return getUserMedia.call(window.navigator, options, success, error);
        },

        detectEndian = function() {
            var buf = new ArrayBuffer(8);
            var data = new Uint32Array(buf);
            data[0] = 0xff000000;
            isLittleEndian = true;
            if (buf[0] === 0xff) {
                isLittleEndian = false;
            }
            return isLittleEndian;
        };

    return {
        URL: URL,
        requestAnimationFrame: requestAnimationFrame,
        cancelAnimationFrame: cancelAnimationFrame,
        getUserMedia: getUserMedia,
        detectEndian: detectEndian,
        isLittleEndian: isLittleEndian
    };
})();