;(function (window, document) {
    'use strict';
    
    window.GumWrapper = function(elements, success, error) {
        // Define our error message
        function sendError(message) {
            if (error) {
                var e = new Error();
                e.message = message;
                error(e);
            } else {
                console.error(message);
            }
        }
        
        // Try to play the media stream
        function play() {
            var video = document.getElementById(elements.video);
            if (!video) {
                sendError('Unable to find the video element.');
                return;
            }

            function successCallback(stream) {
                // Set the source of the video element with the stream from the camera
                if (video.mozSrcObject !== undefined) {
                    video.mozSrcObject = stream;
                } else {
                    video.src = (window.URL && window.URL.createObjectURL(stream)) || stream;
                }
                video.play();
            }

            function errorCallback(error) {
                sendError('Unable to get webcam stream.');
            }

            navigator.getUserMedia = navigator.getUserMedia || navigator.webkitGetUserMedia || navigator.mozGetUserMedia || navigator.msGetUserMedia;
            window.URL = window.URL || window.webkitURL || window.mozURL || window.msURL;

            // Call the getUserMedia method with our callback functions
            if (navigator.getUserMedia) {
                navigator.getUserMedia({video: true}, successCallback, errorCallback);
            } else {
                sendError('Native web camera streaming (getUserMedia) not supported in this browser.');
            }
            
            // Check the video dimensions when data has loaded
            video.addEventListener('loadeddata', function() {
                if (video.videoWidth > 0 && video.videoHeight > 0) {
                    if (success) success(video);
                } else {
                    sendError('Unable to play video stream. Is webcam working?');
                }
            }, false);
        }

        return {
            play: play
        };
    };
})(window, document);
