/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 * The implementation is derived from OpenCV
 */


var videostab = videostab || { version: '0' };


self.console = self.console || {

	info: function () {},
	log: function () {},
	debug: function () {},
	warn: function () {},
	error: function () {}

};

videostab.MM_AFFINE = 0;
videostab.MM_HOMOGRAPHY = 1;

(function(global) {
    "use strict";
    //

	var size_t = (function () {
	    function size_t(w,h) {
	        if (typeof w === "undefined") { w=0; }
	        if (typeof h === "undefined") { h=0; }

	        this.width = w;
	        this.height = h;
	    }
	    return size_t;
	})();

	// RING BUFFER ACCESS
	var get_at = (function () {
        return function(idx, len, items) {
			if( idx < 0 ) idx -= (((idx-len+1)/len)*len);
			if( idx >= len ) idx %= len;
		    return items[idx|0];
		}
    })();
    var get_ring_ind = (function () {
        return function(idx, len) {
			if( idx < 0 ) idx -= (((idx-len+1)/len)*len);
			if( idx >= len ) idx %= len;
		    return idx|0;
		}
    })();

    var __im33 = new jsfeat.matrix_t(3,3,jsfeat.F32C1_t);
    var get_motion = (function () {
        return function(from, to, motions) {
		    var M = __im33;
		    var i=0;
            var n = motions.length;

		    jsfeat.matmath.identity_3x3(M, 1.0);

		    if (to > from) {
		        for (i = from; i < to; ++i) {
		        	jsfeat.matmath.multiply_3x3(M, get_at(i, n, motions), M);
		        }
		    } else if (from > to) {
		        for (i = to; i < from; ++i) {
		            jsfeat.matmath.multiply_3x3(M, get_at(i, n, motions), M);
		        }
		        jsfeat.matmath.invert_3x3(M, M);
		    }
		    return M;
		}
    })();

	global.size_t = size_t;
	global.get_at = get_at;
	global.get_ring_ind = get_ring_ind;
	global.get_motion = get_motion;

})(videostab);

