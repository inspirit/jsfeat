/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 */

(function(global) {
    "use strict";
    //

	var gauss_motion_filter = (function () {

	    function gauss_motion_filter(radius, stdev) {
	        if (typeof radius === "undefined") { radius=15; }
	        if (typeof stdev === "undefined") { stdev=-1.0; }

	        this.radius = radius;
	        this.stdev = stdev > 0 ? stdev : Math.sqrt(radius);
	        this.weight = new Float32Array(2 * radius + 1);

	        this.setup(this.radius, this.stdev);
	    };

	    gauss_motion_filter.prototype.setup = function(radius, stdev) {
	    	if (typeof radius === "undefined") { radius=15; }
	        if (typeof stdev === "undefined") { stdev=-1.0; }

	        this.radius = radius;
	        this.stdev = stdev > 0 ? stdev : Math.sqrt(radius);
	        this.weight = new Float32Array(2 * radius + 1);

	        var sum = 0.0;
	        var stdev2 = this.stdev*this.stdev;
	        var i = -this.radius;
		    for (; i <= this.radius; ++i) {
		        sum += this.weight[this.radius + i] = Math.exp(-i*i/(stdev2));
		    }
		    for (i = -this.radius; i <= this.radius; ++i) {
		        this.weight[this.radius + i] /= sum;
		    }
	    };

	    var m33 = new jsfeat.matrix_t(3,3,jsfeat.F32C1_t);

	    gauss_motion_filter.prototype.stabilize = function(idx, motions, from_idx, to_idx) {
		    var cur, curd;
		    var res = m33;
		    var resd=res.data;
		    var sum = 0.0, val=0.0;
		    var iMin = Math.max(idx - this.radius, from_idx)|0;
		    var iMax = Math.min(idx + this.radius, to_idx)|0;
		    var i=iMin;

		    var res0=0.0,res1=0.0,res2=0.0;
		    var res3=0.0,res4=0.0,res5=0.0;
		    var res6=0.0,res7=0.0,res8=0.0;

		    for (; i <= iMax; ++i) {
		    	val = this.weight[this.radius + i - idx];
		    	cur = videostab.get_motion(idx, i, motions);
		    	curd = cur.data;
		        res0 += val * curd[0]; res1 += val * curd[1]; res2 += val * curd[2];
		        res3 += val * curd[3]; res4 += val * curd[4]; res5 += val * curd[5];
		        res6 += val * curd[6]; res7 += val * curd[7]; res8 += val * curd[8];
		        sum += val;
		    }
		    if(sum > 0.0) {
		    	sum = 1.0 / sum;
		    	resd[0]=res0*sum; resd[1]=res1*sum; resd[2]=res2*sum;
		    	resd[3]=res3*sum; resd[4]=res4*sum; resd[5]=res5*sum;
		    	resd[6]=res6*sum; resd[7]=res7*sum; resd[8]=res8*sum;
		    } else {
		    	jsfeat.matmath.identity_3x3(res, 1.0);
		    }
		    return res;
		};

	    return gauss_motion_filter;
	})();

	global.gauss_motion_filter = gauss_motion_filter;

})(videostab);
