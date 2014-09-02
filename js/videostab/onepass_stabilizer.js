/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 */

(function(global) {
    "use strict";
    //

    // utility
    var areaSign = function(a, b, c) {
		var bax = b.x-a.x, bay = b.y-a.y;
		var cax = c.x-a.x, cay = c.y-a.y;
	    var area = bax*cay - bay*cax;
	    if (area < -1e-5) return -1;
	    if (area > 1e-5) return 1;
	    return 0;
	}


	var segmentsIntersect = function(a, b, c, d) {
	    return areaSign(a, b, c) * areaSign(a, b, d) < 0 &&
	           areaSign(c, d, a) * areaSign(c, d, b) < 0;
	}


	// Checks if rect a (with sides parallel to axis) is inside rect b (arbitrary).
	// Rects must be passed in the [(0,0), (w,0), (w,h), (0,h)] order.
	var isRectInside = function(a, b) {
		var i=0,j=0;
	    for (i = 0; i < 4; ++i)
	        if (b[i].x > a[0].x && b[i].x < a[2].x && b[i].y > a[0].y && b[i].y < a[2].y)
	            return false;
	    for (i = 0; i < 4; ++i)
	    for (j = 0; j < 4; ++j)
	        if (segmentsIntersect(a[i], a[(i+1)&3], b[j], b[(j+1)&3]))
	            return false;
	    return true;
	}


	var isGoodMotion = function(M, w, h, dx, dy) {
	    var pt = [ {'x':0,'y':0}, {'x':w,'y':0}, {'x':w,'y':h}, {'x':0,'y':h} ];
	    var Mpt = [ {'x':0,'y':0}, {'x':0,'y':0}, {'x':0,'y':0}, {'x':0,'y':0} ];
	    var z=0.0, i=0;

	    for (; i < 4; ++i) {
	        Mpt[i].x = M[0]*pt[i].x + M[1]*pt[i].y + M[2];
	        Mpt[i].y = M[3]*pt[i].x + M[4]*pt[i].y + M[5];
	        z = M[6]*pt[i].x + M[7]*pt[i].y + M[8];
	        Mpt[i].x /= z;
	        Mpt[i].y /= z;
	    }

	    pt[0].x = dx, pt[0].y = dy;
	    pt[1].x = w - dx, pt[1].y = dy;
	    pt[2].x = w - dx, pt[2].y = h - dy;
	    pt[3].x = dx, pt[3].y = h - dy;

	    return isRectInside(pt, Mpt);
	}


	var relaxMotion = function(M, t, res) {
	    res[0] = M[0]*(1.0-t) + t;
	    res[1] = M[1]*(1.0-t);
	    res[2] = M[2]*(1.0-t);
	    res[3] = M[3]*(1.0-t);
	    res[4] = M[4]*(1.0-t) + t;
	    res[5] = M[5]*(1.0-t);
	    res[6] = M[6]*(1.0-t);
	    res[7] = M[7]*(1.0-t);
	    res[8] = M[8]*(1.0-t) + t;
	}


	var ensureInclusionConstraint = function(M, size, trimRatio) {
	    var w = size.width;
	    var h = size.height;
	    var dx = Math.floor(w * trimRatio);
	    var dy = Math.floor(h * trimRatio);
	    var md = M.data;
	    var srcM = [ md[0], md[1], md[2],
	                 md[3], md[4], md[5],
	                 md[7], md[7], md[8]];

	    var curM = [];
	    var t = 0.0;
	    relaxMotion(srcM, t, curM);
	    if (isGoodMotion(curM, w, h, dx, dy))
	        return M;

	    var l = 0.0, r = 1.0;
	    while (r - l > 1e-3)
	    {
	        t = (l + r) * 0.5;
	        relaxMotion(srcM, t, curM);
	        if (isGoodMotion(curM, w, h, dx, dy))
	            r = t;
	        else
	            l = t;
	    }

	    l = 1.0 - r;
	    md[0] = md[0]*l+r, md[1] = md[1]*l,   md[2] = md[2]*l;
	    md[3] = md[3]*l,   md[4] = md[4]*l+r, md[5] = md[5]*l;
	    md[6] = md[6]*l,   md[7] = md[7]*l,   md[8] = md[8]*l+r;
	}

	var estimateOptimalTrimRatio = function(M, size) {
	    var w = size.width;
	    var h = size.height;
	    var md=M.data;

	    var pt = [ {'x':0,'y':0}, {'x':w,'y':0}, {'x':w,'y':h}, {'x':0,'y':h} ];
	    var Mpt = [ {'x':0,'y':0}, {'x':0,'y':0}, {'x':0,'y':0}, {'x':0,'y':0} ];
	    var z=0.0;
	    var i=0;

	    for (; i < 4; ++i) {
	        Mpt[i].x = md[0]*pt[i].x + md[1]*pt[i].y + md[2];
	        Mpt[i].y = md[3]*pt[i].x + md[4]*pt[i].y + md[5];
	        z = md[6]*pt[i].x + md[7]*pt[i].y + md[8];
	        Mpt[i].x /= z;
	        Mpt[i].y /= z;
	    }

	    var l = 0, r = 0.5;
	    while (r - l > 1e-3)
	    {
	        var t = (l + r) * 0.5;
	        var dx = Math.floor(w * t);
	        var dy = Math.floor(h * t);
	        pt[0].x = dx, pt[0].y = dy;
	        pt[1].x = w - dx, pt[1].y = dy;
	        pt[2].x = w - dx, pt[2].y = h - dy;
	        pt[3].x = dx, pt[3].y = h - dy;
	        if (isRectInside(pt, Mpt))
	            r = t;
	        else
	            l = t;
	    }

	    return r;
	}
    //

    var onepass_stabilizer = (function () {

	    function onepass_stabilizer() {
	        this.radius = 15;
		    this.trimRatio = 0.0;
		    this.doImage = true;
		    this.doCorrectionForInclusion = false;
		    this.frameSize = new videostab.size_t(0,0);
		    this.frameBufferSize = -1;
		    this.curPos = -1;
		    this.curStabilizedPos = -1;
		    this.preProcessedFrame = null;
		    this.frames = [];
		    this.rgba = []; // original input images
		    this.motions = []; // motions[i] is the motion from i-th to i+1-th frame
		    this.stabilizedFrames = [];
		    this.stabilizationMotions = [];

		    this.motion_estimator = null;
		    this.motion_filter = new videostab.gauss_motion_filter(this.radius);

		    this.reset();
	    };

	    onepass_stabilizer.prototype.reset = function() {
	    	this.frameSize.width = this.frameSize.height = 0;
		    this.curPos = -1;
		    this.curStabilizedPos = -1;
		    this.frameBufferSize = -1;
		    this.preProcessedFrame = null;
		    this.frames.length = 0;
		    this.rgba.length = 0;
		    this.motions.length = 0;
		    this.stabilizedFrames.length = 0;
		    this.stabilizationMotions.length = 0;
	    };

	    var eye3x3 = new jsfeat.matrix_t(3, 3, jsfeat.F32C1_t);
		jsfeat.matmath.identity_3x3(eye3x3, 1.0);

	    onepass_stabilizer.prototype.setup = function(rgbaImageData) {
		    var w = this.frameSize.width = rgbaImageData.width;
		    var h = this.frameSize.height = rgbaImageData.height;

		    var cache_size = 2*this.radius + 1;
		    this.frames.length = cache_size;
		    this.rgba.length = cache_size;
		    this.stabilizedFrames.length = cache_size;
		    this.motions.length = cache_size;
		    this.stabilizationMotions.length = cache_size;

		    this.frameBufferSize = cache_size;

		    var i = 0, j=0;

		    // create objects
		    for(; i < cache_size; ++i) {
		    	this.motions[i] = new jsfeat.matrix_t(3, 3, jsfeat.F32C1_t);
		    	this.frames[i] = new jsfeat.matrix_t(w, h, jsfeat.U8C1_t);

		    	this.stabilizationMotions[i] = new jsfeat.matrix_t(3, 3, jsfeat.F32C1_t);
		    	this.stabilizedFrames[i] = new jsfeat.matrix_t(w, h, jsfeat.U8C1_t);
		    }

		    var gray = new jsfeat.matrix_t(w, h, jsfeat.U8C1_t);
		    jsfeat.imgproc.grayscale(rgbaImageData.data, w, h, gray);

		    for (i = -this.radius; i < 0; ++i) {
		    	j = videostab.get_ring_ind(i, cache_size);
		    	eye3x3.copy_to( this.motions[j] );
		        gray.copy_to( this.frames[j] );
		        this.rgba[i] = rgbaImageData;
		    }

		    gray.copy_to( this.frames[0] );
		    this.rgba[0] = rgbaImageData;

		    if(this.doImage) {
		    	this.preProcessedFrame = new jsfeat.matrix_t(w, h, jsfeat.U8C1_t);
		    }
		};

		onepass_stabilizer.prototype.estimate_motion = function()
		{
			var id0 = videostab.get_ring_ind(this.curPos - 1, this.frameBufferSize)|0;
			var id1 = videostab.get_ring_ind(this.curPos, this.frameBufferSize)|0;
		    return this.motion_estimator.estimate( this.frames[id0], this.frames[id1] );
		};


		onepass_stabilizer.prototype.estimate_stabilization_motion = function()
		{
		    return this.motion_filter.stabilize(this.curStabilizedPos, this.motions, 0, this.curPos);
		};


		onepass_stabilizer.prototype.postprocess_frame = function(frame)
		{
		    return frame;
		};

	    onepass_stabilizer.prototype.next_stabilized_frame = function(frame_rgba) {

		    // check if we've processed all frames already
		    if (this.curStabilizedPos == this.curPos && this.curStabilizedPos != -1)
		    {
		        console.log("no more frames");
		        return null;
		    }

		    var processed = true;
		    do processed = this.do_one_iteration(frame_rgba);
		    while (processed && this.curStabilizedPos == -1);

		    // check if the frame source is empty
		    if (this.curStabilizedPos == -1) {
		    	console.log("frame source is empty");
		        return null;
		    }

		    if(!this.doImage) {
		    	return videostab.get_at(this.curStabilizedPos, this.frameBufferSize, this.stabilizationMotions);
		    }

		    return this.postprocess_frame( videostab.get_at(this.curStabilizedPos, this.frameBufferSize, this.stabilizedFrames) );
		};

		onepass_stabilizer.prototype.do_one_iteration = function(frame_rgba) {

		    var m33;

		    if (frame_rgba) {
		        this.curPos++;

		        if (this.curPos > 0) {

		            var gray_frame = videostab.get_at(this.curPos, this.frameBufferSize, this.frames);
		            jsfeat.imgproc.grayscale(frame_rgba.data, gray_frame.cols, gray_frame.rows, gray_frame);

		            this.rgba[videostab.get_ring_ind(this.curPos, this.frameBufferSize)] = frame_rgba;

		            m33 = this.estimate_motion();
		            m33.copy_to( videostab.get_at(this.curPos-1, this.frameBufferSize, this.motions) );

		            if (this.curPos >= this.radius) {
		                this.curStabilizedPos = this.curPos - this.radius;
		                this.stabilize_frame();
		            }
		        }
		        else
		            this.setup(frame_rgba);

		        return true;
		    }

		    if (this.curStabilizedPos < this.curPos) {
		        this.curStabilizedPos++;

		        videostab.get_at(this.curPos, this.frameBufferSize, this.frames).copy_to( videostab.get_at(this.curStabilizedPos + this.radius, this.frameBufferSize, this.frames) );

		        eye3x3.copy_to( videostab.get_at(this.curStabilizedPos+this.radius-1, this.frameBufferSize, this.motions) );

		        this.rgba[videostab.get_ring_ind(this.curStabilizedPos + this.radius, this.frameBufferSize)] = this.rgba[ videostab.get_ring_ind(this.curPos, this.frameBufferSize) ];

		        this.stabilize_frame();

		        return true;
		    }

		    return false;
		};

		var im33 = new jsfeat.matrix_t(3,3, jsfeat.F32C1_t);

		onepass_stabilizer.prototype.stabilize_frame = function() {

		    var stabilizationMotion = this.estimate_stabilization_motion();

		    if (this.doCorrectionForInclusion) {
		    	ensureInclusionConstraint(stabilizationMotion, this.frameSize, this.trimRatio);
		    }

		    stabilizationMotion.copy_to( videostab.get_at(this.curStabilizedPos, this.frameBufferSize, this.stabilizationMotions) );

		    // apply stabilization transformation

		    if(this.doImage) {

		    	videostab.get_at(this.curStabilizedPos, this.frameBufferSize, this.frames).copy_to( this.preProcessedFrame );

			    jsfeat.matmath.invert_3x3(stabilizationMotion, im33);

			    if(this.motion_estimator.motionModel == videostab.MM_AFFINE) {
				    jsfeat.imgproc.warp_affine(this.preProcessedFrame, 
				    						   videostab.get_at(this.curStabilizedPos, this.frameBufferSize, this.stabilizedFrames), 
				    						   im33, 0);
				} else {
					jsfeat.imgproc.warp_perspective(this.preProcessedFrame, 
				    						   videostab.get_at(this.curStabilizedPos, this.frameBufferSize, this.stabilizedFrames), 
				    						   im33, 0);
				}
			} else {
				this.preProcessedFrame = this.rgba[ videostab.get_ring_ind(this.curStabilizedPos, this.frameBufferSize) ];
			}
		};

	    return onepass_stabilizer;
	})();

	global.onepass_stabilizer = onepass_stabilizer;

})(videostab);
