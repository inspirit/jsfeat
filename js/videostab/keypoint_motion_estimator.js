/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 */

(function(global) {
    "use strict";
    //

    var keypoint_motion_estimator = (function () {

		function keypoint_motion_estimator(motionModel, size) {
			if (typeof motionModel === "undefined") { motionModel=videostab.MM_AFFINE; }

			if(motionModel == videostab.MM_AFFINE) {
				this.mm_kernel = new jsfeat.motion_model.affine2d();
			} else if(motionModel == videostab.MM_HOMOGRAPHY) {
				this.mm_kernel = new jsfeat.motion_model.homography2d();
			}

			this.motionModel = motionModel;

			// ransac params
			this.ransac_param = new jsfeat.ransac_params_t(4, 2, 0.5, 0.99);
			if(motionModel == videostab.MM_AFFINE) {
				this.ransac_param.size = 3;
			}

			// transform matrix
			this.mm3x3 = new jsfeat.matrix_t(3,3,jsfeat.F32C1_t);

			// pyr flow LK
			this.curr_img_pyr = new jsfeat.pyramid_t(3);
	        this.prev_img_pyr = new jsfeat.pyramid_t(3);
	        this.curr_img_pyr.allocate(size.width, size.height, jsfeat.U8C1_t);
	        this.prev_img_pyr.allocate(size.width, size.height, jsfeat.U8C1_t);

	        this.max_points = 100;
	        this.point_count = 0;
	        this.point_status = new Uint8Array(this.max_points);
	        this.prev_xy = new Float32Array(this.max_points*2);
	        this.curr_xy = new Float32Array(this.max_points*2);

	        this.match_mask = new jsfeat.matrix_t(this.max_points,1,jsfeat.U8C1_t);

	        this.corners0 = [];
	        this.corners1 = [];
	        var i = (size.width*size.height) >> 1;
	        while(--i >= 0) {
	            this.corners0[i] = new jsfeat.keypoint_t(0,0,0,0);
	            this.corners1[i] = new jsfeat.keypoint_t(0,0,0,0);
	        }
	    };

	    keypoint_motion_estimator.prototype.estimate = function(frame0, frame1) {

	    	// update pyramids
		    //frame0.copy_to(this.prev_img_pyr.data[0]);
		    //frame1.copy_to(this.curr_img_pyr.data[0]);
		    jsfeat.imgproc.box_blur_gray(frame0, this.prev_img_pyr.data[0], 2, 0);
		    jsfeat.imgproc.box_blur_gray(frame1, this.curr_img_pyr.data[0], 2, 0);

	    	// find keypoints
		    var count = jsfeat.yape06.detect(this.prev_img_pyr.data[0], this.corners0);
		    if(count > this.max_points) {
		    	jsfeat.math.qsort(this.corners0, 0, count-1, function(a,b){return (b.score<a.score);});
		    	count = this.max_points;
		    }

		    //console.log("keypoints: " + count + " [" + this.corners0[0].score + ", " + this.corners0[10].score + "]");

		    // extract coords from keypoints
		    var i=0;
		    for (; i < count; ++i) {
		        this.prev_xy[i<<1] = this.corners0[i].x;
		        this.prev_xy[(i<<1)+1] = this.corners0[i].y;
		    }

		    this.curr_img_pyr.build(this.curr_img_pyr.data[0], true);
		    this.prev_img_pyr.build(this.prev_img_pyr.data[0], true);

		    // find correspondences
		    var lk_win = 20;
			var lk_iters = 30;
		    var lk_epsilon = 0.01;
	        var lk_min_eigen = 0.001;
		    jsfeat.optical_flow_lk.track(this.prev_img_pyr, this.curr_img_pyr, this.prev_xy, this.curr_xy, count, 
		    							lk_win|0, lk_iters|0, this.point_status, lk_epsilon, lk_min_eigen);

		    // leave good correspondences only

		    var good_cnt = 0;
		    for (i = 0; i < count; ++i) {
		        if (this.point_status[i]) {
		        	this.corners0[good_cnt].x = this.corners0[i].x;
		        	this.corners0[good_cnt].y = this.corners0[i].y;

		        	this.corners1[good_cnt].x = this.curr_xy[i<<1];
		        	this.corners1[good_cnt].y = this.curr_xy[(i<<1)+1];

		            good_cnt++;
		        }
		    }

		    //console.log("lk tracked: " + good_cnt);

		    // perform outlier rejection
		    // todo or not todo

		    // estimate motion
		    var ok = false;
		    if(this.motionModel == videostab.MM_AFFINE) {
		    	ok = jsfeat.motion_estimator.ransac(this.ransac_param, this.mm_kernel, 
		    										this.corners0, this.corners1, good_cnt, this.mm3x3, this.match_mask, 1000);
		    } else {
		    	ok = jsfeat.motion_estimator.lmeds(this.ransac_param, this.mm_kernel, 
		    										this.corners0, this.corners1, good_cnt, this.mm3x3, this.match_mask, 1000);
		    }

	    	// extract good matches and re-etimate
	    	if(ok && good_cnt > 15) {
		    	count = good_cnt;
		    	good_cnt = 0;
		    	for(i=0; i < count; ++i) {
		    		if(this.match_mask.data[i]) {
		    			this.corners0[good_cnt].x = this.corners0[i].x;
			        	this.corners0[good_cnt].y = this.corners0[i].y;
			        	this.corners1[good_cnt].x = this.corners1[i].x;
			        	this.corners1[good_cnt].y = this.corners1[i].y;
			        	good_cnt++;
		    		}
		    	}
		    	//console.log("mask: " + good_cnt);
		    	this.mm_kernel.run(this.corners0, this.corners1, this.mm3x3, good_cnt);
		    } else {
		    	jsfeat.matmath.identity_3x3(this.mm3x3, 1.0);
		    }

	    	return this.mm3x3;
	    };

	    return keypoint_motion_estimator;
	})();

	global.keypoint_motion_estimator = keypoint_motion_estimator;

})(videostab);
