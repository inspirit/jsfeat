/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 */

(function(global) {
    "use strict";
    //

    var motion_model = (function() {

    	var sqr = function(x) {
    		return x*x;
    	}

    	// does isotropic normalization
    	var iso_normalize_points = function(from, to, T0, T1, count) {
			var i=0;
		    var cx0=0.0, cy0=0.0, d0=0.0, s0=0.0;
		    var cx1=0.0, cy1=0.0, d1=0.0, s1=0.0;
		    var dx=0.0,dy=0.0;

		    for (; i < count; ++i) {
		        cx0 += from[i].x;
		        cy0 += from[i].y;
		        cx1 += to[i].x;
		        cy1 += to[i].y;
		    }

		    cx0 /= count; cy0 /= count;
		    cx1 /= count; cy1 /= count;

		    for (i = 0; i < count; ++i) {
		        dx = from[i].x - cx0;
		        dy = from[i].y - cy0;
		        d0 += Math.sqrt(dx*dx + dy*dy);
		        dx = to[i].x - cx1;
		        dy = to[i].y - cy1;
		        d1 += Math.sqrt(dx*dx + dy*dy);
		    }

		    d0 /= count; d1 /= count;

		    s0 = Math.SQRT2 / d0; s1 = Math.SQRT2 / d1;

		    T0[0] = T0[4] = s0;
		    T0[2] = -cx0*s0;
		    T0[5] = -cy0*s0;
		    T0[1] = T0[3] = T0[6] = T0[7] = 0.0;
		    T0[8] = 1.0;

		    T1[0] = T1[4] = s1;
		    T1[2] = -cx1*s1;
		    T1[5] = -cy1*s1;
		    T1[1] = T1[3] = T1[6] = T1[7] = 0.0;
		    T1[8] = 1.0;
		}

		var have_collinear_points = function(points, count) {
		    var j=0,k=0,i=(count-1)|0;
		    var dx1=0.0,dy1=0.0,dx2=0.0,dy2=0.0;

		    // check that the i-th selected point does not belong
		    // to a line connecting some previously selected points
		    for(; j < i; ++j) {
		        dx1 = points[j].x - points[i].x;
		        dy1 = points[j].y - points[i].y;
		        for(k = 0; k < j; ++k) {
		            dx2 = points[k].x - points[i].x;
		            dy2 = points[k].y - points[i].y;
		            if( Math.abs(dx2*dy1 - dy2*dx1) <= jsfeat.EPSILON*(Math.abs(dx1) + Math.abs(dy1) + Math.abs(dx2) + Math.abs(dy2)))
		                return true;
		        }
		    }
		    return false;
		}

		var T0 = new jsfeat.matrix_t(3, 3, jsfeat.F32_t|jsfeat.C1_t);
    	var T1 = new jsfeat.matrix_t(3, 3, jsfeat.F32_t|jsfeat.C1_t);
    	var AtA = new jsfeat.matrix_t(6, 6, jsfeat.F32_t|jsfeat.C1_t);
    	var AtB = new jsfeat.matrix_t(6, 1, jsfeat.F32_t|jsfeat.C1_t);
    	
    	var affine2d = (function () {

	        function affine2d() {
	        	// empty constructor
	        }

	        affine2d.prototype.run = function(from, to, model, count) {
	        	var i=0,j=0;
	        	var dt=model.type|jsfeat.C1_t;
	        	var md=model.data, t0d=T0.data, t1d=T1.data;
	        	var pt0,pt1,px=0.0,py=0.0;

	            iso_normalize_points(from, to, t0d, t1d, count);

	            var a_buff = jsfeat.cache.get_buffer((2*count*6)<<3);
                var b_buff = jsfeat.cache.get_buffer((2*count)<<3);

                var a_mt = new jsfeat.matrix_t(6, 2*count, dt, a_buff.data);
                var b_mt = new jsfeat.matrix_t(1, 2*count, dt, b_buff.data);
                var ad=a_mt.data, bd=b_mt.data;

			    for (; i < count; ++i) {
			    	pt0 = from[i];
			        pt1 = to[i];

			        px = t0d[0]*pt0.x + t0d[1]*pt0.y + t0d[2];
			        py = t0d[3]*pt0.x + t0d[4]*pt0.y + t0d[5];

			        j = i*2*6;
			        ad[j]=px, ad[j+1]=py, ad[j+2]=1.0, ad[j+3]=0.0, ad[j+4]=0.0, ad[j+5]=0.0;

			        j += 6;
			        ad[j]=0.0, ad[j+1]=0.0, ad[j+2]=0.0, ad[j+3]=px, ad[j+4]=py, ad[j+5]=1.0;

			        bd[i<<1] = t1d[0]*pt1.x + t1d[1]*pt1.y + t1d[2];
			        bd[(i<<1)+1] = t1d[3]*pt1.x + t1d[4]*pt1.y + t1d[5];
			    }

			    jsfeat.matmath.multiply_AtA(AtA, a_mt);
			    jsfeat.matmath.multiply_AtB(AtB, a_mt, b_mt);

			    jsfeat.linalg.lu_solve(AtA, AtB);

			    md[0] = AtB.data[0], md[1]=AtB.data[1], md[2]=AtB.data[2];
			    md[3] = AtB.data[3], md[4]=AtB.data[4], md[5]=AtB.data[5];
			    md[6] = 0.0, md[7] = 0.0, md[8] = 1.0; // fill last row

			    // denormalize
			    jsfeat.matmath.invert_3x3(T1, T1);
			    jsfeat.matmath.multiply_3x3(model, T1, model);
			    jsfeat.matmath.multiply_3x3(model, model, T0);

			    // free buffer
			    jsfeat.cache.put_buffer(a_buff);
			    jsfeat.cache.put_buffer(b_buff);

			    return 1;
	        }

	        affine2d.prototype.error = function(from, to, model, err, count) {
	        	var i=0;
	        	var pt0,pt1;
	        	var m=model.data;

			    for (; i < count; ++i) {
			        pt0 = from[i];
			        pt1 = to[i];

			        err[i] = sqr(pt1.x - m[0]*pt0.x - m[1]*pt0.y - m[2]) +
			                 sqr(pt1.y - m[3]*pt0.x - m[4]*pt0.y - m[5]);
			    }
	        }

	        affine2d.prototype.check_subset = function(from, to, count) {
	            return true; // all good
	        }

	        return affine2d;
	    })();

	    var mLtL = new jsfeat.matrix_t(9, 9, jsfeat.F32_t|jsfeat.C1_t);
	    var Evec = new jsfeat.matrix_t(9, 9, jsfeat.F32_t|jsfeat.C1_t);

	    var homography2d = (function () {

	        function homography2d() {
	        	// empty constructor
	        	//this.T0 = new jsfeat.matrix_t(3, 3, jsfeat.F32_t|jsfeat.C1_t);
	        	//this.T1 = new jsfeat.matrix_t(3, 3, jsfeat.F32_t|jsfeat.C1_t);
	        	//this.mLtL = new jsfeat.matrix_t(9, 9, jsfeat.F32_t|jsfeat.C1_t);
	        	//this.Evec = new jsfeat.matrix_t(9, 9, jsfeat.F32_t|jsfeat.C1_t);
	        }

	        homography2d.prototype.run = function(from, to, model, count) {
	        	var i=0,j=0;
	        	var md=model.data, t0d=T0.data, t1d=T1.data;
	        	var LtL=mLtL.data, evd=Evec.data;
	        	var x=0.0,y=0.0,X=0.0,Y=0.0;

			    // norm
				var smx=0.0, smy=0.0, cmx=0.0, cmy=0.0, sMx=0.0, sMy=0.0, cMx=0.0, cMy=0.0;

				for(; i < count; ++i) {
				    cmx += to[i].x;
				    cmy += to[i].y;
				    cMx += from[i].x;
				    cMy += from[i].y;
				}

			    cmx /= count; cmy /= count;
			    cMx /= count; cMy /= count;

			    for(i = 0; i < count; ++i)
			    {
				    smx += Math.abs(to[i].x - cmx);
				    smy += Math.abs(to[i].y - cmy);
				    sMx += Math.abs(from[i].x - cMx);
				    sMy += Math.abs(from[i].y - cMy);
				}

			    if( Math.abs(smx) < jsfeat.EPSILON 
			    	|| Math.abs(smy) < jsfeat.EPSILON 
			    	|| Math.abs(sMx) < jsfeat.EPSILON 
			    	|| Math.abs(sMy) < jsfeat.EPSILON ) return 0;

			    smx = count/smx; smy = count/smy;
			    sMx = count/sMx; sMy = count/sMy;

			    t0d[0] = sMx; 	t0d[1] = 0; 	t0d[2] = -cMx*sMx; 
			    t0d[3] = 0; 	t0d[4] = sMy; 	t0d[5] = -cMy*sMy; 
			    t0d[6] = 0; 	t0d[7] = 0; 	t0d[8] = 1;

				t1d[0] = 1.0/smx; 	t1d[1] = 0; 		t1d[2] = cmx;
				t1d[3] = 0; 		t1d[4] = 1.0/smy; 	t1d[5] = cmy;
				t1d[6] = 0; 		t1d[7] = 0; 		t1d[8] = 1;
				//

				// construct system
				i = 81;
				while(--i >= 0) {
					LtL[i] = 0.0;
				}
				for(i = 0; i < count; ++i) {
					x = (to[i].x - cmx) * smx;
					y = (to[i].y - cmy) * smy;
					X = (from[i].x - cMx) * sMx;
					Y = (from[i].y - cMy) * sMy;

					LtL[0] += X*X;
					LtL[1] += X*Y;
					LtL[2] += X;

					LtL[6] += X*-x*X;
					LtL[7] += X*-x*Y;
					LtL[8] += X*-x;
					LtL[10] += Y*Y;
					LtL[11] += Y;

					LtL[15] += Y*-x*X;
					LtL[16] += Y*-x*Y;
					LtL[17] += Y*-x;
					LtL[20] += 1.0;

					LtL[24] += -x*X;
					LtL[25] += -x*Y;
					LtL[26] += -x;
					LtL[30] += X*X;
					LtL[31] += X*Y;
					LtL[32] += X;
					LtL[33] += X*-y*X;
					LtL[34] += X*-y*Y;
					LtL[35] += X*-y;
					LtL[40] += Y*Y;
					LtL[41] += Y;
					LtL[42] += Y*-y*X;
					LtL[43] += Y*-y*Y;
					LtL[44] += Y*-y;
					LtL[50] += 1.0;
					LtL[51] += -y*X;
					LtL[52] += -y*Y;
					LtL[53] += -y;
					LtL[60] += -x*X*-x*X + -y*X*-y*X;
					LtL[61] += -x*X*-x*Y + -y*X*-y*Y;
					LtL[62] += -x*X*-x + -y*X*-y;
					LtL[70] += -x*Y*-x*Y + -y*Y*-y*Y;
					LtL[71] += -x*Y*-x + -y*Y*-y;
					LtL[80] += -x*-x + -y*-y;
				}
				//

				// symmetry
			    for(i = 0; i < 9; ++i) {
			        for(j = 0; j < i; ++j)
			            LtL[i*9+j] = LtL[j*9+i];
			    }

				jsfeat.linalg.eigenVV(mLtL, Evec);

				md[0]=evd[72], md[1]=evd[73], md[2]=evd[74];
			    md[3]=evd[75], md[4]=evd[76], md[5]=evd[77];
			    md[6]=evd[78], md[7]=evd[79], md[8]=evd[80];

				// denormalize
			    jsfeat.matmath.multiply_3x3(model, T1, model);
			    jsfeat.matmath.multiply_3x3(model, model, T0);

			    // set bottom right to 1.0
			    x = 1.0/md[8];
			    md[0] *= x; md[1] *= x; md[2] *= x;
			    md[3] *= x; md[4] *= x; md[5] *= x;
			    md[6] *= x; md[7] *= x; md[8] = 1.0;

			    return 1;
	        }

	        homography2d.prototype.error = function(from, to, model, err, count) {
	        	var i=0;
	        	var pt0,pt1,ww=0.0,dx=0.0,dy=0.0;
	        	var m=model.data;

			    for (; i < count; ++i) {
			        pt0 = from[i];
			        pt1 = to[i];

			        ww = 1.0/(m[6]*pt0.x + m[7]*pt0.y + 1.0);
			        dx = (m[0]*pt0.x + m[1]*pt0.y + m[2])*ww - pt1.x;
			        dy = (m[3]*pt0.x + m[4]*pt0.y + m[5])*ww - pt1.y;
			        err[i] = (dx*dx + dy*dy);
			    }
	        }

	        homography2d.prototype.check_subset = function(from, to, count) {
	        	// seems to reject good subsets actually
	        	//if( have_collinear_points(from, count) || have_collinear_points(to, count) ) {
        			//return false;
        		//}
        		if( count == 4 ) {
			        var negative = 0;

			        var fp0=from[0],fp1=from[1],fp2=from[2],fp3=from[3];
			        var tp0=to[0],tp1=to[1],tp2=to[2],tp3=to[3];

			        // set1
			        var A11=fp0.x, A12=fp0.y, A13=1.0;
			        var A21=fp1.x, A22=fp1.y, A23=1.0;
			        var A31=fp2.x, A32=fp2.y, A33=1.0;

			        var B11=tp0.x, B12=tp0.y, B13=1.0;
			        var B21=tp1.x, B22=tp1.y, B23=1.0;
			        var B31=tp2.x, B32=tp2.y, B33=1.0;

			        var detA = jsfeat.matmath.determinant_3x3(A11,A12,A13, A21,A22,A23, A31,A32,A33);
					var detB = jsfeat.matmath.determinant_3x3(B11,B12,B13, B21,B22,B23, B31,B32,B33);

					if(detA*detB < 0) negative++;

					// set2
					A11=fp1.x, A12=fp1.y;
			        A21=fp2.x, A22=fp2.y;
			        A31=fp3.x, A32=fp3.y;

			        B11=tp1.x, B12=tp1.y;
			        B21=tp2.x, B22=tp2.y;
			        B31=tp3.x, B32=tp3.y;

			        detA = jsfeat.matmath.determinant_3x3(A11,A12,A13, A21,A22,A23, A31,A32,A33);
					detB = jsfeat.matmath.determinant_3x3(B11,B12,B13, B21,B22,B23, B31,B32,B33);

					if(detA*detB < 0) negative++;

					// set3
					A11=fp0.x, A12=fp0.y;
			        A21=fp2.x, A22=fp2.y;
			        A31=fp3.x, A32=fp3.y;

			        B11=tp0.x, B12=tp0.y;
			        B21=tp2.x, B22=tp2.y;
			        B31=tp3.x, B32=tp3.y;

			        detA = jsfeat.matmath.determinant_3x3(A11,A12,A13, A21,A22,A23, A31,A32,A33);
					detB = jsfeat.matmath.determinant_3x3(B11,B12,B13, B21,B22,B23, B31,B32,B33);

					if(detA*detB < 0) negative++;

					// set4
					A11=fp0.x, A12=fp0.y;
			        A21=fp1.x, A22=fp1.y;
			        A31=fp3.x, A32=fp3.y;

			        B11=tp0.x, B12=tp0.y;
			        B21=tp1.x, B22=tp1.y;
			        B31=tp3.x, B32=tp3.y;

			        detA = jsfeat.matmath.determinant_3x3(A11,A12,A13, A21,A22,A23, A31,A32,A33);
					detB = jsfeat.matmath.determinant_3x3(B11,B12,B13, B21,B22,B23, B31,B32,B33);

					if(detA*detB < 0) negative++;

			        if(negative != 0 && negative != 4) {
			        	return false;
			        }
			    }
	            return true; // all good
	        }

	        return homography2d;
	    })();

	    return {

    		affine2d:affine2d,
    		homography2d:homography2d

    	};

    })();

    var ransac_params_t = (function () {
        function ransac_params_t(size, thresh, eps, prob) {
            if (typeof size === "undefined") { size=0; }
            if (typeof thresh === "undefined") { thresh=0.5; }
            if (typeof eps === "undefined") { eps=0.5; }
            if (typeof prob === "undefined") { prob=0.99; }

            this.size = size;
            this.thresh = thresh;
            this.eps = eps;
            this.prob = prob;
        };
        ransac_params_t.prototype.update_iters = function(_eps, max_iters) {
	        var num = Math.log(1 - this.prob);
	        var denom = Math.log(1 - Math.pow(1 - _eps, this.size));
	        return (denom >= 0 || -num >= max_iters*(-denom) ? max_iters : Math.round(num/denom))|0;
        };
        return ransac_params_t;
    })();

    var motion_estimator = (function() {

    	var get_subset = function(kernel, from, to, need_cnt, max_cnt, from_sub, to_sub) {
    		var max_try = 1000;
    		var indices = [];
		    var i=0, j=0, ssiter=0, idx_i=0, ok=false;
		    for(; ssiter < max_try; ++ssiter)  {
		        i = 0;
		        for (; i < need_cnt && ssiter < max_try;) {
		            ok = false;
		            idx_i = 0;
		            while (!ok) {
		                ok = true;
		                idx_i = indices[i] = Math.floor(Math.random() * max_cnt)|0;
		                for (j = 0; j < i; ++j) {
		                    if (idx_i == indices[j])
		                    { ok = false; break; }
		                }
		            }
		            from_sub[i] = from[idx_i];
		            to_sub[i] = to[idx_i];
		            if( !kernel.check_subset( from_sub, to_sub, i+1 ) ) {
		                ssiter++;
		                continue;
		            }
		            ++i;
		        }
		        break;
		    }

		    return (i == need_cnt && ssiter < max_try);
    	}

    	var find_inliers = function(kernel, model, from, to, count, thresh, err, mask) {
    		var numinliers = 0, i=0, f=0;
    		var t = thresh*thresh;

    		kernel.error(from, to, model, err, count);

		    for(; i < count; ++i) {
		        f = err[i] <= t;
		        mask[i] = f;
		        numinliers += f;
		    }
		    return numinliers;
    	}

    	return {

    		ransac: function(params, kernel, from, to, count, model, mask, max_iters) {
    			if (typeof max_iters === "undefined") { max_iters=1000; }

    			if(count < params.size) return false;

    			var model_points = params.size;
			    var niters = max_iters, iter=0;
			    var result = false;

			    var subset0 = [];
			    var subset1 = [];
			    var found = false;

			    var mc=model.cols,mr=model.rows;
                var dt = model.type | jsfeat.C1_t;

			    var m_buff = jsfeat.cache.get_buffer((mc*mr)<<3);
			    var ms_buff = jsfeat.cache.get_buffer(count);
			    var err_buff = jsfeat.cache.get_buffer(count<<2);
			    var M = new jsfeat.matrix_t(mc, mr, dt, m_buff.data);
			    var curr_mask = new jsfeat.matrix_t(count, 1, jsfeat.U8C1_t, ms_buff.data);

			    var inliers_max = -1, numinliers=0;
			    var nmodels = 0;

			    var err = err_buff.f32;

			    // special case
			    if(count == model_points) {
			        if(kernel.run(from, to, M, count) <= 0) {
			        	jsfeat.cache.put_buffer(m_buff);
			        	jsfeat.cache.put_buffer(ms_buff);
			        	jsfeat.cache.put_buffer(err_buff);
			        	return false;
			        }

			        M.copy_to(model);
			        if(mask) {
			        	while(--count >= 0) {
			        		mask.data[count] = 1;
			        	}
			        }
			        jsfeat.cache.put_buffer(m_buff);
			        jsfeat.cache.put_buffer(ms_buff);
			        jsfeat.cache.put_buffer(err_buff);
			        return true;
			    }

			    for (; iter < niters; ++iter) {
			        // generate subset
			        found = get_subset(kernel, from, to, model_points, count, subset0, subset1);
			        if(!found) {
			            if(iter == 0) {
			            	jsfeat.cache.put_buffer(m_buff);
			            	jsfeat.cache.put_buffer(ms_buff);
			            	jsfeat.cache.put_buffer(err_buff);
			                return false;
			            }
			            break;
			        }

			        nmodels = kernel.run( subset0, subset1, M, model_points );
			        if(nmodels <= 0)
			            continue;

			        // TODO handle multimodel output

			        numinliers = find_inliers(kernel, M, from, to, count, params.thresh, err, curr_mask.data);

			        if( numinliers > Math.max(inliers_max, model_points-1) ) {
			            M.copy_to(model);
			            inliers_max = numinliers;
			            if(mask) curr_mask.copy_to(mask);
			            niters = params.update_iters((count - numinliers)/count, niters);
			            result = true;
			        }
			    }

			    jsfeat.cache.put_buffer(m_buff);
			    jsfeat.cache.put_buffer(ms_buff);
			    jsfeat.cache.put_buffer(err_buff);

			    return result;
    		},

    		lmeds: function(params, kernel, from, to, count, model, mask, max_iters) {
    			if (typeof max_iters === "undefined") { max_iters=1000; }

    			if(count < params.size) return false;

    			var model_points = params.size;
			    var niters = max_iters, iter=0;
			    var result = false;

			    var subset0 = [];
			    var subset1 = [];
			    var found = false;

			    var mc=model.cols,mr=model.rows;
                var dt = model.type | jsfeat.C1_t;

			    var m_buff = jsfeat.cache.get_buffer((mc*mr)<<3);
			    var ms_buff = jsfeat.cache.get_buffer(count);
			    var err_buff = jsfeat.cache.get_buffer(count<<2);
			    var M = new jsfeat.matrix_t(mc, mr, dt, m_buff.data);
			    var curr_mask = new jsfeat.matrix_t(count, 1, jsfeat.U8_t|jsfeat.C1_t, ms_buff.data);

			    var numinliers=0;
			    var nmodels = 0;

			    var err = err_buff.f32;
			    var min_median = 1000000000.0, sigma=0.0, median=0.0;

			    params.eps = 0.45;
			    niters = params.update_iters(params.eps, niters);

			    // special case
			    if(count == model_points) {
			        if(kernel.run(from, to, M, count) <= 0) {
			        	jsfeat.cache.put_buffer(m_buff);
			        	jsfeat.cache.put_buffer(ms_buff);
			        	jsfeat.cache.put_buffer(err_buff);
			        	return false;
			        }

			        M.copy_to(model);
			        if(mask) {
			        	while(--count >= 0) {
			        		mask.data[count] = 1;
			        	}
			        }
			        jsfeat.cache.put_buffer(m_buff);
			        jsfeat.cache.put_buffer(ms_buff);
			        jsfeat.cache.put_buffer(err_buff);
			        return true;
			    }

			    for (; iter < niters; ++iter) {
			        // generate subset
			        found = get_subset(kernel, from, to, model_points, count, subset0, subset1);
			        if(!found) {
			            if(iter == 0) {
			            	jsfeat.cache.put_buffer(m_buff);
			            	jsfeat.cache.put_buffer(ms_buff);
			            	jsfeat.cache.put_buffer(err_buff);
			                return false;
			            }
			            break;
			        }

			        nmodels = kernel.run( subset0, subset1, M, model_points );
			        if(nmodels <= 0)
			            continue;

			        // TODO handle multimodel output

			        kernel.error(from, to, M, err, count);
			        median = jsfeat.math.median(err, 0, count-1);

			        if(median < min_median) {
			            min_median = median;
			            M.copy_to(model);
			            result = true;
			        }
			    }

			    if(result) {
			        sigma = 2.5*1.4826*(1 + 5.0/(count - model_points))*Math.sqrt(min_median);
			        sigma = Math.max(sigma, 0.001);

			        numinliers = find_inliers(kernel, model, from, to, count, sigma, err, curr_mask.data);
			        if(mask) curr_mask.copy_to(mask);
			        
			        result = numinliers >= model_points;
			    }

			    jsfeat.cache.put_buffer(m_buff);
			    jsfeat.cache.put_buffer(ms_buff);
			    jsfeat.cache.put_buffer(err_buff);

			    return result;
    		}

    	};

    })();

    global.ransac_params_t = ransac_params_t;
    global.motion_model = motion_model;
    global.motion_estimator = motion_estimator;

})(jsfeat);
