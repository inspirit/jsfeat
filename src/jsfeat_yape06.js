/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 * Copyright 2007 Computer Vision Lab,
 * Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland.
 * @author Vincent Lepetit (http://cvlab.epfl.ch/~lepetit)
 */

(function(global) {
    "use strict";
    //

    var yape06 = (function() {
        
        var compute_laplacian = function(src, dst, w, h, Dxx, Dyy) {
            var y=0,x=0,yrow=(Dxx*w)|0,row=yrow;

            for(y = Dxx; y < h - Dxx; ++y, yrow+=w, row = yrow) {
                for(x = w - Dxx; x >= Dxx; --x) {
                    dst[row] = -4 * src[row] + src[row+Dxx] + src[row-Dxx] + src[row+Dyy] + src[row-Dyy];
                    ++row;
                }
            }
        }

        var hessian_min_eigen_value = function(src, off, tr, Dxx, Dyy, Dxy, Dyx) {
            var Ixx = -2 * src[off] + src[off + Dxx] + src[off - Dxx];
            var Iyy = -2 * src[off] + src[off + Dyy] + src[off - Dyy];
            var Ixy = src[off + Dxy] + src[off - Dxy] - src[off + Dyx] - src[off - Dyx];
            var sqrt_delta = ( Math.sqrt(((Ixx - Iyy) * (Ixx - Iyy) + 4 * Ixy * Ixy) ) )|0;

            return Math.min(Math.abs(tr - sqrt_delta), Math.abs(-(tr + sqrt_delta)));
        }

        return {

            laplacian_threshold: 30,
            min_eigen_value_threshold: 25,

            detect: function(src, points) {
                var x=0,y=0;
                var w=src.cols, h=src.rows, srd_d=src.data;
                var Dxx = 5, Dyy = (5 * w)|0;
                var Dxy = (3 + 3 * w)|0, Dyx = (3 - 3 * w)|0;
                var lap_buf = jsfeat.cache.get_buffer((w*h)<<2);
                var laplacian = lap_buf.i32;
                var lv=0, row=0,rowx=0,min_eigen_value=0,pt;
                var number_of_points = 0;
                var lap_thresh = this.laplacian_threshold;
                var eigen_thresh = this.min_eigen_value_threshold;

                x = w*h;
                while(--x>=0) {laplacian[x]=0;}
                compute_laplacian(srd_d, laplacian, w, h, Dxx, Dyy);

                row = (w+1)|0;
                for(y = 1; y < h-1; ++y, row += w) {
                    for(x = 1, rowx=row; x < w-1; ++x, ++rowx) {

                        lv = laplacian[rowx];
                        if ((lv < -lap_thresh &&
                            lv < laplacian[rowx - 1]      && lv < laplacian[rowx + 1] &&
                            lv < laplacian[rowx - w]     && lv < laplacian[rowx + w] &&
                            lv < laplacian[rowx - w - 1] && lv < laplacian[rowx + w - 1] &&
                            lv < laplacian[rowx - w + 1] && lv < laplacian[rowx + w + 1])
                            ||
                            (lv > lap_thresh &&
                            lv > laplacian[rowx - 1]      && lv > laplacian[rowx + 1] &&
                            lv > laplacian[rowx - w]     && lv > laplacian[rowx + w] &&
                            lv > laplacian[rowx - w - 1] && lv > laplacian[rowx + w - 1] &&
                            lv > laplacian[rowx - w + 1] && lv > laplacian[rowx + w + 1])
                            ) {

                            min_eigen_value = hessian_min_eigen_value(srd_d, rowx, lv, Dxx, Dyy, Dxy, Dyx);
                            if (min_eigen_value > eigen_thresh) {
                                pt = points[number_of_points];
                                pt.x = x, pt.y = y, pt.score = min_eigen_value;
                                ++number_of_points;
                                ++x, ++rowx;
                            }
                        }
                    }
                }

                jsfeat.cache.put_buffer(lap_buf);

                return number_of_points;
            }

        };
    })();

    global.yape06 = yape06;

})(jsfeat);
