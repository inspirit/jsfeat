/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(global) {
    "use strict";
    //

    var imgproc = (function() {

        var RGBA_2_GRAY_LUT = new Int32Array(256*3);

        return {

            init_grayscale_lut: function(blueIdx) {
                var coeffs = new Int32Array(3);
                coeffs[0] = 4899; coeffs[1] = 9617; coeffs[2] = 1868;

                var b = 0, g = 0, r = (1 << (14-1));
                var db = coeffs[blueIdx^2], dg = coeffs[1], dr = coeffs[blueIdx];
                for( var i = 0; i < 256; i++, b += db, g += dg, r += dr ) {
                    RGBA_2_GRAY_LUT[i] = b;
                    RGBA_2_GRAY_LUT[i+256] = g;
                    RGBA_2_GRAY_LUT[i+512] = r;
                }
            },

            grayscale: function(src, dst) {
                var srcLength = src.length | 0;
                var _lut = RGBA_2_GRAY_LUT;
                var j = 0;
                for (var i = 0; i < srcLength; i += 4, ++j) {
                    dst[j] = (_lut[src[i]] + _lut[src[i+1]+256] + _lut[src[i+2]+512]) >> 14;
                    //dst[j] = ((src[i]*77)+(src[i+1]*151)+(src[i+2]*28)) >> 8;
                }
            },

            box_blur: function(src, dst, w, h, hwin) {
                var win = (2*hwin+1)|0;
                var _buf = new Int32Array(w*win+w);
                var sums = (w*win)|0;
                var next_row=0,oldest_row=0;
                var i=0, j=0, s0=0, s1=0;
                var hsum = 0, back=0;

                var input = 0, output = (w*hwin - hwin) | 0;

                for(i = 0; i < w; ++i){ _buf[sums+i] = 0; }

                for (i=0; i<h; ++i) {
                    hsum = 0;
                    back = input;

                    for (j = 0; j < win-1; ++j) hsum += src[input+j];
                    for (; j <= w-2; j+=2) {
                        hsum += src[input+j];
                        _buf[next_row+j] = hsum;
                        _buf[sums+j] += hsum;
                        hsum -= src[back++];
                        //
                        hsum += src[input+j+1];
                        _buf[next_row+j+1] = hsum;
                        _buf[sums+j+1] += hsum;
                        hsum -= src[back++];
                    }
                    for (; j < w; ++j) {
                        hsum += src[input+j];
                        _buf[next_row+j] = hsum;
                        _buf[sums+j] += hsum;
                        hsum -= src[back++];
                    }
                    if (i >= win-1)  {
                        for(j = win-1; j <= w-2; j+=2) {
                            s0 = _buf[sums+j]; s1 = _buf[sums+j+1];
                            dst[output+j] = s0; dst[output+j+1] = s1;
                            _buf[sums+j] = s0 - _buf[oldest_row+j];
                            _buf[sums+j+1] = s1 - _buf[oldest_row+j+1];
                        }
                        for(; j < w; ++j) {
                            s0 = _buf[sums+j];
                            dst[output+j] = s0;
                            _buf[sums+j] = s0 - _buf[oldest_row+j];
                        }
                        output += w;
                        oldest_row += w;
                        oldest_row = (oldest_row < sums) * oldest_row;
                    }
                    input += w;
                    next_row += w;
                    next_row = (next_row < sums) * next_row;
                }
            },

            pyrdown_fast_u8: function(src, dst) {
                var w = src.width | 0;
                var w2 = w >> 1, h2 = src.height >> 1;
                var src_data = src.data_u8;
                var dst_data = dst.data_u8;
                var x=0,y=0,sptr=0,sline=0,dptr=0;

                for(y = 0; y < h2; ++y) {
                    sline = sptr;
                    for(x = 0; x < w2; ++x, ++dptr) {
                        dst_data[dptr] = (src_data[sline] + src_data[sline+1] +
                                            src_data[sline+w] + src_data[sline+w+1] + 2) >> 2;
                        sline += 2;
                    }
                    sptr += w << 1;
                }
            },

            shar_derivatives: function(img, dst, w, h) {
                var dstep = w<<1,x=0,y=0,x1=0,t0=0,t1=0;
                var srow0=0,srow1=0,srow2=0,drow=0;
                var trow0 = new Int32Array(w+2);
                var trow1 = new Int32Array(w+2);

                for(y=0; y < h; ++y) {
                    srow1 = (y*w)|0;
                    srow0 = ((y > 0 ? y-1 : 1)*w)|0;
                    srow2 = ((y < h-1 ? y+1 : h-2)*w)|0;
                    drow = (y*dstep)|0;
                    // do vertical convolution
                    for(x = 0, x1 = 1; x < w; ++x,++x1) {
                        t0 = ( ((img[srow0+x]) + (img[srow2+x]))*3 + (img[srow1+x])*10 );
                        t1 = ( (img[srow2+x]) - (img[srow0+x]) );
                        trow0[x1] = t0;
                        trow1[x1] = t1;
                    }
                    // make border
                    x = (w + 1)|0;
                    trow0[0] = trow0[1]; trow0[x] = trow0[w];
                    trow1[0] = trow1[1]; trow1[x] = trow1[w];
                    // do horizontal convolution, interleave the results and store them
                    for(x = 0; x < w; ++x) {
                        t0 = ( (trow0[x+2] - trow0[x]) );
                        t1 = ( ((trow1[x+2] + trow1[x])*3 + trow1[x+1]*10) );
                        dst[drow++] = t0;
                        dst[drow++] = t1;
                    }
                }
            }
        };
    })();

    global.imgproc = imgproc;

})(jsfeat);
