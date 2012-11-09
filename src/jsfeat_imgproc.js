/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(global) {
    "use strict";
    //

    var imgproc = (function() {

        return {
            // TODO: add support for RGB/BGR order
            grayscale: function(src, dst) {
                var srcLength = src.length|0, srcLength_16 = (srcLength - 16)|0;
                var j = 0;
                for (var i = 0; i <= srcLength_16; i += 16, j += 4) {
                    dst[j] = ((src[i] * 77) + (src[i + 1] * 151) + (src[i + 2] * 28)) >> 8;
                    dst[j + 1] = ((src[i + 4] * 77) + (src[i + 5] * 151) + (src[i + 6] * 28)) >> 8;

                    dst[j + 2] = ((src[i + 8] * 77) + (src[i + 9] * 151) + (src[i + 10] * 28)) >> 8;
                    dst[j + 3] = ((src[i + 12] * 77) + (src[i + 13] * 151) + (src[i + 14] * 28)) >> 8;
                }
                for (; i < srcLength; i += 4, ++j) {
                    dst[j] = ((src[i] * 77) + (src[i + 1] * 151) + (src[i + 2] * 28)) >> 8;
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

                for(; i < w; ++i){ _buf[sums+i] = 0; }

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

            gaussian_blur_u8: function(src, dst, w, h, kernel_size, sigma) {
                if (typeof sigma === "undefined") { sigma = 0.0; }
                if (typeof kernel_size === "undefined") { kernel_size = 0; }
                kernel_size = kernel_size == 0 ? (Math.max(1, (4.0 * sigma + 1.0 - 1e-8)|0) * 2 + 1)|0 : kernel_size;
                var half_kernel = kernel_size >> 1;
                var i=0,j=0,k=0,sum=0,sum1=0,sp=0,dp=0,w2=w<<1;
                var buf = new Uint8Array(kernel_size + Math.max(h, w));
                // we use int based kernel
                var filter = new Int32Array(kernel_size);

                jsfeat.math.get_gaussian_kernel(kernel_size, sigma, filter, jsfeat.U8_t);

                // hor pass
                for (; i < h; ++i) { 
                    sum = src[sp];
                    for (j = 0; j < half_kernel; ++j) {
                        buf[j] = sum;
                    }
                    for (j = 0; j <= w-2; j+=2) {
                        buf[j + half_kernel] = src[sp+j];
                        buf[j + half_kernel+1] = src[sp+j+1];
                    }
                    for (; j < w; ++j) {
                        buf[j + half_kernel] = src[sp+j];
                    }
                    sum = src[sp+w-1];
                    for (j = w; j < half_kernel + w; ++j) {
                        buf[j + half_kernel] = sum;
                    }
                    for (j = 0; j <= w-2; j+=2) {
                        sum = 0, sum1 = 0;
                        for (k = 0; k < kernel_size; ++k) {
                            sum += buf[k + j] * filter[k];
                            sum1 += buf[k + j+1] * filter[k];
                        }
                        dst[dp+j] = sum >> 8;
                        dst[dp+j+1] = sum1 >> 8;
                    }
                    for (; j < w; ++j) {
                        sum = 0;
                        for (k = 0; k < kernel_size; ++k) {
                            sum += buf[k + j] * filter[k];
                        }
                        dst[dp+j] = sum >> 8;
                    }
                    sp += w;
                    dp += w;
                }

                // vert pass
                for (i = 0; i < w; ++i) {
                    sum = dst[i];
                    for (j = 0; j < half_kernel; ++j) {
                        buf[j] = sum;
                    }
                    k = i;
                    for (j = 0; j <= h-2; j+=2, k+=w2) {
                        buf[j+half_kernel] = dst[k];
                        buf[j+half_kernel+1] = dst[k+w];
                    }
                    for (; j < h; ++j, k+=w) {
                        buf[j+half_kernel] = dst[k];
                    }
                    sum = dst[(h-1)*w + i];
                    for (j = h; j < half_kernel + h; ++j) {
                        buf[j + half_kernel] = sum;
                    }
                    dp = i;
                    for (j = 0; j <= h-2; j+=2, dp+=w2) { 
                        sum = 0, sum1=0;
                        for (k = 0; k < kernel_size; ++k) {
                            sum += buf[k + j] * filter[k];
                            sum1 += buf[k + j+1] * filter[k];
                        }
                        dst[dp] = sum >> 8;
                        dst[dp+w] = sum1 >> 8;
                    }
                    for (; j < h; ++j, dp+=w) {
                        sum = 0;
                        for (k = 0; k < kernel_size; ++k) {
                            sum += buf[k + j] * filter[k];
                        }
                        dst[dp] = sum >> 8;
                    }
                }
            },

            pyrdown_fast_u8: function(src, w, h, dst) {
                var w2 = w >> 1, h2 = h >> 1;
                var x=0,y=0,sptr=0,sline=0,dptr=0;

                for(y = 0; y < h2; ++y) {
                    sline = sptr;
                    for(x = 0; x <= w2-2; x+=2, dptr+=2, sline += 4) {
                        dst[dptr] = (src[sline] + src[sline+1] +
                                            src[sline+w] + src[sline+w+1] + 2) >> 2;
                        dst[dptr+1] = (src[sline+2] + src[sline+3] +
                                            src[sline+w+2] + src[sline+w+3] + 2) >> 2;
                    }
                    for(; x < w2; ++x, ++dptr, sline += 2) {
                        dst[dptr] = (src[sline] + src[sline+1] +
                                            src[sline+w] + src[sline+w+1] + 2) >> 2;
                    }
                    sptr += w << 1;
                }
            },

            shar_derivatives: function(img, dst, w, h) {
                var dstep = w<<1,x=0,y=0,x1=0;
                var srow0=0,srow1=0,srow2=0,drow=0;
                var trow0 = new Int32Array(w+2);
                var trow1 = new Int32Array(w+2);

                for(; y < h; ++y, srow1+=w) {
                    srow0 = ((y > 0 ? y-1 : 1)*w)|0;
                    srow2 = ((y < h-1 ? y+1 : h-2)*w)|0;
                    drow = (y*dstep)|0;
                    // do vertical convolution
                    for(x = 0, x1 = 1; x <= w-2; x+=2, x1+=2) {
                        trow0[x1] = ( ((img[srow0+x]) + (img[srow2+x]))*3 + (img[srow1+x])*10 );
                        trow1[x1] = ( (img[srow2+x]) - (img[srow0+x]) );
                        //
                        trow0[x1+1] = ( ((img[srow0+x+1]) + (img[srow2+x+1]))*3 + (img[srow1+x+1])*10 );
                        trow1[x1+1] = ( (img[srow2+x+1]) - (img[srow0+x+1]) );
                    }
                    for(; x < w; ++x, ++x1) {
                        trow0[x1] = ( ((img[srow0+x]) + (img[srow2+x]))*3 + (img[srow1+x])*10 );
                        trow1[x1] = ( (img[srow2+x]) - (img[srow0+x]) );
                    }
                    // make border
                    x = (w + 1)|0;
                    trow0[0] = trow0[1]; trow0[x] = trow0[w];
                    trow1[0] = trow1[1]; trow1[x] = trow1[w];
                    // do horizontal convolution, interleave the results and store them
                    for(x = 0; x <= w-2; x+=2) {
                        dst[drow++] = ( (trow0[x+2] - trow0[x]) );
                        dst[drow++] = ( ((trow1[x+2] + trow1[x])*3 + trow1[x+1]*10) );
                        dst[drow++] = ( (trow0[x+3] - trow0[x+1]) );
                        dst[drow++] = ( ((trow1[x+3] + trow1[x+1])*3 + trow1[x+2]*10) );
                    }
                    for(; x < w; ++x) {
                        dst[drow++] = ( (trow0[x+2] - trow0[x]) );
                        dst[drow++] = ( ((trow1[x+2] + trow1[x])*3 + trow1[x+1]*10) );
                    }
                }
            }
        };
    })();

    global.imgproc = imgproc;

})(jsfeat);
