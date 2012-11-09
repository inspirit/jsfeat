/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

// namespace ?
var jsfeat = jsfeat || { REVISION: 'ALPHA' };

// use typed arrays if possible
self.Int32Array = self.Int32Array || Array;
self.Uint32Array = self.Uint32Array || Array;
self.Uint8Array = self.Uint8Array || Array;
self.Float32Array = self.Float32Array || Array;
/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(global) {
    "use strict";
    //

    var U8_t = 2,
        U32_t = 4,
        I32_t = 8,
        F32_t = 16;

    var img_t = (function () {
        function img_t(w, h, elem_size) {
            this.width = w;
            this.height = h;
            this.buffer = new ArrayBuffer(w*h*elem_size);
            this.data_u8 = new Uint8Array(this.buffer);
            this.data_i32 = new Int32Array(this.buffer);
            this.data_f32 = new Float32Array(this.buffer);
        }
        return img_t;
    })();

    var img_pyr_t = (function () {

        function img_pyr_t(l) {
            this.levels = l;
            this.images = new Array(l);
            this.pyrdown = jsfeat.imgproc.pyrdown_fast_u8;
        }

        img_pyr_t.prototype.allocate = function(start_w, start_h, elem_size) {
            var i = this.levels;
            while(--i >= 0) {
                this.images[i] = new img_t(start_w >> i, start_h >> i, elem_size);
            }
        }

        img_pyr_t.prototype.build = function(src) {
            var i = 2, a = src, b = this.images[1];
            this.pyrdown(a.data_u8, a.width, a.height, b.data_u8);
            for(; i < this.levels; ++i) {
                a = b;
                b = this.images[i];
                this.pyrdown(a.data_u8, a.width, a.height, b.data_u8);
            }
        }

        return img_pyr_t;
    })();

    global.U8_t = U8_t;
    global.U32_t = U32_t;
    global.I32_t = I32_t;
    global.F32_t = F32_t;

    global.img_t = img_t;
    global.img_pyr_t = img_pyr_t;

})(jsfeat);
/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(global) {
    "use strict";
    //

    var math = (function() {

        return {
            get_gaussian_kernel: function(size, sigma, kernel, data_type) {
                var i=0,x=0.0,t=0.0,sigma_x=0.0,scale_2x=0.0;
                var sum = 0.0;
                var _kernel = new Float32Array(size);

                if((size&1) == 1 && size <= 7 && sigma <= 0) {
                    switch(size>>1) {
                        case 0:
                        _kernel[0] = 1.0;
                        sum = 1.0;
                        break;
                        case 1:
                        _kernel[0] = 0.25, _kernel[1] = 0.5, _kernel[2] = 0.25;
                        sum = 0.25+0.5+0.25;
                        break;
                        case 2:
                        _kernel[0] = 0.0625, _kernel[1] = 0.25, _kernel[2] = 0.375, 
                        _kernel[3] = 0.25, _kernel[4] = 0.0625;
                        sum = 0.0625+0.25+0.375+0.25+0.0625;
                        break;
                        case 3:
                        _kernel[0] = 0.03125, _kernel[1] = 0.109375, _kernel[2] = 0.21875, 
                        _kernel[3] = 0.28125, _kernel[4] = 0.21875, _kernel[5] = 0.109375, _kernel[6] = 0.03125;
                        sum = 0.03125+0.109375+0.21875+0.28125+0.21875+0.109375+0.03125;
                        break;
                    }
                } else {
                    sigma_x = sigma > 0 ? sigma : ((size-1)*0.5 - 1.0)*0.3 + 0.8;
                    scale_2x = -0.5/(sigma_x*sigma_x);

                    for( ; i < size; ++i )
                    {
                        x = i - (size-1)*0.5;
                        t = Math.exp(scale_2x*x*x);

                        _kernel[i] = t;
                        sum += t;
                    }
                }

                if(data_type & jsfeat.U8_t) {
                    // int based kernel
                    sum = 256.0/sum;
                    for (i = 0; i < size; ++i) {
                        kernel[i] = (_kernel[i] * sum + 0.5)|0;
                    }
                } else {
                    // classic kernel
                    sum = 1.0/sum;
                    for (i = 0; i < size; ++i) {
                        kernel[i] = _kernel[i] * sum;
                    }
                }
            },

            // The current implementation was derived from *BSD system qsort():
            // Copyright (c) 1992, 1993
            // The Regents of the University of California.  All rights reserved.
            qsort: function(array, low, high, cmp) {
                var isort_thresh = 7;
                var t;
                var sp = 0,left=0,right=0,i=0,n=0,m=0,ptr=0,ptr2=0,d=0;
                var left0=0,left1=0,right0=0,right1=0,pivot=0,a=0,b=0,c=0,swap_cnt=0;

                var stack = new Int32Array(48*2);

                if( (high-low+1) <= 1 ) return;

                stack[0] = low;
                stack[1] = high;

                while( sp >= 0 ) {
                
                    left = stack[sp<<1];
                    right = stack[(sp<<1)+1];
                    sp--;

                    for(;;) {
                        n = (right - left) + 1;

                        if( n <= isort_thresh ) {
                        //insert_sort:
                            for( ptr = left + 1; ptr <= right; ptr++ ) {
                                for( ptr2 = ptr; ptr2 > left && cmp(array[ptr2],array[ptr2-1]); ptr2--) {
                                    t = array[ptr2];
                                    array[ptr2] = array[ptr2-1];
                                    array[ptr2-1] = t;
                                }
                            }
                            break;
                        } else {
                            swap_cnt = 0;

                            left0 = left;
                            right0 = right;
                            pivot = left + (n>>1);

                            if( n > 40 ) {
                                d = n >> 3;
                                a = left, b = left + d, c = left + 2*d;
                                left = cmp(array[a], array[b]) ? (cmp(array[b], array[c]) ? b : (cmp(array[a], array[c]) ? c : a))
                                                  : (cmp(array[c], array[b]) ? b : (cmp(array[a], array[c]) ? a : c));

                                a = pivot - d, b = pivot, c = pivot + d;
                                pivot = cmp(array[a], array[b]) ? (cmp(array[b], array[c]) ? b : (cmp(array[a], array[c]) ? c : a))
                                                  : (cmp(array[c], array[b]) ? b : (cmp(array[a], array[c]) ? a : c));

                                a = right - 2*d, b = right - d, c = right;
                                right = cmp(array[a], array[b]) ? (cmp(array[b], array[c]) ? b : (cmp(array[a], array[c]) ? c : a))
                                                  : (cmp(array[c], array[b]) ? b : (cmp(array[a], array[c]) ? a : c));
                            }

                            a = left, b = pivot, c = right;
                            pivot = cmp(array[a], array[b]) ? (cmp(array[b], array[c]) ? b : (cmp(array[a], array[c]) ? c : a))   
                                               : (cmp(array[c], array[b]) ? b : (cmp(array[a], array[c]) ? a : c));
                            if( pivot != left0 ) {
                                t = array[pivot];
                                array[pivot] = array[left0];
                                array[left0] = t;
                                pivot = left0;
                            }
                            left = left1 = left0 + 1;
                            right = right1 = right0;

                            for(;;) {
                            
                                while( left <= right && !cmp(array[pivot], array[left]) ) {
                                    if( !cmp(array[left], array[pivot]) ) {
                                        if( left > left1 ) {
                                            t = array[left1];
                                            array[left1] = array[left];
                                            array[left] = t;
                                        }
                                        swap_cnt = 1;
                                        left1++;
                                    }
                                    left++;
                                }

                                while( left <= right && !cmp(array[right], array[pivot]) ) {
                                    if( !cmp(array[pivot], array[right]) ) {
                                        if( right < right1 ) {
                                            t = array[right1];
                                            array[right1] = array[right];
                                            array[right] = t;
                                        }
                                        swap_cnt = 1;
                                        right1--;
                                    }
                                    right--;
                                }

                                if( left > right ) break;
                                
                                t = array[left];
                                array[left] = array[right];
                                array[right] = t;
                                swap_cnt = 1;
                                left++;
                                right--;
                            }

                            if( swap_cnt == 0 ) {
                                left = left0, right = right0;
                                //goto insert_sort;
                                for( ptr = left + 1; ptr <= right; ptr++ ) {
                                    for( ptr2 = ptr; ptr2 > left && cmp(array[ptr2],array[ptr2-1]); ptr2--) {
                                        t = array[ptr2];
                                        array[ptr2] = array[ptr2-1];
                                        array[ptr2-1] = t;
                                    }
                                }
                                break;
                            }

                            n = Math.min( (left1 - left0), (left - left1) );
                            for( i = 0; i < n; i++ ) {
                                t = array[left0+i];
                                array[left0+i] = array[left+i-n];
                                array[left+i-n] = t;
                            }

                            n = Math.min( (right0 - right1), (right1 - right) );
                            for( i = 0; i < n; i++ ) {
                                t = array[left+i];
                                array[left+i] = array[right0+i-n+1];
                                array[right0+i-n+1] = t;
                            }
                            n = (left - left1);
                            m = (right1 - right);
                            if( n > 1 ) {
                                if( m > 1 ) {
                                    if( n > m ) {
                                        ++sp;
                                        stack[sp<<1] = left0;
                                        stack[(sp<<1)+1] = left0 + n - 1;
                                        left = right0 - m + 1, right = right0;
                                    } else {
                                        ++sp;
                                        stack[sp<<1] = right0 - m + 1;
                                        stack[(sp<<1)+1] = right0;
                                        left = left0, right = left0 + n - 1;
                                    }
                                } else {
                                    left = left0, right = left0 + n - 1;
                                }
                            }
                            else if( m > 1 )
                                left = right0 - m + 1, right = right0;
                            else
                                break;
                        }
                    }
                }
            },

            median: function(array, low, high) {
                var w;
                var middle=0,ll=0,hh=0,median=(low+high)>>1;
                for (;;) {
                    if (high <= low) return array[median];
                    if (high == (low + 1)) {
                        if (array[low] > array[high]) {
                            w = array[low];
                            array[low] = array[high];
                            array[high] = w;
                        }
                        return array[median];
                    }
                    middle = ((low + high) >> 1);
                    if (array[middle] > array[high]) {
                        w = array[middle];
                        array[middle] = array[high];
                        array[high] = w;
                    }
                    if (array[low] > array[high]) {
                        w = array[low];
                        array[low] = array[high];
                        array[high] = w;
                    }
                    if (array[middle] > array[low]) {
                        w = array[middle];
                        array[middle] = array[low];
                        array[low] = w;
                    }
                    ll = (low + 1);
                    w = array[middle];
                    array[middle] = array[ll];
                    array[ll] = w;
                    hh = high;
                    for (;;) {
                        do ++ll; while (array[low] > array[ll]);
                        do --hh; while (array[hh] > array[low]);
                        if (hh < ll) break;
                        w = array[ll];
                        array[ll] = array[hh];
                        array[hh] = w;
                    }
                    w = array[low];
                    array[low] = array[hh];
                    array[hh] = w;
                    if (hh <= median)
                        low = ll;
                    else if (hh >= median)
                        high = (hh - 1);
                }
                return 0;
            }
        };

    })();

    global.math = math;

})(jsfeat);
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

/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 * This is FAST corner detector, contributed to OpenCV by the author, Edward Rosten.
 */

/*
The references are:
 * Machine learning for high-speed corner detection,
   E. Rosten and T. Drummond, ECCV 2006
 * Faster and better: A machine learning approach to corner detection
   E. Rosten, R. Porter and T. Drummond, PAMI, 2009  
*/

(function(global) {
    "use strict";
    //
    var fast_corners = (function() {

        var offsets16 = new Int32Array([0, 3, 1, 3, 2, 2, 3, 1, 3, 0, 3, -1, 2, -2, 1, -3, 0, -3, -1, -3, -2, -2, -3, -1, -3, 0, -3, 1, -2, 2, -1, 3]);
        var offsets12 = new Int32Array([0, 2, 1, 2, 2, 1, 2, 0, 2, -1, 1, -2, 0, -2, -1, -2, -2, -1, -2, 0, -2, 1, -1, 2]);
        var offsets8 = new Int32Array([0, 1, 1, 1, 1, 0, 1, -1, 0, -1, -1, -1, -1, 0, -1, 1]);

        var threshold_tab = new Uint8Array(512);
        var pixel_off = new Int32Array(25);
        var score_diff = new Int32Array(25);

        // private functions
        var _cmp_offsets = function(pixel, step, pattern_size) {
            var k = 0;
            var offsets = pattern_size == 16 ? offsets16 : (pattern_size == 12 ? offsets12 : offsets8);
            for( ; k < pattern_size; ++k ) {
                pixel[k] = offsets[k<<1] + offsets[(k<<1)+1] * step;
            }
            for( ; k < 25; ++k ) {
                pixel[k] = pixel[k - pattern_size];
            }
        },

        _cmp_score_8 = function(src, off, pixel, d, threshold) {
            var N = 15, k = 0, v = src[off];
            var a0 = threshold,a=0,b0=0,b=0;

            for( ; k < N; ++k ) {
                d[k] = v - src[off+pixel[k]];
            }

            for( k = 0; k < 8; k += 2 ) {
                a = Math.min(d[k+1], d[k+2]);

                if( a <= a0 ) continue;

                a = Math.min(a, d[k+3]);
                a = Math.min(a, d[k+4]);
                a0 = Math.max(a0, Math.min(a, d[k]));
                a0 = Math.max(a0, Math.min(a, d[k+5]));
            }

            b0 = -a0;
            for( k = 0; k < 8; k += 2 ) {
                b = Math.max(d[k+1], d[k+2]);
                b = Math.max(b, d[k+3]);

                if( b >= b0 ) continue;

                b = Math.max(b, d[k+4]);
                b0 = Math.min(b0, Math.max(b, d[k]));
                b0 = Math.min(b0, Math.max(b, d[k+5]));
            }

            return -b0-1;
        },

        _cmp_score_12 = function(src, off, pixel, d, threshold) {
            var N = 19, k = 0, v = src[off];
            var a0 = threshold,a=0,b0=0,b=0;

            for( ; k < N; ++k ) {
                d[k] = v - src[off+pixel[k]];
            }

            for( k = 0; k < 12; k += 2 ) {
                a = Math.min(d[k+1], d[k+2]);

                if( a <= a0 ) continue;

                a = Math.min(a, d[k+3]);
                a = Math.min(a, d[k+4]);
                a = Math.min(a, d[k+5]);
                a = Math.min(a, d[k+6]);
                a0 = Math.max(a0, Math.min(a, d[k]));
                a0 = Math.max(a0, Math.min(a, d[k+7]));
            }

            b0 = -a0;
            for( k = 0; k < 12; k += 2 ) {
                b = Math.max(d[k+1], d[k+2]);
                b = Math.max(b, d[k+3]);
                b = Math.max(b, d[k+4]);

                if( b >= b0 ) continue;

                b = Math.max(b, d[k+5]);
                b = Math.max(b, d[k+6]);
                b0 = Math.min(b0, Math.max(b, d[k]));
                b0 = Math.min(b0, Math.max(b, d[k+7]));
            }

            return -b0-1;
        },

        _cmp_score_16 = function(src, off, pixel, d, threshold) {
            var N = 25, k = 0, v = src[off];
            var a0 = threshold,a=0,b0=0,b=0;

            for( ; k < N; ++k ) {
                d[k] = v - src[off+pixel[k]];
            }

            for( k = 0; k < 16; k += 2 ) {
                a = Math.min(d[k+1], d[k+2]);
                a = Math.min(a, d[k+3]);

                if( a <= a0 ) continue;

                a = Math.min(a, d[k+4]);
                a = Math.min(a, d[k+5]);
                a = Math.min(a, d[k+6]);
                a = Math.min(a, d[k+7]);
                a = Math.min(a, d[k+8]);
                a0 = Math.max(a0, Math.min(a, d[k]));
                a0 = Math.max(a0, Math.min(a, d[k+9]));
            }

            b0 = -a0;
            for( k = 0; k < 16; k += 2 ) {
                b = Math.max(d[k+1], d[k+2]);
                b = Math.max(b, d[k+3]);
                b = Math.max(b, d[k+4]);
                b = Math.max(b, d[k+5]);

                if( b >= b0 ) continue;
                b = Math.max(b, d[k+6]);
                b = Math.max(b, d[k+7]);
                b = Math.max(b, d[k+8]);
                b0 = Math.min(b0, Math.max(b, d[k]));
                b0 = Math.min(b0, Math.max(b, d[k+9]));
            }

            return -b0-1;
        };

        var _threshold = 20;

        return {
            set_threshold: function(threshold) {
                _threshold = Math.min(Math.max(threshold, 0), 255);
                for (var i = -255; i <= 255; ++i) {
                    threshold_tab[(i + 255)] = (i < -_threshold ? 1 : (i > _threshold ? 2 : 0));
                }
                return _threshold;
            },
            
            detect: function(img, w, h, corners, border, pattern_size) {
                if (typeof pattern_size === "undefined") { pattern_size = 16; }
                if (typeof border === "undefined") { border = 3; }

                var K = (pattern_size/2)|0, N = (pattern_size + K + 1)|0;
                var i=0, j=0, k=0, vt=0, x=0, m3=0;
                var buf = new Uint8Array(w*3);
                var cpbuf = new Int32Array((w+1)*3);
                var pixel = pixel_off;
                var sd = score_diff;
                var sy = Math.max(3, border);
                var ey = Math.min((h-2), (h-border));
                var sx = Math.max(3, border);
                var ex = Math.min((w - 3), (w - border));
                var _count = 0, corners_cnt = 0;
                var score_func = pattern_size == 16 ? _cmp_score_16 : (pattern_size == 12 ? _cmp_score_12 : _cmp_score_8);
                var thresh_tab = threshold_tab;
                var threshold = _threshold;

                var v=0,tab=0,d=0,ncorners=0,cornerpos=0,curr=0,ptr=0,prev=0,pprev=0;
                var jp1=0,jm1=0,score=0;

                _cmp_offsets(pixel, w, pattern_size);

                // local vars are faster?
                var pixel0 = pixel[0];
                var pixel1 = pixel[1];
                var pixel2 = pixel[2];
                var pixel3 = pixel[3];
                var pixel4 = pixel[4];
                var pixel5 = pixel[5];
                var pixel6 = pixel[6];
                var pixel7 = pixel[7];
                var pixel8 = pixel[8];
                var pixel9 = pixel[9];
                var pixel10 = pixel[10];
                var pixel11 = pixel[11];
                var pixel12 = pixel[12];
                var pixel13 = pixel[13];
                var pixel14 = pixel[14];
                var pixel15 = pixel[15];

                for(i = 0; i < w*3; ++i) {
                    buf[i] = 0;
                }

                for(i = sy; i < ey; ++i) {
                    ptr = ((i * w) + sx)|0;
                    m3 = (i - 3)%3;
                    curr = (m3*w)|0;
                    cornerpos = (m3*(w+1))|0;
                    for (j = 0; j < w; ++j) buf[curr+j] = 0;
                    ncorners = 0;
                    
                    if( i < (ey - 1) ) {
                        j = sx;
                        
                        for( ; j < ex; ++j, ++ptr ) {
                            v = img[ptr];
                            tab = ( - v + 255 );
                            d = ( thresh_tab[tab+img[ptr+pixel0]] | thresh_tab[tab+img[ptr+pixel8]] );
                            
                            if( d == 0 ) {
                                continue;
                            }
                            
                            d &= ( thresh_tab[tab+img[ptr+pixel2]] | thresh_tab[tab+img[ptr+pixel10]] );
                            d &= ( thresh_tab[tab+img[ptr+pixel4]] | thresh_tab[tab+img[ptr+pixel12]] );
                            d &= ( thresh_tab[tab+img[ptr+pixel6]] | thresh_tab[tab+img[ptr+pixel14]] );
                            
                            if( d == 0 ) {
                                continue;
                            }
                            
                            d &= ( thresh_tab[tab+img[ptr+pixel1]] | thresh_tab[tab+img[ptr+pixel9]] );
                            d &= ( thresh_tab[tab+img[ptr+pixel3]] | thresh_tab[tab+img[ptr+pixel11]] );
                            d &= ( thresh_tab[tab+img[ptr+pixel5]] | thresh_tab[tab+img[ptr+pixel13]] );
                            d &= ( thresh_tab[tab+img[ptr+pixel7]] | thresh_tab[tab+img[ptr+pixel15]] );
                            
                            if( d & 1 ) {
                                vt = (v - threshold);
                                _count = 0;
                                
                                for( k = 0; k < N; ++k ) {
                                    x = img[(ptr+pixel[k])];
                                    if(x < vt) {
                                        ++_count;
                                        if( _count > K ) {
                                            ++ncorners;
                                            cpbuf[cornerpos+ncorners] = j;
                                            buf[curr+j] = score_func(img, ptr, pixel, sd, threshold);
                                            break;
                                        }
                                    }
                                    else {
                                        _count = 0;
                                    }
                                }
                            }
                            
                            if( d & 2 ) {
                                vt = (v + threshold);
                                _count = 0;
                                
                                for( k = 0; k < N; ++k ) {
                                    x = img[(ptr+pixel[k])];
                                    if(x > vt) {
                                        ++_count;
                                        if( _count > K ) {
                                            ++ncorners;
                                            cpbuf[cornerpos+ncorners] = j;
                                            buf[curr+j] = score_func(img, ptr, pixel, sd, threshold);
                                            break;
                                        }
                                    }
                                    else {
                                        _count = 0;
                                    }
                                }
                            }
                        }
                    }
                    
                    cpbuf[cornerpos+w] = ncorners;
            
                    if ( i == sy ) {
                        continue;
                    }
                    
                    m3 = (i - 4 + 3)%3;
                    prev = (m3*w)|0;
                    cornerpos = (m3*(w+1))|0;
                    m3 = (i - 5 + 3)%3;
                    pprev = (m3*w)|0;

                    ncorners = cpbuf[cornerpos+w];
                    
                    for( k = 0; k < ncorners; ++k ) {
                        j = cpbuf[cornerpos+k];
                        jp1 = (j+1)|0;
                        jm1 = (j-1)|0;
                        score = buf[prev+j];
                        if( (score > buf[prev+jp1] && score > buf[prev+jm1] &&
                            score > buf[pprev+jm1] && score > buf[pprev+j] && score > buf[pprev+jp1] &&
                            score > buf[curr+jm1] && score > buf[curr+j] && score > buf[curr+jp1]) ) {
                            // save corner
                            corners[corners_cnt*3] = j;
                            corners[corners_cnt*3+1] = (i-1);
                            corners[corners_cnt*3+2] = score;
                            corners_cnt++;
                        }
                    }
                } // y loop
                return corners_cnt;
            }
        };
    })();

    global.fast_corners = fast_corners;
    fast_corners.set_threshold(20); // set default

})(jsfeat);
/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 * this code is a rewrite from OpenCV's Lucas-Kanade optical flow implementation
 */

(function(global) {
    "use strict";
    //
    var optical_flow_lk = (function() {

        // short link to shar deriv
        var shar_deriv = jsfeat.imgproc.shar_derivatives;

        return {
            track: function(prev_pyr, curr_pyr, prev_xy, curr_xy, count, win_size, max_iter, status, eps, min_eigen_threshold) {
                if (typeof max_iter === "undefined") { max_iter = 30; }
                if (typeof status === "undefined") { status = new Uint8Array(count); }
                if (typeof eps === "undefined") { eps = 0.01; }
                if (typeof min_eigen_threshold === "undefined") { min_eigen_threshold = 0.0001; }

                var half_win = (win_size-1)*0.5;
                var win_area = (win_size*win_size)|0;
                var win_area2 = win_area << 1;
                var prev_imgs = prev_pyr.images, next_imgs = curr_pyr.images;
                var img_prev=prev_imgs[0].data_u8,img_next=next_imgs[0].data_u8;
                var w0 = prev_imgs[0].width, h0 = prev_imgs[0].height,lw=0,lh=0;

                var iwin_buf = new Int32Array(win_area);
                var deriv_iwin = new Int32Array(win_area2);
                var deriv_lev = new Int32Array((w0*h0)*2);

                var dstep=0,src=0,dsrc=0,iptr=0,diptr=0,jptr=0;
                var lev_sc=0.0,prev_x=0.0,prev_y=0.0,next_x=0.0,next_y=0.0;
                var prev_delta_x=0.0,prev_delta_y=0.0,delta_x=0.0,delta_y=0.0;
                var iprev_x=0,iprev_y=0,inext_x=0,inext_y=0;
                var i=0,j=0,x=0,y=0,level=0,ptid=0,iter=0;
                var brd_tl=0,brd_r=0,brd_b=0;
                var a=0.0,b=0.0,b1=0.0,b2=0.0;

                // fixed point math
                var W_BITS14 = 14;
                var W_BITS4 = 14;
                var W_BITS1m5 = W_BITS4 - 5;
                var W_BITS1m51 = (1 << ((W_BITS1m5) - 1));
                var W_BITS14_ = (1 << W_BITS14);
                var W_BITS41 = (1 << ((W_BITS4) - 1));
                var FLT_SCALE = 1.0/(1 << 20);
                var iw00=0,iw01=0,iw10=0,iw11=0,ival=0,ixval=0,iyval=0;
                var A11=0.0,A12=0.0,A22=0.0,D=0.0,min_eig=0.0;

                var FLT_EPSILON = 0.00000011920929;
                eps *= eps;

                // reset status
                for(; i < count; ++i) {
                    status[i] = 1;
                }

                var max_level = (prev_pyr.levels - 1)|0;
                level = max_level;

                for(; level >= 0; --level) {
                    lev_sc = (1.0/(1 << level));
                    lw = w0 >> level;
                    lh = h0 >> level;
                    dstep = lw << 1;
                    img_prev = prev_imgs[level].data_u8;
                    img_next = next_imgs[level].data_u8;
                    
                    brd_r = (lw - win_size)|0;
                    brd_b = (lh - win_size)|0;

                    // calculate level derivatives
                    shar_deriv(img_prev, deriv_lev, lw, lh);

                    // iterate through points
                    for(ptid = 0; ptid < count; ++ptid) {
                        i = ptid << 1;
                        j = i + 1;
                        prev_x = prev_xy[i]*lev_sc;
                        prev_y = prev_xy[j]*lev_sc;

                        if( level == max_level ) {
                            next_x = prev_x;
                            next_y = prev_y;
                        } else {
                            next_x = curr_xy[i]*2.0;
                            next_y = curr_xy[j]*2.0;
                        }
                        curr_xy[i] = next_x;
                        curr_xy[j] = next_y;

                        prev_x -= half_win;
                        prev_y -= half_win;
                        iprev_x = prev_x|0;
                        iprev_y = prev_y|0;

                        // border check
                        x = (iprev_x <= brd_tl)|(iprev_x >= brd_r)|(iprev_y <= brd_tl)|(iprev_y >= brd_b);
                        if( x != 0 ) {
                            if( level == 0 ) {
                                status[ptid] = 0;
                            }
                            continue;
                        }

                        a = prev_x - iprev_x;
                        b = prev_y - iprev_y;
                        iw00 = (((1.0 - a)*(1.0 - b)*W_BITS14_) + 0.5)|0;
                        iw01 = ((a*(1.0 - b)*W_BITS14_) + 0.5)|0;
                        iw10 = (((1.0 - a)*b*W_BITS14_) + 0.5)|0;
                        iw11 = (W_BITS14_ - iw00 - iw01 - iw10);

                        A11 = 0.0, A12 = 0.0, A22 = 0.0;

                        // extract the patch from the first image, compute covariation matrix of derivatives
                        for( y = 0; y < win_size; ++y ) {
                            src = ( (y + iprev_y)*lw + iprev_x )|0;
                            dsrc = src << 1;

                            iptr = (y*win_size)|0;
                            diptr = iptr << 1;
                            for(x = 0 ; x < win_size; ++x, ++src, ++iptr, dsrc += 2) {
                                ival = ( (img_prev[src])*iw00 + (img_prev[src+1])*iw01 +
                                        (img_prev[src+lw])*iw10 + (img_prev[src+lw+1])*iw11 );
                                ival = (((ival) + W_BITS1m51) >> (W_BITS1m5));

                                ixval = ( deriv_lev[dsrc]*iw00 + deriv_lev[dsrc+2]*iw01 +
                                        deriv_lev[dsrc+dstep]*iw10 + deriv_lev[dsrc+dstep+2]*iw11 );
                                ixval = (((ixval) + W_BITS41) >> (W_BITS4));

                                iyval = ( deriv_lev[dsrc+1]*iw00 + deriv_lev[dsrc+3]*iw01 + deriv_lev[dsrc+dstep+1]*iw10 +
                                        deriv_lev[dsrc+dstep+3]*iw11 );
                                iyval = (((iyval) + W_BITS41) >> (W_BITS4));

                                iwin_buf[iptr] = ival;
                                deriv_iwin[diptr++] = ixval;
                                deriv_iwin[diptr++] = iyval;

                                A11 += ixval*ixval;
                                A12 += ixval*iyval;
                                A22 += iyval*iyval;
                            }
                        }

                        A11 *= FLT_SCALE; A12 *= FLT_SCALE; A22 *= FLT_SCALE;

                        D = A11*A22 - A12*A12;
                        min_eig = (A22 + A11 - Math.sqrt((A11-A22)*(A11-A22) + 4.0*A12*A12)) / win_area2;

                        if( min_eig < min_eigen_threshold || D < FLT_EPSILON )
                        {
                            if( level == 0 ) {
                                status[ptid] = 0;
                            }
                            continue;
                        }

                        D = 1.0/D;

                        next_x -= half_win;
                        next_y -= half_win;
                        prev_delta_x = 0.0;
                        prev_delta_y = 0.0;

                        for( iter = 0; iter < max_iter; ++iter ) {
                            inext_x = next_x|0;
                            inext_y = next_y|0;

                            x = (inext_x <= brd_tl)|(inext_x >= brd_r)|(inext_y <= brd_tl)|(inext_y >= brd_b);
                            if( x != 0 ) {
                                if( level == 0 ) {
                                    status[ptid] = 0;
                                }
                                break;
                            }

                            a = next_x - inext_x;
                            b = next_y - inext_y;
                            iw00 = (((1.0 - a)*(1.0 - b)*W_BITS14_) + 0.5)|0;
                            iw01 = ((a*(1.0 - b)*W_BITS14_) + 0.5)|0;
                            iw10 = (((1.0 - a)*b*W_BITS14_) + 0.5)|0;
                            iw11 = (W_BITS14_ - iw00 - iw01 - iw10);
                            b1 = 0.0, b2 = 0.0;

                            for( y = 0; y < win_size; ++y ) {
                                jptr = ( (y + inext_y)*lw + inext_x )|0;

                                iptr = (y*win_size)|0;
                                diptr = iptr << 1;
                                for( x = 0 ; x < win_size; ++x, ++jptr, ++iptr ) {
                                    ival = ( (img_next[jptr])*iw00 + (img_next[jptr+1])*iw01 +
                                            (img_next[jptr+lw])*iw10 + (img_next[jptr+lw+1])*iw11 );
                                    ival = (((ival) + W_BITS1m51) >> (W_BITS1m5));
                                    ival = (ival - iwin_buf[iptr]);

                                    b1 += ival * deriv_iwin[diptr++];
                                    b2 += ival * deriv_iwin[diptr++];
                                }
                            }

                            b1 *= FLT_SCALE;
                            b2 *= FLT_SCALE;

                            delta_x = ((A12*b2 - A22*b1) * D);
                            delta_y = ((A12*b1 - A11*b2) * D);

                            next_x += delta_x;
                            next_y += delta_y;
                            curr_xy[i] = next_x + half_win;
                            curr_xy[j] = next_y + half_win;

                            if( delta_x*delta_x + delta_y*delta_y <= eps ) {
                                break;
                            }

                            if( iter > 0 && Math.abs(delta_x + prev_delta_x) < 0.01 &&
                                            Math.abs(delta_y + prev_delta_y) < 0.01 ) {
                                curr_xy[i] -= delta_x*0.5;
                                curr_xy[j] -= delta_y*0.5;
                                break;
                            }

                            prev_delta_x = delta_x;
                            prev_delta_y = delta_y;
                        }
                    } // points loop
                } // levels loop
            }
        };
    })();

    global.optical_flow_lk = optical_flow_lk;

})(jsfeat);
