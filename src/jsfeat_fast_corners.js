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

        var threshold_tab = new Uint8Array(512);
        var pixel_off = new Int32Array(25);
        var score_diff = new Int32Array(25);

        // private functions
        var _cmp_offsets = function(pixel, step, pattern_size) {
            var k = 0;
            var offsets = offsets16;
            for( ; k < pattern_size; ++k ) {
                pixel[k] = offsets[k<<1] + offsets[(k<<1)+1] * step;
            }
            for( ; k < 25; ++k ) {
                pixel[k] = pixel[k - pattern_size];
            }
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
            
            detect: function(src, corners, border) {
                if (typeof border === "undefined") { border = 3; }

                var K = 8, N = 25;
                var img = src.data, w = src.cols, h = src.rows;
                var i=0, j=0, k=0, vt=0, x=0, m3=0;
                var buf_node = jsfeat.cache.get_buffer(3 * w);
                var cpbuf_node = jsfeat.cache.get_buffer(((w+1)*3)<<2);
                var buf = buf_node.u8;
                var cpbuf = cpbuf_node.i32;
                var pixel = pixel_off;
                var sd = score_diff;
                var sy = Math.max(3, border);
                var ey = Math.min((h-2), (h-border));
                var sx = Math.max(3, border);
                var ex = Math.min((w - 3), (w - border));
                var _count = 0, corners_cnt = 0, pt;
                var score_func = _cmp_score_16;
                var thresh_tab = threshold_tab;
                var threshold = _threshold;

                var v=0,tab=0,d=0,ncorners=0,cornerpos=0,curr=0,ptr=0,prev=0,pprev=0;
                var jp1=0,jm1=0,score=0;

                _cmp_offsets(pixel, w, 16);

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
                            pt = corners[corners_cnt];
                            pt.x = j, pt.y = (i-1), pt.score = score;
                            corners_cnt++;
                        }
                    }
                } // y loop
                jsfeat.cache.put_buffer(buf_node);
                jsfeat.cache.put_buffer(cpbuf_node);
                return corners_cnt;
            }
        };
    })();

    global.fast_corners = fast_corners;
    fast_corners.set_threshold(20); // set default

})(jsfeat);
