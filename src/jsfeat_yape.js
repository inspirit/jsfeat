/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 * Copyright 2007 Computer Vision Lab,
 * Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland.
 */

(function(global) {
    "use strict";
    //

    var yape = (function() {

        var precompute_directions = function(step, dirs, R) {
            var i = 0;
            var x, y;

            x = R;
            for(y = 0; y < x; y++, i++)
            {
                x = (Math.sqrt((R * R - y * y)) + 0.5)|0;
                dirs[i] = (x + step * y);
            }
            for(x-- ; x < y && x >= 0; x--, i++)
            {
                y = (Math.sqrt((R * R - x * x)) + 0.5)|0;
                dirs[i] = (x + step * y);
            }
            for( ; -x < y; x--, i++)
            {
                y = (Math.sqrt((R * R - x * x)) + 0.5)|0;
                dirs[i] = (x + step * y);
            }
            for(y-- ; y >= 0; y--, i++)
            {
                x = (-Math.sqrt((R * R - y * y)) - 0.5)|0;
                dirs[i] = (x + step * y);
            }
            for(; y > x; y--, i++)
            {
                x = (-Math.sqrt((R * R - y * y)) - 0.5)|0;
                dirs[i] = (x + step * y);
            }
            for(x++ ; x <= 0; x++, i++)
            {
                y = (-Math.sqrt((R * R - x * x)) - 0.5)|0;
                dirs[i] = (x + step * y);
            }
            for( ; x < -y; x++, i++)
            {
                y = (-Math.sqrt((R * R - x * x)) - 0.5)|0;
                dirs[i] = (x + step * y);
            }
            for(y++ ; y < 0; y++, i++)
            {
                x = (Math.sqrt((R * R - y * y)) + 0.5)|0;
                dirs[i] = (x + step * y);
            }

            dirs[i] = dirs[0];
            dirs[i + 1] = dirs[1];
            return i;
        }

        var third_check = function (Sb, off, step) {
            var n = 0;
            if(Sb[off+1]   != 0) n++;
            if(Sb[off-1]   != 0) n++;
            if(Sb[off+step]   != 0) n++;
            if(Sb[off+step+1] != 0) n++;
            if(Sb[off+step-1] != 0) n++;
            if(Sb[off-step]   != 0) n++;
            if(Sb[off-step+1] != 0) n++;
            if(Sb[off-step-1] != 0) n++;

            return n;
        }

        var is_local_maxima = function(p, off, v, step, neighborhood) {
            var x, y;

            if (v > 0) {
                off -= step*neighborhood;
                for (y= -neighborhood; y<=neighborhood; ++y) {
                    for (x= -neighborhood; x<=neighborhood; ++x) {
                        if (p[off+x] > v) return false;
                    }
                    off += step;
                }
            } else {
                off -= step*neighborhood;
                for (y= -neighborhood; y<=neighborhood; ++y) {
                    for (x= -neighborhood; x<=neighborhood; ++x) {
                        if (p[off+x] < v) return false;
                    }
                    off += step;
                }
            }
            return true;
        }

        var perform_one_point = function(I, x, Scores, Im, Ip, dirs, opposite, dirs_nb) {
          var score = 0;
          var a = 0, b = (opposite - 1)|0;
          var A=0, B0=0, B1=0, B2=0;
          var state=0;

          // WE KNOW THAT NOT(A ~ I0 & B1 ~ I0):
          A = I[x+dirs[a]];
          if ((A <= Ip)) {
            if ((A >= Im)) { // A ~ I0
              B0 = I[x+dirs[b]];
              if ((B0 <= Ip)) {
                if ((B0 >= Im)) { Scores[x] = 0; return; }
                else {
                  b++; B1 = I[x+dirs[b]];
                  if ((B1 > Ip)) {
                    b++; B2 = I[x+dirs[b]];
                    if ((B2 > Ip)) state = 3;
                    else if ((B2 < Im)) state = 6;
                    else { Scores[x] = 0; return; } // A ~ I0, B2 ~ I0
                  }
                  else/* if ((B1 < Im))*/ {
                    b++; B2 = I[x+dirs[b]];
                    if ((B2 > Ip)) state = 7;
                    else if ((B2 < Im)) state = 2;
                    else { Scores[x] = 0; return; } // A ~ I0, B2 ~ I0
                  }
                  //else { Scores[x] = 0; return; } // A ~ I0, B1 ~ I0
                }
              }
              else { // B0 < I0
                b++; B1 = I[x+dirs[b]];
                if ((B1 > Ip)) {
                  b++; B2 = I[x+dirs[b]];
                  if ((B2 > Ip)) state = 3;
                  else if ((B2 < Im)) state = 6;
                  else { Scores[x] = 0; return; } // A ~ I0, B2 ~ I0
                }
                else if ((B1 < Im)) {
                  b++; B2 = I[x+dirs[b]];
                  if ((B2 > Ip)) state = 7;
                  else if ((B2 < Im)) state = 2;
                  else { Scores[x] = 0; return; } // A ~ I0, B2 ~ I0
                }
                else { Scores[x] = 0; return; } // A ~ I0, B1 ~ I0
              }
            }
            else { // A > I0
              B0 = I[x+dirs[b]];
              if ((B0 > Ip)) { Scores[x] = 0; return; }
                b++; B1 = I[x+dirs[b]];
              if ((B1 > Ip)) { Scores[x] = 0; return; }
                b++; B2 = I[x+dirs[b]];
              if ((B2 > Ip)) { Scores[x] = 0; return; }
                state = 1;
            }
          }
          else // A < I0
          {
            B0 = I[x+dirs[b]];
            if ((B0 < Im)) { Scores[x] = 0; return; }
              b++; B1 = I[x+dirs[b]];
            if ((B1 < Im)) { Scores[x] = 0; return; }
              b++; B2 = I[x+dirs[b]];
            if ((B2 < Im)) { Scores[x] = 0; return; }
              state = 0;
          }

          for(a = 1; a <= opposite; a++)
          {
            A = I[x+dirs[a]];

            switch(state)
            {
            case 0:
              if ((A > Ip)) {
                B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 < Im)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 0; break; };
              }
              if ((A < Im)) {
                if ((B1 > Ip)) { Scores[x] = 0; return; }
                  if ((B2 > Ip)) { Scores[x] = 0; return; }
                    B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 > Ip)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 8; break; };
              } 
              // A ~ I0
              if ((B1 <= Ip)) { Scores[x] = 0; return; }
                if ((B2 <= Ip)) { Scores[x] = 0; return; }
                  B1 = B2; b++; B2 = I[x+dirs[b]];
              if ((B2 > Ip)) { score -= A + B1; state = 3; break; };
              if ((B2 < Im)) { score -= A + B1; state = 6; break; };
              { Scores[x] = 0; return; }

            case 1:
              if ((A < Im)) {
                B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 > Ip)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 1; break; };
              }
              if ((A > Ip)) {
                if ((B1 < Im)) { Scores[x] = 0; return; }
                  if ((B2 < Im)) { Scores[x] = 0; return; }
                    B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 < Im)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 9; break; };
              }
              // A ~ I0
              if ((B1 >= Im)) { Scores[x] = 0; return; }
                if ((B2 >= Im)) { Scores[x] = 0; return; }
                  B1 = B2; b++; B2 = I[x+dirs[b]];
              if ((B2 < Im)) { score -= A + B1; state = 2; break; };
              if ((B2 > Ip)) { score -= A + B1; state = 7; break; };
              { Scores[x] = 0; return; }

            case 2:
              if ((A > Ip)) { Scores[x] = 0; return; }
                B1 = B2; b++; B2 = I[x+dirs[b]];
              if ((A < Im))
              {
                if ((B2 > Ip)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 4; break; };
              } 
              // A ~ I0
              if ((B2 > Ip)) { score -= A + B1; state = 7; break; };
              if ((B2 < Im)) { score -= A + B1; state = 2; break; };
              { Scores[x] = 0; return; } // A ~ I0, B2 ~ I0

            case 3:
              if ((A < Im)) { Scores[x] = 0; return; }
                B1 = B2; b++; B2 = I[x+dirs[b]];
              if ((A > Ip)) {
                if ((B2 < Im)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 5; break; };
              }
              // A ~ I0
              if ((B2 > Ip)) { score -= A + B1; state = 3; break; };
              if ((B2 < Im)) { score -= A + B1; state = 6; break; };
              { Scores[x] = 0; return; }

            case 4:
              if ((A > Ip)) { Scores[x] = 0; return; }
                if ((A < Im)) {
                  B1 = B2; b++; B2 = I[x+dirs[b]];
                  if ((B2 > Ip)) { Scores[x] = 0; return; }
                    { score -= A + B1; state = 1; break; };
                }
                if ((B2 >= Im)) { Scores[x] = 0; return; }
                  B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 < Im)) { score -= A + B1; state = 2; break; };
                if ((B2 > Ip)) { score -= A + B1; state = 7; break; };
                { Scores[x] = 0; return; }

            case 5:
              if ((A < Im)) { Scores[x] = 0; return; }
                if ((A > Ip)) {
                  B1 = B2; b++; B2 = I[x+dirs[b]];
                  if ((B2 < Im)) { Scores[x] = 0; return; }
                    { score -= A + B1; state = 0; break; };
                }
                // A ~ I0
                if ((B2 <= Ip)) { Scores[x] = 0; return; }
                  B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 > Ip)) { score -= A + B1; state = 3; break; };
                if ((B2 < Im)) { score -= A + B1; state = 6; break; };
                { Scores[x] = 0; return; }

            case 7:
              if ((A > Ip)) { Scores[x] = 0; return; }
                if ((A < Im)) { Scores[x] = 0; return; }
                  B1 = B2; b++; B2 = I[x+dirs[b]];
              // A ~ I0
              if ((B2 > Ip)) { score -= A + B1; state = 3; break; };
              if ((B2 < Im)) { score -= A + B1; state = 6; break; };
              { Scores[x] = 0; return; } // A ~ I0, B2 ~ I0

            case 6:
              if ((A > Ip)) { Scores[x] = 0; return; }
                if ((A < Im)) { Scores[x] = 0; return; }
                  B1 = B2; b++; B2 = I[x+dirs[b]];
              // A ~ I0
              if ((B2 < Im)) { score -= A + B1; state = 2; break; };
              if ((B2 > Ip)) { score -= A + B1; state = 7; break; };
              { Scores[x] = 0; return; } // A ~ I0, B2 ~ I0

            case 8:
              if ((A > Ip)) {
                if ((B2 < Im)) { Scores[x] = 0; return; }
                  B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 < Im)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 9; break; };
              }
              if ((A < Im)) {
                B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 > Ip)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 1; break; };
              }
              { Scores[x] = 0; return; }

            case 9:
              if ((A < Im)) {
                if ((B2 > Ip)) { Scores[x] = 0; return; }
                  B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 > Ip)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 8; break; };
              }
              if ((A > Ip)) {
                B1 = B2; b++; B2 = I[x+dirs[b]];
                if ((B2 < Im)) { Scores[x] = 0; return; }
                  { score -= A + B1; state = 0; break; };
              }
              { Scores[x] = 0; return; }

            default:
              //"PB default";
              break;
            } // switch(state)
          } // for(a...)

          Scores[x] = (score + dirs_nb * I[x]);
        }

        var lev_table_t = (function () {
            function lev_table_t(w, h, r) {
                this.dirs = new Int32Array(1024);
                this.dirs_count = precompute_directions(w, this.dirs, r)|0;
                this.scores = new Int32Array(w*h);
                this.radius = r|0;
            }
            return lev_table_t;
        })();
        
        return {

            level_tables: [],
            tau: 7,

            init: function(width, height, radius, pyramid_levels) {
                if (typeof pyramid_levels === "undefined") { pyramid_levels = 1; }
                var i;
                radius = Math.min(radius, 7);
                radius = Math.max(radius, 3);
                for(i = 0; i < pyramid_levels; ++i) {
                    this.level_tables[i] = new lev_table_t(width>>i, height>>i, radius);
                }
            },

            detect: function(src, points, border) {
                if (typeof border === "undefined") { border = 4; }
                var t = this.level_tables[0];
                var R = t.radius|0, Rm1 = (R-1)|0;
                var dirs = t.dirs;
                var dirs_count = t.dirs_count|0;
                var opposite = dirs_count >> 1;
                var img = src.data, w=src.cols|0, h=src.rows|0,hw=w>>1;
                var scores = t.scores;
                var x=0,y=0,row=0,rowx=0,ip=0,im=0,abs_score=0, score=0;
                var tau = this.tau|0;
                var number_of_points = 0, pt;

                var sx = Math.max(R+1, border)|0;
                var sy = Math.max(R+1, border)|0;
                var ex = Math.min(w-R-2, w-border)|0;
                var ey = Math.min(h-R-2, h-border)|0;

                row = (sy*w+sx)|0;
                for(y = sy; y < ey; ++y, row+=w) {
                    for(x = sx, rowx = row; x < ex; ++x, ++rowx) {
                        ip = img[rowx] + tau, im = img[rowx] - tau;

                        if (im<img[rowx+R] && img[rowx+R]<ip && im<img[rowx-R] && img[rowx-R]<ip) {
                            scores[rowx] = 0;
                        } else {
                            perform_one_point(img, rowx, scores, im, ip, dirs, opposite, dirs_count);
                        }
                    }
                }

                // local maxima
                row = (sy*w+sx)|0;
                for(y = sy; y < ey; ++y, row+=w) {
                    for(x = sx, rowx = row; x < ex; ++x, ++rowx) {
                        score = scores[rowx];
                        abs_score = Math.abs(score);
                        if(abs_score < 5) {
                            // if this pixel is 0, the next one will not be good enough. Skip it.
                            ++x, ++rowx;
                        } else {
                            if(third_check(scores, rowx, w) >= 3 && is_local_maxima(scores, rowx, score, hw, R)) {
                                pt = points[number_of_points];
                                pt.x = x, pt.y = y, pt.score = abs_score;
                                ++number_of_points;

                                x += Rm1, rowx += Rm1;
                            }
                        }
                    }
                }

                return number_of_points;
            }
        };

    })();

    global.yape = yape;

})(jsfeat);