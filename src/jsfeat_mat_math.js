/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 */

(function(global) {
    "use strict";
    //

    var matmath = (function() {
        
        return {
            transpose: function(At, A) {
                var i=0,j=0,nrows=A.rows,ncols=A.cols;
                var Ai=0,Ati=0,pAt=0;
                var ad=A.data,atd=At.data;

                for (; i < nrows; Ati += 1, Ai += ncols, i++) {
                    pAt = Ati;
                    for (j = 0; j < ncols; pAt += nrows, j++) atd[pAt] = ad[Ai+j];
                }
            },

            // C = A * B
            multiply: function(C, A, B) {
                var i=0,j=0,k=0;
                var Ap=0,pA=0,pB=0,p_B=0,Cp=0;
                var ncols=A.cols,nrows=A.rows,mcols=B.cols;
                var ad=A.data,bd=B.data,cd=C.data;
                var sum=0.0;

                for (; i < nrows; Ap += ncols, i++) {
                    for (p_B = 0, j = 0; j < mcols; Cp++, p_B++, j++) {
                        pB = p_B;
                        pA = Ap;
                        sum = 0.0;
                        for (k = 0; k < ncols; pA++, pB += mcols, k++) {
                            sum += ad[pA] * bd[pB];
                        }
                        cd[Cp] = sum;
                    }
                }
            },

            // C = A * B'
            multiply_ABt: function(C, A, B) {
                var i=0,j=0,k=0;
                var Ap=0,pA=0,pB=0,Cp=0;
                var ncols=A.cols,nrows=A.rows,mrows=B.rows;
                var ad=A.data,bd=B.data,cd=C.data;
                var sum=0.0;

                for (; i < nrows; Ap += ncols, i++) {
                    for (pB = 0, j = 0; j < mrows; Cp++, j++) {
                        pA = Ap;
                        sum = 0.0;
                        for (k = 0; k < ncols; pA++, pB++, k++) {
                            sum += ad[pA] * bd[pB];
                        }
                        cd[Cp] = sum;
                    }
                }
            },

            // C = A' * B
            multiply_AtB: function(C, A, B) {
                var i=0,j=0,k=0;
                var Ap=0,pA=0,pB=0,p_B=0,Cp=0;
                var ncols=A.cols,nrows=A.rows,mcols=B.cols;
                var ad=A.data,bd=B.data,cd=C.data;
                var sum=0.0;

                for (; i < ncols; Ap++, i++) {
                    for (p_B = 0, j = 0; j < mcols; Cp++, p_B++, j++) {
                        pB = p_B;
                        pA = Ap;
                        sum = 0.0;
                        for (k = 0; k < nrows; pA += ncols, pB += mcols, k++) {
                            sum += ad[pA] * bd[pB];
                        }
                        cd[Cp] = sum;
                    }
                }
            },

            // C = A * A'
            multiply_AAt: function(C, A) {
                var i=0,j=0,k=0;
                var pCdiag=0,p_A=0,pA=0,pB=0,pC=0,pCt=0;
                var ncols=A.cols,nrows=A.rows;
                var ad=A.data,cd=C.data;
                var sum=0.0;

                for (; i < nrows; pCdiag += nrows + 1, p_A = pA, i++) {
                    pC = pCdiag;
                    pCt = pCdiag;
                    pB = p_A; 
                    for (j = i; j < nrows; pC++, pCt += nrows, j++) {
                        pA = p_A;
                        sum = 0.0;
                        for (k = 0; k < ncols; k++) {
                            sum += ad[pA++] * ad[pB++];
                        }
                        cd[pC] = sum
                        cd[pCt] = sum;
                    }
                }
            },

            // C = A' * A
            multiply_AtA: function(C, A) {
                var i=0,j=0,k=0;
                var p_A=0,pA=0,pB=0,p_C=0,pC=0,p_CC=0;
                var ncols=A.cols,nrows=A.rows;
                var ad=A.data,cd=C.data;
                var sum=0.0;

                for (; i < ncols; p_C += ncols, i++) {
                    p_A = i;
                    p_CC = p_C + i;
                    pC = p_CC;
                    for (j = i; j < ncols; pC++, p_CC += ncols, j++) {
                        pA = p_A;
                        pB = j;
                        sum = 0.0;
                        for (k = 0; k < nrows; pA += ncols, pB += ncols, k++) {
                            sum += ad[pA] * ad[pB];
                        }
                        cd[pC] = sum
                        cd[p_CC] = sum;
                    }
                }
            }
        };

    })();

    global.matmath = matmath;

})(jsfeat);