/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 */

(function(global) {
    "use strict";
    //

    var matmath = (function() {
        
        return {
            identity: function(M, value) {
                if (typeof value === "undefined") { value=1; }
                var src=M.data;
                var rows=M.rows, cols=M.cols, cols_1=(cols+1)|0;
                var len = rows * cols;
                var k = len;
                while(--len >= 0) src[len] = 0.0;
                len = k;
                k = 0;
                while(k < len)  {
                    src[k] = value;
                    k = k + cols_1;
                }
            },

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
            },

            // various small matrix operations
            identity_3x3: function(M, value) {
                if (typeof value === "undefined") { value=1; }
                var dt=M.data;
                dt[0] = dt[4] = dt[8] = value;
                dt[1] = dt[2] = dt[3] = 0;
                dt[5] = dt[6] = dt[7] = 0;
            },

            invert_3x3: function(from, to) {
                var A = from.data, invA = to.data;
                var t1 = A[4];
                var t2 = A[8];
                var t4 = A[5];
                var t5 = A[7];
                var t8 = A[0];

                var t9 = t8*t1;
                var t11 = t8*t4;
                var t13 = A[3];
                var t14 = A[1];
                var t15 = t13*t14;
                var t17 = A[2];
                var t18 = t13*t17;
                var t20 = A[6];
                var t21 = t20*t14;
                var t23 = t20*t17;
                var t26 = 1.0/(t9*t2-t11*t5-t15*t2+t18*t5+t21*t4-t23*t1);
                invA[0] = (t1*t2-t4*t5)*t26;
                invA[1] = -(t14*t2-t17*t5)*t26;
                invA[2] = -(-t14*t4+t17*t1)*t26;
                invA[3] = -(t13*t2-t4*t20)*t26;
                invA[4] = (t8*t2-t23)*t26;
                invA[5] = -(t11-t18)*t26;
                invA[6] = -(-t13*t5+t1*t20)*t26;
                invA[7] = -(t8*t5-t21)*t26;
                invA[8] = (t9-t15)*t26;
            },
            // C = A * B
            multiply_3x3: function(C, A, B) {
                var Cd=C.data, Ad=A.data, Bd=B.data;
                var m1_0 = Ad[0], m1_1 = Ad[1], m1_2 = Ad[2];
                var m1_3 = Ad[3], m1_4 = Ad[4], m1_5 = Ad[5];
                var m1_6 = Ad[6], m1_7 = Ad[7], m1_8 = Ad[8];

                var m2_0 = Bd[0], m2_1 = Bd[1], m2_2 = Bd[2];
                var m2_3 = Bd[3], m2_4 = Bd[4], m2_5 = Bd[5];
                var m2_6 = Bd[6], m2_7 = Bd[7], m2_8 = Bd[8];

                Cd[0] = m1_0 * m2_0 + m1_1 * m2_3 + m1_2 * m2_6;
                Cd[1] = m1_0 * m2_1 + m1_1 * m2_4 + m1_2 * m2_7;
                Cd[2] = m1_0 * m2_2 + m1_1 * m2_5 + m1_2 * m2_8;
                Cd[3] = m1_3 * m2_0 + m1_4 * m2_3 + m1_5 * m2_6;
                Cd[4] = m1_3 * m2_1 + m1_4 * m2_4 + m1_5 * m2_7;
                Cd[5] = m1_3 * m2_2 + m1_4 * m2_5 + m1_5 * m2_8;
                Cd[6] = m1_6 * m2_0 + m1_7 * m2_3 + m1_8 * m2_6;
                Cd[7] = m1_6 * m2_1 + m1_7 * m2_4 + m1_8 * m2_7;
                Cd[8] = m1_6 * m2_2 + m1_7 * m2_5 + m1_8 * m2_8;
            },

            mat3x3_determinant: function(M) {
                var md=M.data;
                return  md[0] * md[4] * md[8] -
                        md[0] * md[5] * md[7] -
                        md[3] * md[1] * md[8] +
                        md[3] * md[2] * md[7] +
                        md[6] * md[1] * md[5] -
                        md[6] * md[2] * md[4];
            },

            determinant_3x3: function(M11, M12, M13, 
                                      M21, M22, M23, 
                                      M31, M32, M33) {
                return  M11 * M22 * M33 - M11 * M23 * M32 -
                          M21 * M12 * M33 + M21 * M13 * M32 +
                          M31 * M12 * M23 - M31 * M13 * M22;
            }
        };

    })();

    global.matmath = matmath;

})(jsfeat);