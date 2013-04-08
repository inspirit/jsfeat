/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 */

(function(global) {
    "use strict";
    //

    var linalg = (function() {

        var swap = function(A, i0, i1, t) {
            t = A[i0];
            A[i0] = A[i1];
            A[i1] = t;
        }

        var hypot = function(a, b) {
            a = Math.abs(a);
            b = Math.abs(b);
            if( a > b ) {
                b /= a;
                return a*Math.sqrt(1.0 + b*b);
            }
            if( b > 0 ) {
                a /= b;
                return b*Math.sqrt(1.0 + a*a);
            }
            return 0.0;
        }

        var JacobiImpl = function(A, astep, W, V, vstep, n) {
            var eps = jsfeat.EPSILON;
            var i=0,j=0,k=0,m=0,l=0,idx=0,_in=0,_in2=0;
            var iters=0,max_iter=n*n*30;
            var mv=0.0,val=0.0,p=0.0,y=0.0,t=0.0,s=0.0,c=0.0,a0=0.0,b0=0.0;

            var indR_buff = jsfeat.cache.get_buffer(n<<2);
            var indC_buff = jsfeat.cache.get_buffer(n<<2);
            var indR = indR_buff.i32;
            var indC = indC_buff.i32;

            if(V) {
                for(; i < n; i++) {
                    k = i*vstep;
                    for(j = 0; j < n; j++) {
                        V[k + j] = 0.0;
                    }
                    V[k + i] = 1.0;
                }
            }

            for(k = 0; k < n; k++) {
                W[k] = A[(astep + 1)*k];
                if(k < n - 1) {
                    for(m = k+1, mv = Math.abs(A[astep*k + m]), i = k+2; i < n; i++) {
                        val = Math.abs(A[astep*k+i]);
                        if(mv < val)
                            mv = val, m = i;
                    }
                    indR[k] = m;
                }
                if(k > 0) {
                    for(m = 0, mv = Math.abs(A[k]), i = 1; i < k; i++) {
                        val = Math.abs(A[astep*i+k]);
                        if(mv < val)
                            mv = val, m = i;
                    }
                    indC[k] = m;
                }
            }

            if(n > 1) for( ; iters < max_iter; iters++) {
                // find index (k,l) of pivot p
                for(k = 0, mv = Math.abs(A[indR[0]]), i = 1; i < n-1; i++) {
                    val = Math.abs(A[astep*i + indR[i]]);
                    if( mv < val )
                        mv = val, k = i;
                }
                l = indR[k];
                for(i = 1; i < n; i++) {
                    val = Math.abs(A[astep*indC[i] + i]);
                    if( mv < val )
                        mv = val, k = indC[i], l = i;
                }
                
                p = A[astep*k + l];

                if(Math.abs(p) <= eps) break;

                y = (W[l] - W[k])*0.5;
                t = Math.abs(y) + hypot(p, y);
                s = hypot(p, t);
                c = t/s;
                s = p/s; t = (p/t)*p;
                if(y < 0)
                    s = -s, t = -t;
                A[astep*k + l] = 0;
                
                W[k] -= t;
                W[l] += t;
                
                // rotate rows and columns k and l
                for (i = 0; i < k; i++) {
                    _in = (astep * i + k);
                    _in2 = (astep * i + l);
                    a0 = A[_in];
                    b0 = A[_in2];
                    A[_in] = a0 * c - b0 * s;
                    A[_in2] = a0 * s + b0 * c;
                }
                for (i = (k + 1); i < l; i++) {
                    _in = (astep * k + i);
                    _in2 = (astep * i + l);
                    a0 = A[_in];
                    b0 = A[_in2];
                    A[_in] = a0 * c - b0 * s;
                    A[_in2] = a0 * s + b0 * c;
                }
                i = l + 1;
                _in = (astep * k + i);
                _in2 = (astep * l + i);
                for (; i < n; i++, _in++, _in2++) {
                    a0 = A[_in];
                    b0 = A[_in2];
                    A[_in] = a0 * c - b0 * s;
                    A[_in2] = a0 * s + b0 * c;
                }
                
                // rotate eigenvectors
                if (V) {
                    _in = vstep * k;
                    _in2 = vstep * l;
                    for (i = 0; i < n; i++, _in++, _in2++) {
                        a0 = V[_in];
                        b0 = V[_in2];
                        V[_in] = a0 * c - b0 * s;
                        V[_in2] = a0 * s + b0 * c;
                    }
                }
                
                for(j = 0; j < 2; j++) {
                    idx = j == 0 ? k : l;
                    if(idx < n - 1) {
                        for(m = idx+1, mv = Math.abs(A[astep*idx + m]), i = idx+2; i < n; i++) {
                            val = Math.abs(A[astep*idx+i]);
                            if( mv < val )
                                mv = val, m = i;
                        }
                        indR[idx] = m;
                    }
                    if(idx > 0) {
                        for(m = 0, mv = Math.abs(A[idx]), i = 1; i < idx; i++) {
                            val = Math.abs(A[astep*i+idx]);
                            if( mv < val )
                                mv = val, m = i;
                        }
                        indC[idx] = m;
                    }
                }
            }

            // sort eigenvalues & eigenvectors
            for(k = 0; k < n-1; k++) {
                m = k;
                for(i = k+1; i < n; i++) {
                    if(W[m] < W[i])
                        m = i;
                }
                if(k != m) {
                    swap(W, m, k, mv);
                    if(V) {
                        for(i = 0; i < n; i++) {
                            swap(V, vstep*m + i, vstep*k + i, mv);
                        }
                    }
                }
            }


            jsfeat.cache.put_buffer(indR_buff);
            jsfeat.cache.put_buffer(indC_buff);
        }

        var JacobiSVDImpl = function(At, astep, _W, Vt, vstep, m, n, n1) {
            var eps = jsfeat.EPSILON * 2.0;
            var minval = jsfeat.FLT_MIN;
            var i=0,j=0,k=0,iter=0,max_iter=Math.max(m, 30);
            var Ai=0,Aj=0,Vi=0,Vj=0,changed=0;
            var c=0.0, s=0.0, t=0.0;
            var t0=0.0,t1=0.0,sd=0.0,beta=0.0,gamma=0.0,delta=0.0,a=0.0,p=0.0,b=0.0;
            var seed = 0x1234;
            var val=0.0,val0=0.0,asum=0.0;

            var W_buff = jsfeat.cache.get_buffer(n<<3);
            var W = W_buff.f64;
            
            for(; i < n; i++) {
                for(k = 0, sd = 0; k < m; k++) {
                    t = At[i*astep + k];
                    sd += t*t;
                }
                W[i] = sd;
                
                if(Vt) {
                    for(k = 0; k < n; k++) {
                        Vt[i*vstep + k] = 0;
                    }
                    Vt[i*vstep + i] = 1;
                }
            }
            
            for(; iter < max_iter; iter++) {
                changed = 0;
                
                for(i = 0; i < n-1; i++) {
                    for(j = i+1; j < n; j++) {
                        Ai = (i*astep)|0, Aj = (j*astep)|0;
                        a = W[i], p = 0, b = W[j];
                        
                        k = 2;
                        p += At[Ai]*At[Aj];
                        p += At[Ai+1]*At[Aj+1];

                        for(; k < m; k++)
                            p += At[Ai+k]*At[Aj+k];
                        
                        if(Math.abs(p) <= eps*Math.sqrt(a*b)) continue;
                        
                        p *= 2.0;
                        beta = a - b, gamma = hypot(p, beta);
                        if( beta < 0 ) {
                            delta = (gamma - beta)*0.5;
                            s = Math.sqrt(delta/gamma);
                            c = (p/(gamma*s*2.0));
                        } else {
                            c = Math.sqrt((gamma + beta)/(gamma*2.0));
                            s = (p/(gamma*c*2.0));
                        }
                        
                        a=0.0, b=0.0;
                        
                        k = 2; // unroll
                        t0 = c*At[Ai] + s*At[Aj];
                        t1 = -s*At[Ai] + c*At[Aj];
                        At[Ai] = t0; At[Aj] = t1;
                        a += t0*t0; b += t1*t1;

                        t0 = c*At[Ai+1] + s*At[Aj+1];
                        t1 = -s*At[Ai+1] + c*At[Aj+1];
                        At[Ai+1] = t0; At[Aj+1] = t1;
                        a += t0*t0; b += t1*t1;

                        for( ; k < m; k++ )
                        {
                            t0 = c*At[Ai+k] + s*At[Aj+k];
                            t1 = -s*At[Ai+k] + c*At[Aj+k];
                            At[Ai+k] = t0; At[Aj+k] = t1;
                            
                            a += t0*t0; b += t1*t1;
                        }
                        
                        W[i] = a; W[j] = b;
                        
                        changed = 1;
                        
                        if(Vt) {
                            Vi = (i*vstep)|0, Vj = (j*vstep)|0;

                            k = 2;
                            t0 = c*Vt[Vi] + s*Vt[Vj];
                            t1 = -s*Vt[Vi] + c*Vt[Vj];
                            Vt[Vi] = t0; Vt[Vj] = t1;

                            t0 = c*Vt[Vi+1] + s*Vt[Vj+1];
                            t1 = -s*Vt[Vi+1] + c*Vt[Vj+1];
                            Vt[Vi+1] = t0; Vt[Vj+1] = t1;

                            for(; k < n; k++) {
                                t0 = c*Vt[Vi+k] + s*Vt[Vj+k];
                                t1 = -s*Vt[Vi+k] + c*Vt[Vj+k];
                                Vt[Vi+k] = t0; Vt[Vj+k] = t1;
                            }
                        }
                    }
                }
                if(changed == 0) break;
            }
            
            for(i = 0; i < n; i++) {
                for(k = 0, sd = 0; k < m; k++) {
                    t = At[i*astep + k];
                    sd += t*t;
                }
                W[i] = Math.sqrt(sd);
            }
            
            for(i = 0; i < n-1; i++) {
                j = i;
                for(k = i+1; k < n; k++) {
                    if(W[j] < W[k])
                        j = k;
                }
                if(i != j) {
                    swap(W, i, j, sd);
                    if(Vt) {
                        for(k = 0; k < m; k++) {
                            swap(At, i*astep + k, j*astep + k, t);
                        }
                        
                        for(k = 0; k < n; k++) {
                            swap(Vt, i*vstep + k, j*vstep + k, t);
                        }
                    }
                }
            }
            
            for(i = 0; i < n; i++) {
                _W[i] = W[i];
            }
            
            if(!Vt) {
                jsfeat.cache.put_buffer(W_buff);
                return;
            }

            for(i = 0; i < n1; i++) {

                sd = i < n ? W[i] : 0;
                
                while(sd <= minval) {
                    // if we got a zero singular value, then in order to get the corresponding left singular vector
                    // we generate a random vector, project it to the previously computed left singular vectors,
                    // subtract the projection and normalize the difference.
                    val0 = (1.0/m);
                    for(k = 0; k < m; k++) {
                        seed = (seed * 214013 + 2531011);
                        val = (((seed >> 16) & 0x7fff) & 256) != 0 ? val0 : -val0;
                        At[i*astep + k] = val;
                    }
                    for(iter = 0; iter < 2; iter++) {
                        for(j = 0; j < i; j++) {
                            sd = 0;
                            for(k = 0; k < m; k++) {
                                sd += At[i*astep + k]*At[j*astep + k];
                            }
                            asum = 0.0;
                            for(k = 0; k < m; k++) {
                                t = (At[i*astep + k] - sd*At[j*astep + k]);
                                At[i*astep + k] = t;
                                asum += Math.abs(t);
                            }
                            asum = asum ? 1.0/asum : 0;
                            for(k = 0; k < m; k++) {
                                At[i*astep + k] *= asum;
                            }
                        }
                    }
                    sd = 0;
                    for(k = 0; k < m; k++) {
                        t = At[i*astep + k];
                        sd += t*t;
                    }
                    sd = Math.sqrt(sd);
                }
                
                s = (1.0/sd);
                for(k = 0; k < m; k++) {
                    At[i*astep + k] *= s;
                }
            }

            jsfeat.cache.put_buffer(W_buff);
        }
        
        return {

            lu_solve: function(A, B) {
                var i=0,j=0,k=0,p=1,astep=A.cols;
                var ad=A.data, bd=B.data;
                var t,alpha,d,s;

                for(i = 0; i < astep; i++) {
                    k = i;                    
                    for(j = i+1; j < astep; j++) {
                        if(Math.abs(ad[j*astep + i]) > Math.abs(ad[k*astep+i])) {
                            k = j;
                        }
                    }
                    
                    if(Math.abs(ad[k*astep+i]) < jsfeat.EPSILON) {
                        return 0; // FAILED
                    }
                    
                    if(k != i) {
                        for(j = i; j < astep; j++ ) {
                            swap(ad, i*astep+j, k*astep+j, t);
                        }
                        
                        swap(bd, i, k, t);
                        p = -p;
                    }
                    
                    d = -1.0/ad[i*astep+i];
                    
                    for(j = i+1; j < astep; j++) {
                        alpha = ad[j*astep+i]*d;
                        
                        for(k = i+1; k < astep; k++) {
                            ad[j*astep+k] += alpha*ad[i*astep+k];
                        }
                        
                        bd[j] += alpha*bd[i];
                    }
                    
                    ad[i*astep+i] = -d;
                }
                
                for(i = astep-1; i >= 0; i--) {
                    s = bd[i];
                    for(k = i+1; k < astep; k++) {
                        s -= ad[i*astep+k]*bd[k];
                    }
                    bd[i] = s*ad[i*astep+i];
                }

                return 1; // OK
            },

            cholesky_solve: function(A, B) {
                var col=0,row=0,col2=0,cs=0,rs=0,i=0,j=0;
                var size = A.cols;
                var ad=A.data, bd=B.data;
                var val,inv_diag;

                for (col = 0; col < size; col++) {
                    inv_diag = 1.0;
                    cs = (col * size);
                    rs = cs;
                    for (row = col; row < size; row++)
                    {
                        // correct for the parts of cholesky already computed
                        val = ad[(rs+col)];
                        for (col2 = 0; col2 < col; col2++) {
                            val -= ad[(col2*size+col)] * ad[(rs+col2)];
                        }
                        if (row == col) {
                            // this is the diagonal element so don't divide
                            ad[(rs+col)] = val;
                            if(val == 0) {
                                return 0;
                            }
                            inv_diag = 1.0 / val;
                        } else {
                            // cache the value without division in the upper half
                            ad[(cs+row)] = val;
                            // divide my the diagonal element for all others
                            ad[(rs+col)] = val * inv_diag;
                        }
                        rs = (rs + size);
                    }
                }

                // first backsub through L
                cs = 0;
                for (i = 0; i < size; i++) {
                    val = bd[i];
                    for (j = 0; j < i; j++) {
                        val -= ad[(cs+j)] * bd[j];
                    }
                    bd[i] = val;
                    cs = (cs + size);
                }
                // backsub through diagonal
                cs = 0;
                for (i = 0; i < size; i++) {
                    bd[i] /= ad[(cs + i)];
                    cs = (cs + size);
                }
                // backsub through L Transpose
                i = (size-1);
                for (; i >= 0; i--) {
                    val = bd[i];
                    j = (i + 1);
                    cs = (j * size);
                    for (; j < size; j++) {
                        val -= ad[(cs + i)] * bd[j];
                        cs = (cs + size);
                    }
                    bd[i] = val;
                }

                return 1;
            },

            svd_decompose: function(A, W, U, V, options) {
                if (typeof options === "undefined") { options = 0; };
                var at=0,i=0,j=0,_m=A.rows,_n=A.cols,m=_m,n=_n;
                var dt = A.type | jsfeat.C1_t; // we only work with single channel

                if(m < n) {
                    at = 1;
                    i = m;
                    m = n;
                    n = i;
                }

                var a_buff = jsfeat.cache.get_buffer((m*m)<<3);
                var w_buff = jsfeat.cache.get_buffer(n<<3);
                var v_buff = jsfeat.cache.get_buffer((n*n)<<3);

                var a_mt = new jsfeat.matrix_t(m, m, dt, a_buff.data);
                var w_mt = new jsfeat.matrix_t(1, n, dt, w_buff.data);
                var v_mt = new jsfeat.matrix_t(n, n, dt, v_buff.data);

                if(at == 0) {
                    // transpose
                    jsfeat.matmath.transpose(a_mt, A);
                } else {
                    for(i = 0; i < _n*_m; i++) {
                        a_mt.data[i] = A.data[i];
                    }
                    for(; i < n*m; i++) {
                        a_mt.data[i] = 0;
                    }
                }

                JacobiSVDImpl(a_mt.data, m, w_mt.data, v_mt.data, n, m, n, m);

                if(W) {
                    for(i=0; i < n; i++) {
                        W.data[i] = w_mt.data[i];
                    }
                    for(; i < _n; i++) {
                        W.data[i] = 0;
                    }
                }

                if (at == 0) {
                    if(U && (options & jsfeat.SVD_U_T)) {
                        i = m*m;
                        while(--i >= 0) {
                            U.data[i] = a_mt.data[i];
                        }
                    } else if(U) {
                        jsfeat.matmath.transpose(U, a_mt);
                    }

                    if(V && (options & jsfeat.SVD_V_T)) {
                        i = n*n;
                        while(--i >= 0) {
                            V.data[i] = v_mt.data[i];
                        }
                    } else if(V) {
                        jsfeat.matmath.transpose(V, v_mt);
                    }
                } else {
                    if(U && (options & jsfeat.SVD_U_T)) {
                        i = n*n;
                        while(--i >= 0) {
                            U.data[i] = v_mt.data[i];
                        }
                    } else if(U) {
                        jsfeat.matmath.transpose(U, v_mt);
                    }

                    if(V && (options & jsfeat.SVD_V_T)) {
                        i = m*m;
                        while(--i >= 0) {
                            V.data[i] = a_mt.data[i];
                        }
                    } else if(V) {
                        jsfeat.matmath.transpose(V, a_mt);
                    }
                }

                jsfeat.cache.put_buffer(a_buff);
                jsfeat.cache.put_buffer(w_buff);
                jsfeat.cache.put_buffer(v_buff);

            },

            svd_solve: function(A, X, B) {
                var i=0,j=0,k=0;
                var pu=0,pv=0;
                var nrows=A.rows,ncols=A.cols;
                var sum=0.0,xsum=0.0,tol=0.0;
                var dt = A.type | jsfeat.C1_t;

                var u_buff = jsfeat.cache.get_buffer((nrows*nrows)<<3);
                var w_buff = jsfeat.cache.get_buffer(ncols<<3);
                var v_buff = jsfeat.cache.get_buffer((ncols*ncols)<<3);

                var u_mt = new jsfeat.matrix_t(nrows, nrows, dt, u_buff.data);
                var w_mt = new jsfeat.matrix_t(1, ncols, dt, w_buff.data);
                var v_mt = new jsfeat.matrix_t(ncols, ncols, dt, v_buff.data);

                var bd = B.data, ud = u_mt.data, wd = w_mt.data, vd = v_mt.data;

                this.svd_decompose(A, w_mt, u_mt, v_mt, 0);

                tol = jsfeat.EPSILON * wd[0] * ncols;

                for (; i < ncols; i++, pv += ncols) {
                    xsum = 0.0;
                    for(j = 0; j < ncols; j++) {
                        if(wd[j] > tol) {
                            for(k = 0, sum = 0.0, pu = 0; k < nrows; k++, pu += ncols) {
                                sum += ud[pu + j] * bd[k];
                            }
                            xsum += sum * vd[pv + j] / wd[j];
                        }
                    }
                    X.data[i] = xsum;
                }

                jsfeat.cache.put_buffer(u_buff);
                jsfeat.cache.put_buffer(w_buff);
                jsfeat.cache.put_buffer(v_buff);
            },

            svd_invert: function(Ai, A) {
                var i=0,j=0,k=0;
                var pu=0,pv=0,pa=0;
                var nrows=A.rows,ncols=A.cols;
                var sum=0.0,tol=0.0;
                var dt = A.type | jsfeat.C1_t;

                var u_buff = jsfeat.cache.get_buffer((nrows*nrows)<<3);
                var w_buff = jsfeat.cache.get_buffer(ncols<<3);
                var v_buff = jsfeat.cache.get_buffer((ncols*ncols)<<3);

                var u_mt = new jsfeat.matrix_t(nrows, nrows, dt, u_buff.data);
                var w_mt = new jsfeat.matrix_t(1, ncols, dt, w_buff.data);
                var v_mt = new jsfeat.matrix_t(ncols, ncols, dt, v_buff.data);

                var id = Ai.data, ud = u_mt.data, wd = w_mt.data, vd = v_mt.data;

                this.svd_decompose(A, w_mt, u_mt, v_mt, 0);

                tol = jsfeat.EPSILON * wd[0] * ncols;

                for (; i < ncols; i++, pv += ncols) {
                    for (j = 0, pu = 0; j < nrows; j++, pa++) {
                        for (k = 0, sum = 0.0; k < ncols; k++, pu++) {
                            if (wd[k] > tol) sum += vd[pv + k] * ud[pu] / wd[k];
                        }
                        id[pa] = sum;
                    }
                }

                jsfeat.cache.put_buffer(u_buff);
                jsfeat.cache.put_buffer(w_buff);
                jsfeat.cache.put_buffer(v_buff);
            },

            eigenVV: function(A, vects, vals) {
                var n=A.cols,i=n*n;
                var dt = A.type | jsfeat.C1_t;

                var a_buff = jsfeat.cache.get_buffer((n*n)<<3);
                var w_buff = jsfeat.cache.get_buffer(n<<3);
                var a_mt = new jsfeat.matrix_t(n, n, dt, a_buff.data);
                var w_mt = new jsfeat.matrix_t(1, n, dt, w_buff.data);

                while(--i >= 0) {
                    a_mt.data[i] = A.data[i];
                }

                JacobiImpl(a_mt.data, n, w_mt.data, vects ? vects.data : null, n, n);

                if(vals) {
                    while(--n >= 0) {
                        vals.data[n] = w_mt.data[n];
                    }
                }

                jsfeat.cache.put_buffer(a_buff);
                jsfeat.cache.put_buffer(w_buff);
            }

        };

    })();

    global.linalg = linalg;

})(jsfeat);