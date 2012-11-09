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