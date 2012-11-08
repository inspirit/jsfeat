/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(global) {
    "use strict";
    //

    var math = (function() {

        return {
            imin: function(value0, value1) {
                return value1 ^ ((value0 ^ value1) & -(value0 < value1));
            },
            imax: function(value0, value1) {
                return value0 ^ ((value0 ^ value1) & -(value0 < value1));
            },
            imod3: function(a) {
                a = ((a >> 16) + (a & 0xFFFF)); /* sum base 2**16 digits a <= 0x1FFFE */
                a = ((a >>  8) + (a & 0xFF));   /* sum base 2**8 digits a <= 0x2FD */
                a = ((a >>  4) + (a & 0xF));    /* sum base 2**4 digits a <= 0x3C; worst case 0x3B */
                a = ((a >>  4) + (a & 0xF));    /* sum base 2**4 digits a <= 0x3C; worst case 0x3B */
                a = ((a >>  2) + (a & 0x3));    /* sum base 2**2 digits a <= 0x1D; worst case 0x1B */
                a = ((a >>  2) + (a & 0x3));    /* sum base 2**2 digits a <= 0x9; worst case 0x7 */
                //a = ((a >>  2) + (a & 0x3));    /* sum base 2**2 digits a <= 0x4 */
                //if (a > 2) a = (a - 3);
                return (a - 3*(a > 2));
            },
            fabs: function (value) {
                return value * (1.0 - ((value < 0.0) << 1));
            }
        };

    })();

    global.math = math;

})(jsfeat);