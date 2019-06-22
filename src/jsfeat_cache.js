/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 */

(function(global) {
    "use strict";
    //

    var cache = (function() {

        // very primitive array cache, still need testing if it helps
        // of course V8 has its own powerful cache sys but i'm not sure
        // it caches several multichannel 640x480 buffer creations each frame

        var _pool_node_t = (function () {
            function _pool_node_t(size_in_bytes) {
                this.next = null;
                this.data = new jsfeat.data_t(size_in_bytes);
                this.size = this.data.size;
                this.buffer = this.data.buffer;
                this.u8 = this.data.u8;
                this.i32 = this.data.i32;
                this.f32 = this.data.f32;
                this.f64 = this.data.f64;
            }
            _pool_node_t.prototype.resize = function(size_in_bytes) {
                delete this.data;
                this.data = new jsfeat.data_t(size_in_bytes);
                this.size = this.data.size;
                this.buffer = this.data.buffer;
                this.u8 = this.data.u8;
                this.i32 = this.data.i32;
                this.f32 = this.data.f32;
                this.f64 = this.data.f64;
            }
            return _pool_node_t;
        })();

        var _pool_head, _pool_tail;
        var _pool_size = 0;

        return {

            allocate: function(capacity, data_size) {
                _pool_head = _pool_tail = new _pool_node_t(data_size);
                for (var i = 0; i < capacity; ++i) {
                    var node = new _pool_node_t(data_size);
                    _pool_tail = _pool_tail.next = node;

                    _pool_size++;
                }
            },

            get_buffer: function(size_in_bytes) {
                // assume we have enough free nodes
                var node = _pool_head;
                _pool_head = _pool_head.next;
                _pool_size--;

                if(size_in_bytes > node.size) {
                    node.resize(size_in_bytes);
                }

                return node;
            },

            put_buffer: function(node) {
                _pool_tail = _pool_tail.next = node;
                _pool_size++;
            }
        };
    })();

    global.cache = cache;
    // for now we dont need more than 30 buffers
    // if having cache sys really helps we can add auto extending sys
    cache.allocate(30, 640*4);

})(jsfeat);
