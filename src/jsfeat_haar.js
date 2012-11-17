/**
 * @author Eugene Zatepyakin / http://inspirit.ru/
 *
 * this code is a rewrite from https://github.com/mtschirs/js-objectdetect implementation
 * @author Martin Tschirsich / http://www.tu-darmstadt.de/~m_t
 */

(function(global) {
    "use strict";
    //
    var haar = (function() {

        var _group_func = function(r1, r2) {
            var distance = (r1.width * 0.25 + 0.5)|0;

            return r2.x <= r1.x + distance &&
                   r2.x >= r1.x - distance &&
                   r2.y <= r1.y + distance &&
                   r2.y >= r1.y - distance &&
                   r2.width <= (r1.width * 1.5 + 0.5)|0 &&
                   (r2.width * 1.5 + 0.5)|0 >= r1.width;
        }
        
        return {

            edges_density: 0.07,

            detect_single_scale: function(int_sum, int_sqsum, int_tilted, int_canny_sum, width, height, scale, classifier) {
                var win_w = (classifier.size[0] * scale)|0,
                    win_h = (classifier.size[1] * scale)|0,
                    step_x = (0.5 * scale + 1.5)|0,
                    step_y = step_x;
                var i,j,k,x,y,ex=(width-win_w)|0,ey=(height-win_h)|0;
                var w1=(width+1)|0,edge_dens,mean,variance,std;
                var inv_area = 1.0 / (win_w * win_h);
                var stages,stage,trees,tree,sn,tn,fn,found=true,stage_thresh,stage_sum,tree_sum,feature,features;
                var fi_a,fi_b,fi_c,fi_d,fw,fh;

                var ii_a=0,ii_b=win_w,ii_c=win_h*w1,ii_d=ii_c+win_w;
                var edges_thresh = ((win_w*win_h) * 0xff * this.edges_density)|0;
                // if too much gradient we also can skip
                //var edges_thresh_high = ((win_w*win_h) * 0xff * 0.3)|0;

                var rects = [];
                for(y = 0; y < ey; y += step_y) {
                    ii_a = y * w1;
                    for(x = 0; x < ex; x += step_x, ii_a += step_x) {

                        mean =    int_sum[ii_a] 
                                - int_sum[ii_a+ii_b]
                                - int_sum[ii_a+ii_c]
                                + int_sum[ii_a+ii_d];

                        // canny prune
                        if(int_canny_sum) {
                            edge_dens = (int_canny_sum[ii_a] 
                                        - int_canny_sum[ii_a+ii_b]
                                        - int_canny_sum[ii_a+ii_c]
                                        + int_canny_sum[ii_a+ii_d]);
                            if(edge_dens < edges_thresh || mean < 20) {
                                x += step_x, ii_a += step_x;
                                continue;
                            }
                        }

                        mean *= inv_area;
                        variance = (int_sqsum[ii_a] 
                                    - int_sqsum[ii_a+ii_b]
                                    - int_sqsum[ii_a+ii_c]
                                    + int_sqsum[ii_a+ii_d]) * inv_area - mean * mean;

                        std = variance > 0. ? Math.sqrt(variance) : 1;

                        stages = classifier.complexClassifiers;
                        sn = stages.length;
                        found =  true;
                        for(i = 0; i < sn; ++i) {
                            stage = stages[i];
                            stage_thresh = stage.threshold;
                            trees = stage.simpleClassifiers;
                            tn = trees.length;
                            stage_sum = 0;
                            for(j = 0; j < tn; ++j) {
                                tree = trees[j];
                                tree_sum = 0;
                                features = tree.features;
                                fn = features.length;
                                if(tree.tilted === 1) {
                                    for(k=0; k < fn; ++k) {
                                        feature = features[k];
                                        fi_a = ~~(x + feature[0] * scale) + ~~(y + feature[1] * scale) * w1;
                                        fw = ~~(feature[2] * scale);
                                        fh = ~~(feature[3] * scale);
                                        fi_b = fw * w1;
                                        fi_c =  fh * w1;

                                        tree_sum += (int_tilted[fi_a]
                                                    - int_tilted[fi_a + fw + fi_b]
                                                    - int_tilted[fi_a - fh + fi_c]
                                                    + int_tilted[fi_a + fw - fh + fi_b + fi_c]) * feature[4];
                                    }
                                } else {
                                    for(k=0; k < fn; ++k) {
                                        feature = features[k];
                                        fi_a = ~~(x + feature[0] * scale) + ~~(y + feature[1] * scale) * w1;
                                        fw = ~~(feature[2] * scale);
                                        fh = ~~(feature[3] * scale);
                                        fi_c = fh * w1;

                                        tree_sum += (int_sum[fi_a] 
                                                    - int_sum[fi_a+fw]
                                                    - int_sum[fi_a+fi_c]
                                                    + int_sum[fi_a+fi_c+fw]) * feature[4];
                                    }
                                }
                                stage_sum += (tree_sum * inv_area < tree.threshold * std) ? tree.left_val : tree.right_val;
                            }
                            if (stage_sum < stage_thresh) {
                                found = false;
                                break;
                            }
                        }
                        
                        if(found) {
                            rects.push({"x" : x,
                                        "y" : y,
                                        "width" : win_w,
                                        "height" : win_h,
                                        "neighbor" : 1,
                                        "confidence" : stage_sum});
                            x += step_x, ii_a += step_x;
                        }
                    }
                }
                return rects;
            },

            detect_multi_scale: function(int_sum, int_sqsum, int_tilted, int_canny_sum, width, height, classifier, scale_factor, scale_min) {
                if (typeof scale_factor === "undefined") { scale_factor = 1.2; }
                if (typeof scale_min === "undefined") { scale_min = 1.0; }
                var win_w = classifier.size[0];
                var win_h = classifier.size[1];
                var rects = [];
                while (scale_min * win_w < width && scale_min * win_h < height) {
                    rects = rects.concat(this.detect_single_scale(int_sum, int_sqsum, int_tilted, int_canny_sum, width, height, scale_min, classifier));
                    scale_min *= scale_factor;
                }
                return rects;
            },

            // OpenCV method to group detected rectangles
            group_rectangles: function(rects, min_neighbors) {
                if (typeof min_neighbors === "undefined") { min_neighbors = 1; }
                var i, j, n = rects.length;
                var node = [];
                for (i = 0; i < n; ++i) {
                    node[i] = {"parent" : -1,
                               "element" : rects[i],
                               "rank" : 0};
                }
                for (i = 0; i < n; ++i) {
                    if (!node[i].element)
                        continue;
                    var root = i;
                    while (node[root].parent != -1)
                        root = node[root].parent;
                    for (j = 0; j < n; ++j) {
                        if( i != j && node[j].element && _group_func(node[i].element, node[j].element)) {
                            var root2 = j;

                            while (node[root2].parent != -1)
                                root2 = node[root2].parent;

                            if(root2 != root) {
                                if(node[root].rank > node[root2].rank)
                                    node[root2].parent = root;
                                else {
                                    node[root].parent = root2;
                                    if (node[root].rank == node[root2].rank)
                                    node[root2].rank++;
                                    root = root2;
                                }

                                /* compress path from node2 to the root: */
                                var temp, node2 = j;
                                while (node[node2].parent != -1) {
                                    temp = node2;
                                    node2 = node[node2].parent;
                                    node[temp].parent = root;
                                }

                                /* compress path from node to the root: */
                                node2 = i;
                                while (node[node2].parent != -1) {
                                    temp = node2;
                                    node2 = node[node2].parent;
                                    node[temp].parent = root;
                                }
                            }
                        }
                    }
                }
                var idx_seq = [];
                var class_idx = 0;
                for(i = 0; i < n; i++) {
                    j = -1;
                    var node1 = i;
                    if(node[node1].element) {
                        while (node[node1].parent != -1)
                            node1 = node[node1].parent;
                        if(node[node1].rank >= 0)
                            node[node1].rank = ~class_idx++;
                        j = ~node[node1].rank;
                    }
                    idx_seq[i] = j;
                }
                
                var comps = [];
                for (i = 0; i < class_idx+1; ++i) {
                    comps[i] = {"neighbors" : 0,
                                "x" : 0,
                                "y" : 0,
                                "width" : 0,
                                "height" : 0,
                                "confidence" : 0};
                }

                // count number of neighbors
                for(i = 0; i < n; ++i) {
                    var r1 = rects[i];
                    var idx = idx_seq[i];

                    if (comps[idx].neighbors == 0)
                        comps[idx].confidence = r1.confidence;

                    ++comps[idx].neighbors;

                    comps[idx].x += r1.x;
                    comps[idx].y += r1.y;
                    comps[idx].width += r1.width;
                    comps[idx].height += r1.height;
                    comps[idx].confidence = Math.max(comps[idx].confidence, r1.confidence);
                }

                var seq2 = [];
                // calculate average bounding box
                for(i = 0; i < class_idx; ++i) {
                    n = comps[i].neighbors;
                    if (n >= min_neighbors)
                        seq2.push({"x" : (comps[i].x * 2 + n) / (2 * n),
                                   "y" : (comps[i].y * 2 + n) / (2 * n),
                                   "width" : (comps[i].width * 2 + n) / (2 * n),
                                   "height" : (comps[i].height * 2 + n) / (2 * n),
                                   "neighbors" : comps[i].neighbors,
                                   "confidence" : comps[i].confidence});
                }

                var result_seq = [];
                n = seq2.length;
                // filter out small face rectangles inside large face rectangles
                for(i = 0; i < n; ++i) {
                    var r1 = seq2[i];
                    var flag = true;
                    for(j = 0; j < n; ++j) {
                        var r2 = seq2[j];
                        var distance = (r2.width * 0.25 + 0.5)|0;

                        if(i != j &&
                           r1.x >= r2.x - distance &&
                           r1.y >= r2.y - distance &&
                           r1.x + r1.width <= r2.x + r2.width + distance &&
                           r1.y + r1.height <= r2.y + r2.height + distance &&
                           (r2.neighbors > Math.max(3, r1.neighbors) || r1.neighbors < 3)) {
                            flag = false;
                            break;
                        }
                    }

                    if(flag)
                        result_seq.push(r1);
                }
                return result_seq;
            }
        };

    })();

    global.haar = haar;

})(jsfeat);
