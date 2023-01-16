
export function buildIntervalTree(start, end, { slice = null } = {}) {
    // Converting all start/end positions into ranks.
    let n = (slice == null ? this._start.length : slice.length);
    let positions = new Int32Array(n * 2);

    {
        let counter = 0;
        let fillIndex = i => {
            let at = counter * 2;
            let next = at + 1;
            positions[at] = start[i];
            positions[next] = end[i];
            counter++;
        };

        if (slice === null) {
            for (var i = 0; i < n; i++) {
                fillIndex(i);                                
            }
        } else {
            for (const i of slice) {
                fillIndex(i);
            }
        }
    }

    let order = utils.createSequence(positions.length);
    order.sort((i, j) => positions[i] - positions[j]);

    let ranks2position = [];
    let position2ranks = {};
    let last = -Number.POSITIVE_INFINITY;
    for (const i of order) {
        let pos = positions[i];
        if (pos != last) {
            position2ranks[pos] = ranks2position.length;
            ranks2position.push(pos);
            last = pos;
        }
    }
    
    // Now, building an nicely balanced interval tree based on the ranks.
    let tree = [ create_node(0, ranks2position.length - 1) ];

    {
        let build_tree = i => {
            let r1 = position2ranks[start[i]];
            let r2 = position2ranks[end[i]];
            recursive_build_tree(r1, r2, i, tree, 0);
        };

        if (slice === null) {
            for (var i = 0; i < n; i++) {
                build_tree(i);
            }
        } else {
            for (const i of slice) {
                build_tree(i);
            }
        }
    }

    // Running a clean-up operation to convert ranks back to positions.
    for (const x of tree) {
        x.left_bound = ranks2position[x.left_bound];
        x.right_bound = ranks2position[x.right_bound];
        x.center = ranks2position[x.center];

        // Also sorting ranges by increasing start and DECREASING end positions.
        let start_overlaps_sorted = x.overlaps.slice().sort((a, b) => start[a] - start[b]);
        let end_overlaps_sorted = x.overlaps.sort((a, b) => end[b] - end[a]) // reversed order - deliberate!
        x.overlaps = {
            start: start_overlaps_sorted.map(i => [start[i], i]),
            end: end_overlaps_sorted.map(i => [end[i], i])
        };
    }

    return tree;
}

function create_node(left_bound, right_bound) {
    return { 
        left_bound: left_bound,
        right_bound: right_bound,
        center: left_bound + Math.floor((right_bound - left_bound) / 2),
        left_node: null,
        right_node: null,
        overlaps: []
    };
}

function recursive_build_tree(start, end, index, tree, node) {
    let current = tree[node];

    if (start > current.center) {
        if (current.right_node === null) {
            current.right_node = tree.length;
            tree.push(create_node(current.center, current.right_bound));
        }
        recursive_build_tree(start, end, index, tree, current.right_node);

    } else if (end < current.center || (end == current.center && end > start)) { // Let 0-length ranges fall through to the next clause if they lie exactly on the center.
        if (current.left_node === null) {
            current.left_node = tree.length;
            tree.push(create_node(current.left_bound, center));
        }
        recursive_build_tree(start, end, index, tree, current.left_node);

    } else {
        // At some point, everyone ends up here, because the new nodes are created
        // such that left_bound == center upon successive halving; so every range
        // will end up overlapping a center defined at its start position.
        current.overlaps.push(index);
    }
}

export function queryIntervalTree(start, end, tree) {
    let results = [];
    recursive_query_tree(start, end, tree, 0, results):    
    return results;
}

function recursive_query_tree(start, end, tree, node, results) {
    let current = tree[node];

    if (start > current.center) {
        for (const overlap of current.overlaps.end) {
            if (overlap[0] > start) {
                results.push(overlap[1]);
            } else {
                break;
            }
        }
        if (current.right_node !== null) {
            recursive_query_tree(start, end, tree, current.right_node, results);
        }

    } else if (end < current.center || (end == current.center && end > start)) { // Again, let zero-length ranges fall through if they lie directly on the center.
        for (const overlap of current.overlaps.start) {
            if (overlap[0] < end) {
                results.push(overlap[1]);
            } else {
                break;
            }
        }
        if (current.left_node !== null) {
            recursive_query_tree(start, end, tree, current.left_node, results);
        }

    } else {
        for (const overlap of current.overlaps.start) {
            results.push(overlap[1]);
        }

        if (end > current.center) {
            if (current.right_node !== null) {
                recursive_query_tree(start, end, tree, current.right_node, results);
            }
        }
        if (start < current.center) {
            if (current.left_node !== null) {
                recursive_query_tree(start, end, tree, current.left_node, results);
            }
        }
    }
}
