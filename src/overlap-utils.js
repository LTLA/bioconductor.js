import * as utils from "./utils.js";

export function convertPositionToRank(start, end, { slice = null } = {}) {
    let n = (slice == null ? start.length : slice.length);

    let positions = new Int32Array(n * 2);
    let add = new Uint8Array(n * 2);
    let index = new Int32Array(n * 2);

    {
        let counter = 0;
        let fillIndex = i => {
            let at = counter * 2;
            let next = at + 1;
            positions[at] = start[i];
            positions[next] = end[i];
            add[at] = 1;
            add[next] = 0;
            index[at] = counter;
            index[next] = counter;
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

    let rank2position = [];
    let new_starts = new Int32Array(n);
    let new_ends = new Int32Array(n);

    let last = null;
    for (const i of order) {
        let pos = positions[i];
        let idx = index[i];

        if (pos !== last) {
            rank2position.push(pos);
            last = pos;
        }

        if (add[i]) {
            new_starts[idx] = rank2position.length - 1;
        } else {
            new_ends[idx] = rank2position.length - 1;
        }
    }

    return { rank2position, startRanks: new_starts, endRanks: new_ends };
}

export function buildIntervalTree(start, end, { slice = null } = {}) {
    let { rank2position, startRanks, endRanks } = convertPositionToRank(start, end, { slice });

    // Now, building an nicely balanced interval tree based on the ranks.
    let tree = [ create_node(0, rank2position.length) ];
    if (slice === null) {
        for (var i = 0; i < startRanks.length; i++) {
            recursive_build_tree(startRanks[i], endRanks[i], i, tree, 0);
        }
    } else {
        for (var i = 0; i < startRanks.length; i++) {
            recursive_build_tree(startRanks[i], endRanks[i], slice[i], tree, 0);
        }
    }

    // Running a clean-up operation to convert ranks back to positions.
    let one_past_the_end = (rank2position.length > 0 ? rank2position[rank2position.length - 1] + 1 : 1);
    rank2position.push(one_past_the_end);

    for (const x of tree) {
        x.left_bound = rank2position[x.left_bound];
        x.right_bound = rank2position[x.right_bound];
        x.center = rank2position[x.center];

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
            tree.push(create_node(current.left_bound, current.center));
        }
        recursive_build_tree(start, end, index, tree, current.left_node);

    } else {
        // At some point, every range ends up here. This is because left_bound
        // == center upon successive halving to create new nodes, so every
        // range will eventually overlap a center at its own start position.
        current.overlaps.push(index);
    }
}

export function queryIntervalTree(start, end, tree) {
    let results = [];
    if (start > tree.right_bound) {
        return results;
    }

    if (end < tree.left_bound || (end == tree.left_bound && end > start)) { // Still allow 0-length ranges to fall through for search.
        return results;
    }

    recursive_query_tree(start, end, tree, 0, results); 
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
            if (overlap[0] < end || (overlap[0] == end && start == end)) { // handle zero-length ranges directly on the start position of the center-overlapping range.
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
