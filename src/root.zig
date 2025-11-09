//! By convention, root.zig is the root source file when making a library. If
//! you are making an executable, the convention is to delete this file and
//! start with main.zig instead.
const builtin = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const testing = std.testing;

fn computePow2Mask(power_of_two: usize) usize {
    std.debug.assert(std.math.isPowerOfTwo(power_of_two));

    var mask: usize = 0;
    var bit: usize = power_of_two >> 1;
    while (bit != 0) {
        mask |= bit;
        bit = bit >> 1;
    }
    return mask;
}

pub fn ChunkList(comptime T: type, comptime num_elements_per_chunk: usize) type {
    comptime if (num_elements_per_chunk == 0) {
        @compileError("num_elements_per_chunk must be greater than 0");
    };
    comptime if (std.math.isPowerOfTwo(num_elements_per_chunk) == false) {
        @compileError("num_elements_per_chunk must be a power of two");
    };

    return struct {
        const Self = @This();

        const inner_mask: usize = computePow2Mask(num_elements_per_chunk);
        const chunk_shift: usize = std.math.log2_int(usize, num_elements_per_chunk);

        chunks: [][]T,
        num_elements_in_last_chunk: usize,
        num_full_chunks: usize,

        const InternalIndex = struct {
            inner: usize,
            chunk: usize,

            fn fromAbsolute(index: usize) InternalIndex {
                return .{
                    .inner = index & inner_mask,
                    .chunk = index >> chunk_shift,
                };
            }

            fn next(self: *const InternalIndex) InternalIndex {
                if (self.inner < num_elements_per_chunk - 1) {
                    return .{
                        .inner = self.inner + 1,
                        .chunk = self.chunk,
                    };
                } else {
                    return .{
                        .inner = 0,
                        .chunk = self.chunk + 1,
                    };
                }
            }

            fn lessThan(self: *const InternalIndex, other: *const InternalIndex) bool {
                if (self.chunk < other.chunk) {
                    return true;
                } else if (self.chunk == other.chunk) {
                    return self.inner < other.inner;
                } else {
                    return false;
                }
            }
        };

        const Iterator = struct {
            current: InternalIndex,
            end: InternalIndex,
            list: *Self,

            fn next(self: *Iterator) ?*T {
                var ptr: ?*T = null;
                const next_index = self.current.next();
                if (next_index.lessThan(&self.end)) {
                    ptr = &self.list.chunks[self.current.chunk][self.current.inner];
                    self.current = next_index;
                }
                return ptr;
            }
        };

        fn init() Self {
            return .{
                .chunks = &.{},
                .num_elements_in_last_chunk = 0,
                .num_full_chunks = 0,
            };
        }

        fn deinit(self: *Self, allocator: Allocator) void {
            for (self.chunks) |slice| {
                allocator.free(slice);
            }
            allocator.free(self.chunks);
            self.* = init();
        }

        fn iterate(self: *Self) Iterator {
            const begin: InternalIndex = .{
                .inner = 0,
                .chunk = 0,
            };
            const end: InternalIndex = .{
                .inner = self.num_elements_in_last_chunk,
                .chunk = self.num_full_chunks,
            };
            return .{
                .current = begin,
                .end = end,
                .list = self,
            };
        }

        fn range(self: *Self, start_index: usize, end_index: usize) Iterator {
            const _len = self.len();
            std.debug.assert(start_index < _len);
            std.debug.assert(end_index < _len);

            const begin: InternalIndex = .fromAbsolute(start_index);
            const end: InternalIndex = .fromAbsolute(end_index);

            return .{
                .current = begin,
                .end = end,
                .list = self,
            };
        }

        fn capacity(self: *const Self) usize {
            return self.chunks.len * num_elements_per_chunk;
        }

        fn unusedCapacity(self: *const Self) usize {
            const c = self.capacity();
            const l = self.len();
            return c - l;
        }

        fn len(self: *const Self) usize {
            return self.num_full_chunks * num_elements_per_chunk + self.num_elements_in_last_chunk;
        }

        /// allows indexing this data structure like a contiguous array
        fn at(self: *Self, absolute_index: usize) *T {
            const index = InternalIndex.fromAbsolute(absolute_index);

            std.debug.assert(index.chunk <= self.num_full_chunks);
            std.debug.assert(index.chunk < self.num_full_chunks or index.inner < self.num_elements_in_last_chunk);

            return &self.chunks[index.chunk][index.inner];
        }

        fn resize(self: *Self, new_size: usize, allocator: Allocator) Allocator.Error!void {
            try self.ensureTotalCapacity(new_size, allocator);
            self.num_elements_in_last_chunk = new_size & inner_mask;
            self.num_full_chunks = new_size >> chunk_shift;
        }

        fn ensureTotalCapacity(self: *Self, requested_capacity: usize, allocator: Allocator) Allocator.Error!void {
            const old_capacity = self.capacity();
            if (requested_capacity <= old_capacity) {
                return;
            }

            const num_old_chunks = self.chunks.len;
            const num_new_chunks = (requested_capacity + num_elements_per_chunk - 1) / num_elements_per_chunk;
            if (num_new_chunks <= num_old_chunks) {
                return;
            }

            self.chunks = try allocator.realloc(self.chunks, num_new_chunks);

            const new_chunks = self.chunks[num_old_chunks..num_new_chunks];

            // initialize all the chunk slices so in case one of the allocations fail, it won't be in an undefined state
            for (new_chunks) |*chunk| {
                chunk.* = &.{};
            }

            for (new_chunks) |*chunk| {
                chunk.* = try allocator.alloc(T, num_elements_per_chunk);
            }
        }

        fn ensureUnusedCapacity(self: *Self, unused_capacity: usize, allocator: Allocator) Allocator.Error!void {
            const _len = self.len();
            try self.ensureTotalCapacity(_len + unused_capacity, allocator);
        }

        fn addOne(self: *Self, allocator: Allocator) Allocator.Error!*T {
            try ensureUnusedCapacity(1, allocator);
            return self.addOneAssumeCapacity();
        }

        fn addOneAssumeCapacity(self: *Self) *T {
            std.debug.assert(self.num_full_chunks < self.chunks.len);

            const item: *T = &self.chunks[self.num_full_chunks][self.num_elements_in_last_chunk];
            self.num_elements_in_last_chunk += 1;

            if (self.num_elements_in_last_chunk == num_elements_per_chunk) {
                self.num_elements_in_last_chunk = 0;
                self.num_full_chunks += 1;
            }

            return item;
        }

        fn append(self: *Self, item: T, allocator: Allocator) Allocator.Error!void {
            try self.ensureUnusedCapacity(1, allocator);
            self.appendAssumeCapacity(item);
        }

        fn appendAssumeCapacity(self: *Self, item: T) void {
            const new_item: *T = self.addOneAssumeCapacity();
            new_item.* = item;
        }

        fn appendNTimes(self: *Self, value: T, n: usize, allocator: Allocator) Allocator.Error!void {
            try self.ensureUnusedCapacity(n, allocator);
            self.appendNTimesAssumeCapacity(value, n);
        }

        fn appendNTimesAssumeCapacity(self: *Self, value: T, n: usize) void {
            var num_inserted: usize = 0;
            var dest_begin = self.num_elements_in_last_chunk;

            while (num_inserted < n) {
                const num_to_copy = @min(num_elements_per_chunk - dest_begin, n - num_inserted);

                const dest: []T = self.chunks[self.num_full_chunks][dest_begin .. dest_begin + num_to_copy];
                @memset(dest, value);

                std.debug.assert(dest_begin + num_to_copy <= num_elements_per_chunk);
                dest_begin = if (dest_begin + num_to_copy == num_elements_per_chunk) 0 else dest_begin + num_to_copy;
                self.num_full_chunks += if (dest_begin == 0) 1 else 0;
                num_inserted += num_to_copy;
            }

            self.num_elements_in_last_chunk = dest_begin;
        }

        fn appendSlice(self: *Self, items: []const T, allocator: Allocator) Allocator.Error!void {
            try self.ensureUnusedCapacity(items.len, allocator);
            self.appendSliceAssumeCapacity(items);
        }

        fn appendSliceAssumeCapacity(self: *Self, items: []const T) void {
            // std.debug.print("[before] self.num_full_chunks: {}, self.num_elements_in_last_chunk: {}\n", .{ self.num_full_chunks, self.num_elements_in_last_chunk });

            var src_begin: usize = 0;
            var dest_begin = self.num_elements_in_last_chunk;

            while (src_begin < items.len) {
                const num_to_copy = @min(num_elements_per_chunk - dest_begin, items.len - src_begin);

                const dest: []T = self.chunks[self.num_full_chunks][dest_begin .. dest_begin + num_to_copy];
                const src: []const T = items[src_begin .. src_begin + num_to_copy];
                // std.debug.print("\tnum to copy : {}, items.len: {}, src.len: {}, dest.len: {}, src_begin: {}\n", .{ num_to_copy, items.len, src.len, dest.len, src_begin });
                @memcpy(dest, src);

                std.debug.assert(dest_begin + src.len <= num_elements_per_chunk);
                dest_begin = if (dest_begin + src.len == num_elements_per_chunk) 0 else dest_begin + src.len;
                self.num_full_chunks += if (dest_begin == 0) 1 else 0;
                src_begin += src.len;
            }

            self.num_elements_in_last_chunk = dest_begin;

            // std.debug.print("[after] self.num_full_chunks: {}, self.num_elements_in_last_chunk: {}\n", .{ self.num_full_chunks, self.num_elements_in_last_chunk });
        }

        fn clearRetainingCapacity(self: *Self) void {
            self.num_elements_in_last_chunk = 0;
            self.num_full_chunks = 0;
        }

        fn concat(self: *const Self, allocator: Allocator) Allocator.Error![]T {
            var all: []T = try allocator.alloc(T, self.len());
            var begin: usize = 0;

            const last_chunk: usize = self.num_full_chunks + (if (self.num_elements_in_last_chunk > 0) @as(usize, 1) else 0);
            std.debug.assert(last_chunk <= self.chunks.len);

            for (self.chunks[0..last_chunk], 0..) |chunk, chunk_index| {
                const num_elements_in_chunk = if (chunk_index == self.num_full_chunks) self.num_elements_in_last_chunk else num_elements_per_chunk;
                const src: []const T = chunk[0..num_elements_in_chunk];
                const dest: []T = all[begin .. begin + src.len];
                @memcpy(dest, src);
                begin += src.len;
            }

            return all;
        }
    };
}

test "pow2 mask" {
    try std.testing.expectEqual(0b01, computePow2Mask(2));
    try std.testing.expectEqual(0b011, computePow2Mask(4));
    try std.testing.expectEqual(0b0111, computePow2Mask(8));
    try std.testing.expectEqual(0b01111, computePow2Mask(16));
    try std.testing.expectEqual(0b011111, computePow2Mask(32));
}

const TestList = ChunkList(usize, 4);

test "list init/deinit" {
    var a = TestList.init();

    try std.testing.expectEqual(0, a.capacity());
    try std.testing.expectEqual(0, a.len());

    a.deinit(std.testing.allocator);

    try std.testing.expectEqual(0, a.capacity());
    try std.testing.expectEqual(0, a.len());
}

test "add/append" {
    const allocator = std.testing.allocator;

    var a = TestList.init();
    defer a.deinit(allocator);

    // add 5 items to get 2 chunks allocated
    for (0..5) |i| {
        try a.append(i, allocator);

        try std.testing.expectEqual((i / 4 + 1) * 4, a.capacity());
        try std.testing.expectEqual(i + 1, a.len());
    }

    try std.testing.expect(a.capacity() > 0);
    try std.testing.expectEqual(5, a.len());

    for (0..5) |i| {
        try std.testing.expectEqual(i, a.at(i).*);
    }
}

test "concat" {
    const allocator = std.testing.allocator;

    var a = TestList.init();
    defer a.deinit(allocator);

    for (0..10) |i| {
        try a.append(i, allocator);
    }

    const all = try a.concat(allocator);
    defer allocator.free(all);

    const expected: []const usize = &.{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    try std.testing.expectEqualSlices(usize, expected, all);
}

test "appendSlice" {
    const allocator = std.testing.allocator;

    var a = TestList.init();
    defer a.deinit(allocator);

    const slice: []const usize = &.{ 0, 1, 2, 3, 4 };

    for (0..3) |_| {
        try a.appendSlice(slice, allocator);
    }

    try std.testing.expectEqual(16, a.capacity());
    try std.testing.expectEqual(15, a.len());

    try a.append(0, allocator);
    try std.testing.expectEqual(16, a.len());
    try a.appendSlice(&.{ 0, 1 }, allocator);
    try std.testing.expectEqual(20, a.capacity());
    try std.testing.expectEqual(18, a.len());

    try a.appendSlice(&.{ 0, 1 }, allocator);
    try std.testing.expectEqual(20, a.capacity());
    try std.testing.expectEqual(20, a.len());

    try a.appendSlice(&.{ 0, 1, 2, 4 }, allocator);
    try std.testing.expectEqual(24, a.capacity());
    try std.testing.expectEqual(24, a.len());

    const all: []const usize = &.{ 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 0, 1, 0, 1, 0, 1, 2, 4 };
    const concat = try a.concat(allocator);
    defer allocator.free(concat);
    try std.testing.expectEqualSlices(usize, all, concat);
}

test "appendNTimes" {
    const allocator = std.testing.allocator;

    var a = TestList.init();
    defer a.deinit(allocator);

    try a.appendNTimes(0, 1, allocator);
    try std.testing.expectEqual(4, a.capacity());
    try std.testing.expectEqual(1, a.len());

    try a.appendNTimes(1, 3, allocator);
    try std.testing.expectEqual(4, a.capacity());
    try std.testing.expectEqual(4, a.len());

    try a.appendNTimes(2, 2, allocator);
    try a.appendNTimes(3, 2, allocator);
    try std.testing.expectEqual(8, a.capacity());
    try std.testing.expectEqual(8, a.len());

    try a.appendNTimes(4, 9, allocator);
    try std.testing.expectEqual(20, a.capacity());
    try std.testing.expectEqual(17, a.len());

    try a.appendNTimes(5, 3, allocator);
    try std.testing.expectEqual(20, a.capacity());
    try std.testing.expectEqual(20, a.len());

    try a.appendNTimes(6, 4, allocator);
    try std.testing.expectEqual(24, a.capacity());
    try std.testing.expectEqual(24, a.len());

    const all: []const usize = &.{ 0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6 };
    const concat = try a.concat(allocator);
    defer allocator.free(concat);
    try std.testing.expectEqualSlices(usize, all, concat);
}

test "resize" {
    const allocator = std.testing.allocator;

    var a = TestList.init();
    defer a.deinit(allocator);

    for (0..5) |i| {
        try a.append(i, allocator);
    }

    try a.resize(0, allocator);
    try std.testing.expect(a.capacity() > 0);
    try std.testing.expectEqual(0, a.len());

    try a.resize(12, allocator);
    try std.testing.expectEqual(12, a.capacity());
    try std.testing.expectEqual(12, a.len());

    try a.resize(13, allocator);
    try std.testing.expectEqual(16, a.capacity());
    try std.testing.expectEqual(13, a.len());

    try a.resize(16, allocator);
    try std.testing.expectEqual(16, a.capacity());
    try std.testing.expectEqual(16, a.len());
}

test "iterate" {
    const allocator = std.testing.allocator;

    var a = TestList.init();
    defer a.deinit(allocator);

    for (0..10) |i| {
        try a.append(i, allocator);
    }

    var counter: usize = 0;
    var iter = a.iterate();
    while (iter.next()) |ptr| {
        try std.testing.expectEqual(counter, ptr.*);
        counter += 1;
    }
}

test "range" {
    const allocator = std.testing.allocator;

    var a = TestList.init();
    defer a.deinit(allocator);

    for (0..100) |i| {
        try a.append(i, allocator);
    }

    var counter: usize = 57;
    var iter = a.range(57, 87);
    while (iter.next()) |ptr| {
        try std.testing.expectEqual(counter, ptr.*);
        counter += 1;
    }

    iter = a.range(57, 57);
    try std.testing.expectEqual(null, iter.next());

    iter = a.range(57, 56);
    try std.testing.expectEqual(null, iter.next());
}

test "clear" {
    const allocator = std.testing.allocator;

    var a = TestList.init();
    defer a.deinit(allocator);

    for (0..16) |i| {
        try a.append(i, allocator);
    }

    try std.testing.expectEqual(16, a.capacity());
    try std.testing.expectEqual(16, a.len());

    a.clearRetainingCapacity();
    try std.testing.expectEqual(16, a.capacity());
    try std.testing.expectEqual(0, a.len());
}
