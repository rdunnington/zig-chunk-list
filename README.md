# zig-chunk-list
`ChunkList` is an address-stable container similar to `std.SegmentedList` that uses a user-specified, consistent chunk size. Allows random access to elements via `at()` in constant time. Since elements are allocated in similarly-sized chunks, allocation overhead is consistent and predictable. The container is not a single contiguous array, but when new chunks are added, no recopies of the data are needed.

Ideal use cases are for large arrays that grow over time, or for address-stable elements if pointers to elements are needed. Since users can specify the chunk size, if adding a lot of elements, it's recommended to choose a size that aligns closely to page boundaries based on the size of the type. 

Typical usage is similar to `std.ArrayList`, except for random access, which must be through `at()`, and inserting/removing elements is only allowed at the end.
