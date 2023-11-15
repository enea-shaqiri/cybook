from posix.time cimport clock_gettime, timespec, CLOCK_REALTIME
cdef timespec ts


cdef class Order:
    __slots__ = ['uid', 'is_bid', 'size', 'price', 'timestamp',
                 'next_item', 'previous_item', 'root']
    def __cinit__(self,
                  int uid,
                  bint is_bid,
                  double size,
                  double price,
                  unsigned long timestamp):
        # Data Values
        self.uid = uid
        self.is_bid = is_bid
        self.price = price
        self.size = size
        clock_gettime(CLOCK_REALTIME, &ts)
        # Conversion to milliseconds
        self.timestamp = timestamp if timestamp != 0 else ts.tv_nsec // 1_000_000

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str((self.uid, self.is_bid, self.price, self.size, self.timestamp))
