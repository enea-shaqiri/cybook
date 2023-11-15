cdef class Order:
    cdef public int uid
    cdef public bint is_bid
    cdef public double price
    cdef public double size
    cdef public unsigned long timestamp