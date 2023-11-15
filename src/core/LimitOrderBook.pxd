from Order cimport Order

cdef class LimitOrderBook:

    # Ask prices increase with the indices.
    cdef public prices_ask
    cdef public sizes_ask
    # Bid prices decrease with the indices.
    cdef public prices_bid
    cdef public sizes_bid
    # Number of prices to keep.
    cdef public int max_size
    cdef public int index_best_ask
    cdef public int index_worst_ask
    cdef public int index_best_bid
    cdef public int index_worst_bid
    cdef public double best_price_ask
    cdef public double best_price_bid
    cdef public double max_price_ask
    cdef public double min_price_bid
    cdef public double tick_size
    cdef public int decimals
    cdef public double time_update
    cdef public double time_remove
    cdef public int n_update
    cdef public int n_remove
    cdef public int temp_bench_1
    cdef public int n_bench_1
    cdef public double time_get_new_bids
    cdef public int n_new_bids

    cdef void _first_update(self, Order order)
    cpdef void process_orders(self, list orders)
    cpdef void process(self, Order order)
    cpdef process_exchange_messages(self, dict messages)
    cpdef list get_cum_volume_ask(self, int levels)
    cpdef list get_cum_volume_bid(self, int levels)
    cdef void update(self, Order order)
    cdef void remove(self, Order order)
    cdef void balance_bids(self)
    cdef void balance_asks(self)
    cdef (double, double) get_best_ask(self)
    cdef double get_best_price_ask(self)
    cdef double get_ask_price_at_index(self, int index)
    cdef (double, double) get_ask_at_index(self, int index)
    cdef double get_best_price_bid(self)
    cdef (double, double) get_best_bid(self)
    cdef double get_bid_price_at_index(self, int index)
    cdef (double, double) get_bid_at_index(self, int index)
    cpdef list get_new_updates_bid(self, list old_updates, double limit)
    cpdef list get_new_updates_ask(self, list old_updates, double limit)