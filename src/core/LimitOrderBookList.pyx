from sortedcontainers import SortedDict
from posix.time cimport clock_gettime, timespec, CLOCK_REALTIME
import cython
cimport cython
from libc.stdio cimport printf

import numpy as np
#cimport numpy as np

#DTYPE = np.float64
#ctypedef np.float64_t DTYPE_t

from Cython.Compiler.Options import get_directive_defaults
directive_defaults = get_directive_defaults()

directive_defaults['linetrace'] = True
directive_defaults['binding'] = True

cdef timespec ts

cdef class Order:

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
        # Fix this
        self.timestamp = timestamp if timestamp != 0 else clock_gettime(CLOCK_REALTIME, &ts)


    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str((self.uid, self.is_bid, self.price, self.size, self.timestamp))






cdef class LimitOrderBookList:


    def __cinit__(self, int max_size=2000, tick_size=0.5, decimals=2):
        self.prices_bid = []
        self.sizes_bid = []
        self.sizes_ask = []
        self.prices_ask = []
        self.best_price_ask
        self.best_price_bid
        self.max_size = max_size
        self.index_best_ask = 0
        self.index_worst_ask = max_size
        self.index_best_bid = max_size
        self.index_worst_bid = 0
        self.tick_size = tick_size
        self.decimals = decimals
        self.time_update = 0
        self.n_update = 0
        self.time_remove = 0
        self.n_remove = 0
        self.temp_bench_1 = 0
        self.n_bench_1 = 0
        self.time_get_new_bids = 0
        self.n_new_bids = 0

    cdef void _first_update(self, Order order):
        min_price = max(order.price - self.tick_size * 0.5 * self.max_size, self.tick_size)
        max_price = order.price + self.tick_size * 0.5 * self.max_size
        prices = list(np.round(np.linspace(min_price, max_price, int((max_price - min_price) / self.tick_size) + 1), self.decimals))
        if order.is_bid:
            self.prices_bid = prices
            self.sizes_bid = [0] * len(self.prices_bid)
            self.min_price_bid = self.prices_bid[0]
            self.index_best_bid = int(self.max_size / 2)
            self.best_price_bid = self.get_best_price_bid()
            self.sizes_bid[self.index_best_bid] = order.size
        else:
            self.prices_ask = prices
            self.sizes_ask = [0] * len(self.prices_ask)
            self.max_price_ask = self.prices_ask[-1]
            self.index_best_ask = int(self.max_size / 2)
            self.best_price_ask = self.get_best_price_ask()
            self.sizes_ask[self.index_best_ask] = order.size

    cdef void process_orders(self, list orders):
        cdef int j
        for j in range(len(orders)):
            self.process(orders[j])

    cdef void process(self, Order order):
        cdef double current
        if len(self.prices_bid) == 0 and order.is_bid:
                self._first_update(order)
                return
        elif len(self.prices_ask) == 0 and not order.is_bid:
                self._first_update(order)
                return
        if order.size == 0:
            clock_gettime(CLOCK_REALTIME, &ts)
            current = ts.tv_nsec
            self.remove(order)
            clock_gettime(CLOCK_REALTIME, &ts)
            self.time_remove += ts.tv_nsec - current
            self.n_remove += 1
        else:
            clock_gettime(CLOCK_REALTIME, &ts)
            current = ts.tv_nsec
            self.update(order)
            clock_gettime(CLOCK_REALTIME, &ts)
            self.time_update = ts.tv_nsec - current
            self.n_update += 1

    cpdef list get_cum_volume_ask(self, int levels):
        cdef int cumsum = self.sizes_ask[self.index_best_ask]
        cdef list cumsums = [cumsum]
        cdef int j
        for j in range(1, levels):
            if self.index_best_ask + j > len(self.sizes_ask):
                break
            cumsums.append(cumsums[-1] + self.sizes_ask[self.index_best_ask + j])
        return cumsums

    cpdef list get_cum_volume_bid(self, int levels):
        cdef int cumsum = self.sizes_bid[self.index_best_bid]
        cdef list cumsums = [cumsum]
        cdef int j
        for j in range(1, levels):
            if self.index_best_bid - j < 0:
                break
            cumsums.append(cumsums[-1] + self.sizes_bid[self.index_best_bid - j])
        return cumsums

    cdef void update(self, Order order):
        cdef int current_index
        cdef double bid_ask_diff
        cdef int n_iterations
        cdef int j

        if order.is_bid:
            if self.min_price_bid > order.price or self.prices_bid[-1] < order.price:
                return
            current_index = self.index_best_bid + int((order.price - self.best_price_bid) / self.tick_size)
            self.sizes_bid[current_index] = order.size
            if current_index > self.index_best_bid:
                self.index_best_bid = current_index
                self.best_price_bid = self.get_best_price_bid()
                if self.best_price_ask == 0.0:
                    return
                bid_ask_diff = self.best_price_ask - self.best_price_bid
                # Fix the book if the ask price is below the bid price
                if bid_ask_diff <= 0:
                    n_iterations = int(- bid_ask_diff / self.tick_size) + 1
                    for j in range(n_iterations):
                        self.sizes_ask[self.index_best_ask] = 0
                        self.index_best_ask += 1
                    self.best_price_ask = self.get_best_price_ask()
        else:
            if self.max_price_ask < order.price or self.prices_ask[0] > order.price:
                return
            current_index = self.index_best_ask + int((order.price - self.best_price_ask) / self.tick_size)
            self.sizes_ask[current_index] = order.size
            if current_index < self.index_best_ask:
                self.index_best_ask = current_index
                self.best_price_ask = self.get_best_price_ask()
                if self.best_price_bid == 0.0:
                    return
                bid_ask_diff = self.best_price_ask - self.best_price_bid
                # Fix the book if the ask price is below the bid price
                if bid_ask_diff <= 0:
                    n_iterations = int(- bid_ask_diff / self.tick_size) + 1
                    for j in range(n_iterations):
                        self.sizes_bid[self.index_best_bid] = 0
                        self.index_best_bid -= 1
                    self.best_price_bid = self.get_best_price_bid()



    cdef void remove(self, Order order):
        cdef int j
        if order.is_bid:
            if self.min_price_bid > order.price: return
            if order.price > self.best_price_bid: return
            current_index = self.index_best_bid + int((order.price - self.best_price_bid) / self.tick_size)
            try:
                self.sizes_bid[current_index] = 0
            except IndexError:
                print(order.price, self.best_price_bid, self.index_best_bid, order.size)
                #print(current_index)
            if order.price == self.best_price_bid:
                # find new index_best_bid
                for j in range(self.index_best_bid - 1, -1, -1):
                    if self.sizes_bid[j] > 0:
                        self.index_best_bid = j
                        self.best_price_bid = self.get_best_price_bid()
                        return
        else:
            if self.max_price_ask < order.price: return
            if order.price < self.best_price_ask: return
            current_index = self.index_best_ask + int((order.price - self.best_price_ask) / self.tick_size)
            self.sizes_ask[current_index] = 0
            if order.price == self.best_price_ask:
                # find new index_best_ask
                for j in range(self.index_best_ask + 1, len(self.sizes_ask)):
                    if self.sizes_ask[j] > 0:
                        self.index_best_ask = j
                        self.best_price_ask = self.get_best_price_ask()

                        return

    cdef void balance_bids(self):
        cdef double best_price
        cdef double best_size
        cdef int j
        best_price, best_size = self.get_best_bid()
        min_price = max(best_price - self.tick_size * 0.5 * self.max_size, self.tick_size)
        max_price = best_price + self.tick_size * 0.5 * self.max_size
        prices = np.round(np.linspace(min_price, max_price, int((max_price - min_price) / self.tick_size) + 1), self.decimals)
        temp_sorted_dict = SortedDict()
        for j in range(len(prices)):
            volume = self.bids.get(prices[j], None)
            if volume is None:
                temp_sorted_dict[prices[j]] = 0
            else:
                temp_sorted_dict[prices[j]] = volume
        self.bids = temp_sorted_dict
        self.bid_indices = list(self.bids.keys())
        self.min_price_bid = self.bid_indices[0]

    cdef void balance_asks(self):
        cdef double best_price
        cdef double best_size
        best_price, best_size = self.get_best_ask()
        min_price = max(best_price - self.tick_size * 0.5 * self.max_size, self.tick_size)
        max_price = best_price + self.tick_size * 0.5 * self.max_size
        prices = np.round(np.linspace(min_price, max_price, int((max_price - min_price) / self.tick_size) + 1), self.decimals)
        temp_sorted_dict = SortedDict()
        for j in range(len(prices)):
            volume = self.asks.get(prices[j], 0)
            temp_sorted_dict[prices[j]] = volume
        self.asks = temp_sorted_dict
        self.ask_indices = list(self.bids.keys())
        self.max_price_ask = self.ask_indices[-1]

    cdef (double, double) get_best_ask(self):
        return self.get_best_price_ask(), self.sizes_ask[self.index_best_ask]

    cdef double get_best_price_ask(self):
        return self.prices_ask[self.index_best_ask]

    cdef double get_ask_price_at_index(self, int index):
        return self.prices_ask[index]

    cdef (double, double) get_ask_at_index(self, int index):
        return self.get_ask_price_at_index(index), self.sizes_ask[index]

    cdef double get_best_price_bid(self):
        clock_gettime(CLOCK_REALTIME, &ts)
        current = ts.tv_nsec
        cdef double temp = self.prices_bid[self.index_best_bid]
        clock_gettime(CLOCK_REALTIME, &ts)
        self.temp_bench_1 += ts.tv_nsec - current
        self.n_bench_1 += 1
        return temp

    cdef (double, double) get_best_bid(self):
        return self.get_best_price_bid(), self.sizes_bid[self.index_best_bid]

    cdef double get_bid_price_at_index(self, int index):
        return self.prices_bid[index]
    cdef (double, double) get_bid_at_index(self, int index):
        return self.get_bid_price_at_index(index), self.sizes_bid[index]


    cpdef list get_new_updates_bid(self, list old_updates, double limit):
        """

        """

        #cdef np.ndarray[DTYPE, ndim=2] new_updates
        #cdef int [:, ::1] new_updates
        cdef double current
        clock_gettime(CLOCK_REALTIME, &ts)
        current = ts.tv_nsec
        cdef list new_updates = []
        cdef double temp
        cdef list old_update
        cdef int j
        for j in range(len(old_updates) - 1, -1, -1):
            old_update = old_updates[j]
            #
            temp = float(old_update[0])
            new_updates.append([temp, float(old_update[1])])
            self.n_new_bids += 1
            if temp < limit:
                break
        clock_gettime(CLOCK_REALTIME, &ts)
        self.time_get_new_bids += ts.tv_nsec - current

        return new_updates

    cpdef list get_new_updates_ask(self, list old_updates, double limit):
        cdef list new_updates = []
        cdef double temp
        cdef list old_update
        for j in range(len(old_updates)):
            old_update = old_updates[j]
            temp = float(old_update[0])
            new_updates.append([temp, float(old_update[1])])
            if temp < limit:
                break
        return new_updates


