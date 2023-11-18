
from sortedcontainers import SortedDict
from posix.time cimport clock_gettime, timespec, CLOCK_REALTIME
from Order cimport Order

import numpy as np
#cimport numpy as np
import cython
#DTYPE = np.float64
#ctypedef np.float64_t DTYPE_t

from Cython.Compiler.Options import get_directive_defaults
directive_defaults = get_directive_defaults()

directive_defaults['linetrace'] = True
directive_defaults['binding'] = True

cdef timespec ts



cdef class LimitOrderBook:
    """
    Implementation of the orderbook using lists. Two lists are used for the ask side (prices and volumes) and
    two lists are used for the bid side. All the lists have the same length. The lists containing the sizes are
    initialized to zeros when processing the first order.
    The lists containing the prices are also initialized when the first order is processed.
    Contrary to common implementations, prices_ask contains prices that are below the best ask and prices_bid contains
    prices that are above the best bid. To keep track of the current best bid and best ask, we use two pointers.
    For the orderbook to be consistent, it is necessary that the size is zero for all the indices larger (smaller)
    than the pointer for the ask (bid) side.
    """
    # @TODO Finish docs. Add optional self centering. Add finding empty levels. Add measures of fair price.
    # @TODO Remove benchmarking. Expand book to include adaptors.


    __slots__ = ['prices_ask', 'sizes_ask', 'prices_bid', 'sizes_bid', 'best_price_ask', 'best_price_bid', 'max_size',
                 'index_best_bid', 'index_best_ask', 'index_worst_ask', 'index_worst_bid', 'tick_size', 'decimals']

    def __cinit__(self, int max_size=2000, double tick_size=0.5, int decimals=2):
        # Bid prices decrease with the indices.
        self.prices_bid = []
        self.sizes_bid = []
        self.sizes_ask = []
        # Ask prices increase with the indices.
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
        # The timestamp of the last book update in millisecond. Useful to check book consistency
        self.timestamp_last_updated = 0
        self.time_update = 0
        self.n_update = 0
        self.time_remove = 0
        self.n_remove = 0
        self.temp_bench_1 = 0
        self.n_bench_1 = 0
        self.time_get_new_bids = 0
        self.n_new_bids = 0
        self.time_first_condition = 0
        self.n_first_condition = 0
        self.time_second_condition = 0
        self.n_second_condition = 0
        self.time_third_condition = 0
        self.n_third_condition = 0


    cdef void _first_update(self, Order order):
        min_price = max(order.price - self.tick_size * 0.5 * self.max_size, self.tick_size)
        max_price = order.price + self.tick_size * 0.5 * self.max_size
        prices = list(np.round(np.linspace(min_price, max_price, int((max_price - min_price) / self.tick_size) + 1),
                               self.decimals))
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

    cpdef void process_orders(self, list orders):
        """
        This function process a list of orders and calls process_order to update the book. 
        """
        cdef int j
        for j in range(len(orders)):
            self.process_order(orders[j])

    cpdef process_exchange_messages(self, dict messages):
        pass

    cpdef void process_order(self, Order order):
        cdef double current
        if len(self.prices_bid) == 0 and order.is_bid:
            self._first_update(order)
            return
        elif len(self.prices_ask) == 0 and not order.is_bid:
            self._first_update(order)
            return
        if order.size == 0:
            #clock_gettime(CLOCK_REALTIME, &ts)
            #current = ts.tv_nsec
            self.remove(order)
            #clock_gettime(CLOCK_REALTIME, &ts)
            #self.time_remove += ts.tv_nsec - current
            #self.n_remove += 1
        else:
            #clock_gettime(CLOCK_REALTIME, &ts)
            #current = ts.tv_nsec
            self.update(order)
            #clock_gettime(CLOCK_REALTIME, &ts)
            #self.time_update += ts.tv_nsec - current
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void update(self, Order order):
        cdef int current_index
        cdef double bid_ask_diff
        cdef int n_iterations
        cdef int j
        cdef long current
        if order.is_bid:
            current_index = self.index_best_bid + int((order.price - self.best_price_bid) / self.tick_size)
            self.sizes_bid[current_index] = order.size
            if current_index > self.index_best_bid:
                self.index_best_bid = current_index
                self.best_price_bid = self.get_best_price_bid()
                bid_ask_diff = self.best_price_ask - self.best_price_bid
                self.sizes_bid[current_index] = order.size
                # Fix the book if the ask price is below the bid price
                if bid_ask_diff <= 0:
                    n_iterations = int(- bid_ask_diff / self.tick_size) + 1
                    for j in range(n_iterations):
                        self.sizes_ask[self.index_best_ask] = 0
                        self.index_best_ask += 1
                    self.best_price_ask = self.get_best_price_ask()

        else:
            current_index = self.index_best_ask + int((order.price - self.best_price_ask) / self.tick_size)
            self.sizes_ask[current_index] = order.size
            if current_index < self.index_best_ask:
                self.index_best_ask = current_index
                self.best_price_ask = self.get_best_price_ask()
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
        prices = np.round(np.linspace(min_price, max_price, int((max_price - min_price) / self.tick_size) + 1),
                          self.decimals)
        # @TODO Why is the sortedDict being used?
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
        prices = np.round(np.linspace(min_price, max_price, int((max_price - min_price) / self.tick_size) + 1),
                          self.decimals)
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

    cpdef list get_new_updates_ask(self, list old_updates, double limit):
        """
        This function creates orders from exchange update messages. It will be overriden in child classes
        """
        return []

    cpdef list get_new_updates_bid(self, list old_updates, double limit):
        """
        This function creates orders from exchange update messages. It will be overriden in child classes
        """
        return []

    cpdef list get_ask_orders_from_exchange_message(self, list message):
        """
        This function creates orders from exchange update messages. It will be overriden in child classes
        """
        return []

    cpdef list get_bid_orders_from_exchange_message(self, list message):
        """
        This function creates orders from exchange update messages. It will be overriden in child classes
        """
        return []

cdef class LimitOrderBookBinance(LimitOrderBook):
    def __cinit__(self, int max_size=2000, tick_size=0.1, decimals=1):
        super().__init__(max_size, tick_size, decimals)

    cpdef list get_bid_orders_from_exchange_message(self, list message):
        cdef double current
        clock_gettime(CLOCK_REALTIME, &ts)
        current = ts.tv_nsec
        cdef list orders = []
        cdef double temp
        cdef list old_update
        cdef int j
        for j in range(len(message) - 1, -1, -1):
            old_update = message[j]
            temp = float(old_update[0])
            orders.append([temp, float(old_update[1])])
            self.n_new_bids += 1
            if temp < self.min_price_bid:
                break
        clock_gettime(CLOCK_REALTIME, &ts)
        self.time_get_new_bids += ts.tv_nsec - current
        return orders

    cpdef list get_ask_orders_from_exchange_message(self, list message):
        cdef list new_updates = []
        cdef double temp
        cdef int j
        cdef list old_update
        for j in range(len(message)):
            old_update = message[j]
            temp = float(old_update[0])
            new_updates.append([temp, float(old_update[1])])
            if temp > self.max_price_ask:
                break
        return new_updates

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def process_exchange_messages(self, dict messages):
        # Creating one order. The attributes will be overriden.
        cdef Order placeholder_order = Order(0, 0, 0, 0, 0)
        cdef list bid_messages = messages["data"]["b"]
        cdef list ask_messages = messages["data"]["a"]
        cdef int j
        cdef int bid_length = len(bid_messages)
        # Bid levels are sorted from lowest to highest price
        for j in range(bid_length):
            placeholder_order.price = float(bid_messages[bid_length - j - 1][0])
            # Do not process further if the order is not in the range we care about
            if self.min_price_bid and placeholder_order.price < self.min_price_bid:
                break
            placeholder_order.size = float(bid_messages[bid_length - j - 1][1])
            placeholder_order.is_bid = True
            placeholder_order.timestamp = messages["data"]["E"]
            self.process_order(placeholder_order)
        cdef int ask_length = len(ask_messages)
        # Ask levels are sorted from highest to lowest price
        for j in range(ask_length):
            placeholder_order.price = float(ask_messages[j][0])
            # Do not process further if the order is not in the range we care about
            if self.max_price_ask and placeholder_order.price > self.max_price_ask:
                break
            placeholder_order.size = float(ask_messages[j][1])
            placeholder_order.is_bid = False
            placeholder_order.timestamp = messages["data"]["E"]
            self.process_order(placeholder_order)




