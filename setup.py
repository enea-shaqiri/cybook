import setuptools
import Cython.Build as cb
setuptools.setup(name='My Cython Project',
                 ext_modules=cb.cythonize(['src/core/Order.pyx', 'src/core/LimitOrderBook.pyx', 'src/core/binance_connection_test.pyx']))