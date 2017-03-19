from new_data_type import DoubleList
from random import random
import timeit

v = [random() for _ in range(10000)]

vd = DoubleList(v)
