

class Util(object):

    @staticmethod
    def dot_to_none(func):
        """
        Checks if value is '.' and returns None if true.
        Otherwise, func will be called.
        """
        def wrapper(val):
            if val == '.':
                return None
            return func(val)
        return wrapper

    @staticmethod
    def conv_to_number(value):
        """
        Tries to convert a string to a number, silently returning the originally value if it fails.
        """
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        return value

    @staticmethod
    def split_and_convert(conv_func, split_max=-1, extract_single=False):
        """
        Performs a normal split() on a string, with support for converting the values and extraction of single values.

        :param conv_func: Function for converting the values
        :type conv_func: functions
        :param split_max: Maximum number of splits to perform. Default: No limit.
        :type split_max: int
        :param extract_single: If value ends up being a single value, do not return a list. Default: False
        :type extract_single: bool
        """
        def inner(x):
            l = [conv_func(i) for i in x.split(',', split_max)]
            if len(l) == 1 and extract_single:
                l = l[0]
            return l
        return inner
