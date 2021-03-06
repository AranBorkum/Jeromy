""" Module to hold classes for various test statistics that can be used
for fitting.

.. note:: All forms of chi-squared are as defined in:

  * S. Baker & R. D. Cousins, Nucl. Inst. and Meth. in Phys. Res. 211,
    437-442 (1984)
"""
import numpy
import abc


class TestStatistic(object):
    """ Base class for the calculation of a test statistic.

    The calculation of any test statistic is based on one spectrum
    containing observed events and one containing expected events. It
    is assumed that the observed events form the "data" spectrum and
    the expected events form the spectrum predicted by the model.

    Args:
      name (string): Name of test statistic.
      per_bin (bool): If True (default) the statistic in each bin is
        returned as a :class:`numpy.array`. If False one value for the
        statistic is returned for the entire array.

    Attributes:
      _name (string): Name of test statistic.
      _per_bin (bool): If True the statistic in each bin is returned as
        a :class:`numpy.array`. If False one value for the statistic is
        returned for the entire array.
    """
    __metaclass__ = abc.ABCMeta  # Only required for python 2

    def __init__(self, name, per_bin=False):
        self._name = name
        self._per_bin = per_bin

    def get_name(self):
        """
        Returns:
          string: Name of test statistic, stored in :attr:`_name`
        """
        return self._name

    def compute_statistic(self, observed, expected):
        """ Compute the value of the test statistic.

        Args:
          observed (:class:`numpy.ndarray`): 1D array containing the
            observed data points.
          expected (:class:`numpy.ndarray`): 1D array containing the
            expected values predicted by the model.

        Returns:
          (float or :class:`numpy.array`): Computed value(s) of test
            statistic. A float is returned if :attr:`_per_bin` is
            False. Otherwise an array of computed test statistic
            values, for each bin, is returned.
        """
        if len(observed.shape) != 1:
            raise TypeError("Incompatible shape %s for observed array, "
                            "expected 1-D array" % str(observed.shape))
        if len(expected.shape) != 1:
            raise TypeError("Incompatible shape %s for expected array, "
                            "expected 1-D array" % str(expected.shape))
        if len(observed) != len(expected):
            raise ValueError(
                "Number of bins mismatch, for observed with %d bins, "
                "and expected with %d bins" % (len(observed), len(expected)))
        if not self._per_bin:
            return self._compute(observed, expected)
        else:
            return self._get_stats(observed, expected)

    # Method adapted from http://codereview.stackexchange.com/a/47115
    @abc.abstractmethod
    def _compute(self, observed, expected):
        """ Calculates the test statistic.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          float: Calculated test statistic.
        """
        return None

    @abc.abstractmethod
    def _get_stats(self, observed, expected):
        """ Gets the test statistic for each bin.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Raises:
          ValueError: If arrays are different lengths.

        Returns:
          :class:`numpy.array`: Calculated chi squared for each bin.
        """
        return None

    @abc.abstractmethod
    def get_penalty_term(self, current_value, prior, sigma):
        """ Calculates a penalty term value, for a given fit parameter,
        for this test statistic.

        Args:
          current_value (float): current value of a given fit parameter
          prior (float): Prior value of a given fit parameter
          sigma (float): Sigma value of a given fit parameter

        Returns:
          float: Value of the penalty term
        """
        return None


class BakerCousinsChi(TestStatistic):
    """ Test statistic class for calculating the Baker-Cousins
    chi-squared (poisson likelihood) test statistic.

    Args:
      name (string): Name of test statistic.
      per_bin (bool): If True (default) the statistic in each bin is
        returned as a :class:`numpy.array`. If False one value for the
        statistic is returned for the entire array.

    Attributes:
      _name (string): Name of test statistic.
      _per_bin (bool): If True the statistic in each bin is returned as
        a :class:`numpy.array`. If False one value for the statistic is
        returned for the entire array.
    """
    def __init__(self, per_bin=False):
        super(BakerCousinsChi, self).__init__("baker_cousins", per_bin)

    @classmethod
    def _compute(self, observed, expected):
        """ Calculates the chi-squared.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          float: Calculated Baker-Cousins chi squared
        """
        # Convert to arrays, if observed and expected are floats
        if isinstance(observed, float):
            observed = numpy.array([observed])
        if isinstance(expected, float):
            expected = numpy.array([expected])
        observed = observed.astype('float')
        expected = expected.astype('float')
        epsilon = 1e-34  # In the limit of zero
        total = 0
        for exp, obs in zip(expected, observed):
            if exp < epsilon:
                exp = epsilon
            if obs < epsilon:
                bin_value = exp
            else:
                bin_value = exp - obs + obs * numpy.log(obs / exp)
            total += bin_value
        return 2. * total

    @classmethod
    def _get_stats(self, observed, expected):
        """ Gets chi squared for each bin.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Raises:
          ValueError: If arrays are different lengths.

        Returns:
          :class:`numpy.array`: Calculated chi squared for each bin
        """
        not_per_bin = self._compute(observed, expected)
        observed = observed.astype('float')
        expected = expected.astype('float')
        epsilon = 1e-34  # In the limit of zero
        stats = []
        for exp, obs in zip(expected, observed):
            if exp < epsilon:
                exp = epsilon
            if obs < epsilon:
                bin_value = exp
            else:
                bin_value = exp - obs + obs * numpy.log(obs / exp)
            stats.append(bin_value)
        stats = 2. * numpy.array(stats)
        if not numpy.allclose(numpy.sum(stats), not_per_bin):
            raise ValueError("Sum of chi squared array and value returned by "
                             "_compute method for same observed and expected "
                             "do not match!")
        return stats

    @classmethod
    def get_penalty_term(self, current_value, prior, sigma):
        """ Calculates a penalty term value, for a given fit parameter,
        for the BakerCousinsChi test statistic.

        Args:
          current_value (float): current value of a given fit parameter
          prior (float): Prior value of a given fit parameter
          sigma (float): Sigma value of a given fit parameter

        Returns:
          float: Value of the penalty term
        """
        return ((current_value - prior)/float(sigma)) ** 2


class BakerCousinsChiUp(TestStatistic):
    """ Test statistic class for calculating the Baker-Cousins
    chi-squared (poisson likelihood) test statistic. This test statistic only
    counts upward fluctuations.

    Args:
      name (string): Name of test statistic.
      per_bin (bool): If True (default) the statistic in each bin is
        returned as a :class:`numpy.array`. If False one value for the
        statistic is returned for the entire array.

    Attributes:
      _name (string): Name of test statistic.
      _per_bin (bool): If True the statistic in each bin is returned as
        a :class:`numpy.array`. If False one value for the statistic is
        returned for the entire array.
    """
    def __init__(self, per_bin=False):
        super(BakerCousinsChiUp, self).__init__("baker_cousins_up", per_bin)

    @classmethod
    def _compute(self, observed, expected):
        """ Calculates the chi-squared.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          float: Calculated Baker-Cousins chi squared
        """
        # Convert to arrays, if observed and expected are floats
        if isinstance(observed, float):
            observed = numpy.array([observed])
        if isinstance(expected, float):
            expected = numpy.array([expected])
        observed = observed.astype('float')
        expected = expected.astype('float')
        epsilon = 1e-34  # In the limit of zero
        total = 0
        for exp, obs in zip(expected, observed):
            if obs > exp:
                continue
            if exp < epsilon:
                exp = epsilon
            if obs < epsilon:
                bin_value = exp
            else:
                bin_value = exp - obs + obs * numpy.log(obs / exp)
            total += bin_value
        return 2. * total

    @classmethod
    def _get_stats(self, observed, expected):
        """ Gets chi squared for each bin.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Raises:
          ValueError: If arrays are different lengths.

        Returns:
          :class:`numpy.array`: Calculated chi squared for each bin
        """
        not_per_bin = self._compute(observed, expected)
        observed = observed.astype('float')
        expected = expected.astype('float')
        epsilon = 1e-34  # In the limit of zero
        stats = []
        for exp, obs in zip(expected, observed):
            if exp < epsilon:
                exp = epsilon
            if obs < epsilon:
                bin_value = exp
            else:
                bin_value = exp - obs + obs * numpy.log(obs / exp)
            stats.append(bin_value)
        stats = 2. * numpy.array(stats)
        if not numpy.allclose(numpy.sum(stats), not_per_bin):
            raise ValueError("Sum of chi squared array and value returned by "
                             "_compute method for same observed and expected "
                             "do not match!")
        return stats

    @classmethod
    def get_penalty_term(self, current_value, prior, sigma):
        """ Calculates a penalty term value, for a given fit parameter,
        for the BakerCousinsChi test statistic.

        Args:
          current_value (float): current value of a given fit parameter
          prior (float): Prior value of a given fit parameter
          sigma (float): Sigma value of a given fit parameter

        Returns:
          float: Value of the penalty term
        """
        return ((current_value - prior)/float(sigma)) ** 2


class BakerCousinsLL(TestStatistic):
    """ Test statistic class for calculating the Baker-Cousins log likelihood
      ratio test statistic.

    Args:
      name (string): Name of test statistic.
      per_bin (bool): If True (default) the statistic in each bin is
        returned as a :class:`numpy.array`. If False one value for the
        statistic is returned for the entire array.

    Attributes:
      _name (string): Name of test statistic.
      _per_bin (bool): If True the statistic in each bin is returned as
        a :class:`numpy.array`. If False one value for the statistic is
        returned for the entire array.
    """
    def __init__(self, per_bin=False):
        super(BakerCousinsLL, self).__init__("baker_cousins_ll", per_bin)

    @classmethod
    def _compute(self, observed, expected):
        """ Calculates the log likelihood.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          float: Calculated log-likelihood value
        """
        # Convert to arrays, if observed and expected are floats
        if isinstance(observed, float):
            observed = numpy.array([observed])
        if isinstance(expected, float):
            expected = numpy.array([expected])
        epsilon = 1e-34  # In the limit of zero
        total = 0
        for i in range(len(observed)):
            if expected[i] < epsilon:
                expected[i] = epsilon
            if observed[i] < epsilon:
                bin_value = expected[i]
            else:
                bin_value = expected[i] - observed[i] + observed[i] * numpy.log(observed[i] / expected[i])
            total += bin_value
        return total

    @classmethod
    def _get_stats(self, observed, expected):
        """ Gets chi squared for each bin.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          :class:`numpy.array`: Calculated log-likelihood for each bin
        """
        epsilon = 1e-34  # In the limit of zero
        stats = []
        for i in range(len(observed)):
            if expected[i] < epsilon:
                expected[i] = epsilon
            if observed[i] < epsilon:
                bin_value = expected[i]
            else:
                bin_value = expected[i] - observed[i] + observed[i] *\
                    numpy.log(observed[i] / expected[i])
            stats.append(bin_value)
        return numpy.array(stats)

    @classmethod
    def get_penalty_term(self, current_value, prior, sigma):
        """ Calculates a penalty term value, for a given fit parameter,
        for the BakerCousinsLL test statistic.

        Args:
          current_value (float): current value of a given fit parameter
          prior (float): Prior value of a given fit parameter
          sigma (float): Sigma value of a given fit parameter

        Returns:
          float: Value of the penalty term
        """
        return 0.5 * ((current_value - prior)/sigma) ** 2


class Neyman(TestStatistic):
    """ Test statistic class for calculating the Neyman chi-squared test
    statistic.

    Args:
      name (string): Name of test statistic.
      per_bin (bool): If True (default) the statistic in each bin is
        returned as a :class:`numpy.array`. If False one value for the
        statistic is returned for the entire array.

    Attributes:
      _name (string): Name of test statistic.
      _per_bin (bool): If True the statistic in each bin is returned as
        a :class:`numpy.array`. If False one value for the statistic is
        returned for the entire array.
    """
    def __init__(self, per_bin=False):
        super(Neyman, self).__init__("neyman", per_bin)

    @classmethod
    def _compute(self, observed, expected):
        """ Calculates chi squared.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          float: Calculated Neyman's chi squared
        """
        # Convert to arrays, if observed and expected are floats
        if isinstance(observed, float):
            observed = numpy.array([observed])
        if isinstance(expected, float):
            expected = numpy.array([expected])
        # Chosen due to backgrounds with low rates in ROI
        epsilon = 1e-34  # In the limit of zero
        total = 0
        for i in range(len(observed)):
            if observed[i] < epsilon:
                observed[i] = epsilon
            if expected[i] < epsilon:
                bin_value = observed[i]
            else:
                bin_value = (expected[i] - observed[i])**2 / observed[i]
            total += bin_value
        return total

    @classmethod
    def _get_stats(self, observed, expected):
        """ Gets chi squared for each bin.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          :class:`numpy.array`: Calculated chi squared for each bin
        """
        # Chosen due to backgrounds with low rates in ROI
        epsilon = 1e-34  # In the limit of zero
        stats = []
        for i in range(len(observed)):
            if observed[i] < epsilon:
                observed[i] = epsilon
            if expected[i] < epsilon:
                bin_value = observed[i]
            else:
                bin_value = (expected[i] - observed[i])**2 / observed[i]
            stats.append(bin_value)
        return numpy.array(stats)

    @classmethod
    def get_penalty_term(self, current_value, prior, sigma):
        """ Calculates a penalty term value, for a given fit parameter,
        for the Neyman chi squared test statistic.

        Args:
          current_value (float): current value of a given fit parameter
          prior (float): Prior value of a given fit parameter
          sigma (float): Sigma value of a given fit parameter

        Returns:
          float: Value of the penalty term
        """
        return ((current_value - prior)/sigma) ** 2


class Pearson(TestStatistic):
    """ Test statistic class for calculating the Pearson chi-squared test
    statistic.

    Args:
      name (string): Name of test statistic.
      per_bin (bool): If True (default) the statistic in each bin is
        returned as a :class:`numpy.array`. If False one value for the
        statistic is returned for the entire array.

    Attributes:
      _name (string): Name of test statistic.
      _per_bin (bool): If True the statistic in each bin is returned as
        a :class:`numpy.array`. If False one value for the statistic is
        returned for the entire array.
    """
    def __init__(self, per_bin=False):
        super(Pearson, self).__init__("pearson", per_bin)

    @classmethod
    def _compute(self, observed, expected):
        """ Calculates chi squared.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Raises:
          ValueError: If arrays are different lengths.

        Returns:
          float: Calculated Pearson's chi squared
        """
        # Convert to arrays, if observed and expected are floats
        if isinstance(observed, float):
            observed = numpy.array([observed])
        if isinstance(expected, float):
            expected = numpy.array([expected])
        # Chosen due to backgrounds with low rates in ROI
        epsilon = 1e-34  # Limit of zero
        total = 0
        for i in range(len(observed)):
            if expected[i] < epsilon:
                expected[i] = epsilon
            if observed[i] < epsilon:
                bin_value = expected[i]
            else:
                bin_value = (observed[i] - expected[i])**2 / expected[i]
            total += bin_value
        return total

    @classmethod
    def _get_stats(self, observed, expected):
        """ Gets chi squared for each bin.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          :class:`numpy.array`: Calculated chi squared for each bin
        """
        # Chosen due to backgrounds with low rates in ROI
        epsilon = 1e-34  # In the limit of zero
        stats = []
        for i in range(len(observed)):
            if expected[i] < epsilon:
                expected[i] = epsilon
            if observed[i] < epsilon:
                bin_value = expected[i]
            else:
                bin_value = (observed[i] - expected[i])**2 / expected[i]
            stats.append(bin_value)
        return numpy.array(stats)

    @classmethod
    def get_penalty_term(self, current_value, prior, sigma):
        """ Calculates a penalty term value, for a given fit parameter,
        for the Pearson chi squared test statistic.

        Args:
          current_value (float): current value of a given fit parameter
          prior (float): Prior value of a given fit parameter
          sigma (float): Sigma value of a given fit parameter

        Returns:
          float: Value of the penalty term
        """
        return ((current_value - prior)/sigma) ** 2


class ExtendedLL(TestStatistic):
    """ Test statistic class for calculating the Extended Log Likelihood test
    statistic.

    Args:
      name (string): Name of test statistic.
      per_bin (bool): If True (default) the statistic in each bin is
        returned as a :class:`numpy.array`. If False one value for the
        statistic is returned for the entire array.

    Attributes:
      _name (string): Name of test statistic.
      _per_bin (bool): If True the statistic in each bin is returned as
        a :class:`numpy.array`. If False one value for the statistic is
        returned for the entire array.
    """
    def __init__(self, per_bin=False):
        super(ExtendedLL, self).__init__("extended_ll", per_bin)

    @classmethod
    def _compute(self, observed, expected):
        """ Calculates the Extended Log Likelihood value.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          float: Calculated Extended Log Likelihood value.
        """
        # Convert to arrays, if observed and expected are floats
        if isinstance(observed, float):
            observed = numpy.array([observed])
        if isinstance(expected, float):
            expected = numpy.array([expected])
        # Chosen due to backgrounds with low rates in ROI
        epsilon = 1e-34  # In the limit of zero
        total = 0
        for exp, obs in zip(expected, observed):
            if exp < epsilon:
                exp = epsilon
            bin_value = obs * numpy.log(exp)
            total += bin_value
        return expected.sum() - total

    @classmethod
    def _get_stats(self, observed, expected):
        """ Gets the Extended Log Likelihood value for each bin.

        Args:
          observed (:class:`numpy.array` or float): Number of observed
            events
          expected (:class:`numpy.array` or float): Number of expected
            events

        Returns:
          :class:`numpy.array`: Calculated Extended Log Likelihood value
          for each bin
        """
        # Chosen due to backgrounds with low rates in ROI
        epsilon = 1e-34  # In the limit of zero
        stats = []
        for exp, obs in zip(expected, observed):
            if exp < epsilon:
                exp = epsilon
            bin_value = obs * numpy.log(exp)
            stats.append(exp - bin_value)
        return numpy.array(stats)

    @classmethod
    def get_penalty_term(self, current_value, prior, sigma):
        """ Calculates a penalty term value, for a given fit parameter,
        for the Extended Log Likelihood test statistic.

        Args:
          current_value (float): current value of a given fit parameter
          prior (float): Prior value of a given fit parameter
          sigma (float): Sigma value of a given fit parameter

        Returns:
          float: Value of the penalty term
        """
        return ((current_value - prior)/sigma) ** 2
