"""A module to multiply polynomials using the Fast Fourier Transform (FFT), Number Theoretic Transform (NTT),
and Fermat Theoretic Transform (FTT). See https://rijndael.ece.vt.edu/schaum/pdf/papers/2013hostb.pdf.
"""
"""A module to multiply polynomials using the Fast Fourier Transform (FFT), 
Number Theoretic Transform (NTT), and Fermat Theoretic Transform (FTT). 
See https://rijndael.ece.vt.edu/schaum/pdf/papers/2013hostb.pdf.
"""

from math import log, pi, cos, sin
import util.number_theory as nbtheory
from util.bit_operations import bit_reverse_vec, reverse_bits
from util.barrette import BarrettReducer  # NEW IMPORT (if using Montgomery: from util.montgomery import MontgomeryReducer)

class NTTContext:
    """
    An instance of Number/Fermat Theoretic Transform parameters.
    ...
    """
    def __init__(self, poly_degree, coeff_modulus, root_of_unity=None):
        assert (poly_degree & (poly_degree - 1)) == 0, \
            "Polynomial degree must be a power of 2. d = " + str(poly_degree) + " is not."
        self.coeff_modulus = coeff_modulus
        self.degree = poly_degree
        self.reducer = BarrettReducer(coeff_modulus)  # NEW: instantiate your reducer

        if not root_of_unity:
            root_of_unity = nbtheory.root_of_unity(order=2 * poly_degree, modulus=coeff_modulus)
        self.precompute_ntt(root_of_unity)

    def precompute_ntt(self, root_of_unity):
        # Find powers of root of unity.
        self.roots_of_unity = [1] * self.degree
        for i in range(1, self.degree):
            self.roots_of_unity[i] = self.reducer.reduce(self.roots_of_unity[i - 1] * root_of_unity)  # CHANGED

        # Find powers of inverse root of unity.
        root_of_unity_inv = self.reducer.inverse(root_of_unity)  # CHANGED
        self.roots_of_unity_inv = [1] * self.degree
        for i in range(1, self.degree):
            self.roots_of_unity_inv[i] = self.reducer.reduce(self.roots_of_unity_inv[i - 1] * root_of_unity_inv)  # CHANGED

        # Compute precomputed array of reversed bits for iterated NTT.
        self.reversed_bits = [0] * self.degree
        width = int(log(self.degree, 2))
        for i in range(self.degree):
            self.reversed_bits[i] = reverse_bits(i, width) % self.degree

    def ntt(self, coeffs, rou):
        num_coeffs = len(coeffs)
        assert len(rou) == num_coeffs, \
            "Length of the roots of unity is too small. Length is " + str(len(rou))
        result = bit_reverse_vec(coeffs)
        log_num_coeffs = int(log(num_coeffs, 2))
        for logm in range(1, log_num_coeffs + 1):
            for j in range(0, num_coeffs, (1 << logm)):
                for i in range(1 << (logm - 1)):
                    index_even = j + i
                    index_odd = j + i + (1 << (logm - 1))
                    rou_idx = (i << (1 + log_num_coeffs - logm))
                    omega_factor = self.reducer.reduce(rou[rou_idx] * result[index_odd])  # CHANGED
                    butterfly_plus = self.reducer.reduce(result[index_even] + omega_factor)  # CHANGED
                    butterfly_minus = self.reducer.reduce(result[index_even] - omega_factor)  # CHANGED
                    result[index_even] = butterfly_plus
                    result[index_odd] = butterfly_minus
        return result

    def ftt_fwd(self, coeffs):
        num_coeffs = len(coeffs)
        assert num_coeffs == self.degree, "ftt_fwd: input length does not match context degree"
        ftt_input = [self.reducer.reduce(int(coeffs[i]) * self.roots_of_unity[i]) for i in range(num_coeffs)]  # CHANGED
        return self.ntt(coeffs=ftt_input, rou=self.roots_of_unity)

    def ftt_inv(self, coeffs):
        num_coeffs = len(coeffs)
        assert num_coeffs == self.degree, "ntt_inv: input length does not match context degree"
        to_scale_down = self.ntt(coeffs=coeffs, rou=self.roots_of_unity_inv)
        poly_degree_inv = self.reducer.inverse(self.degree)  # CHANGED
        result = [self.reducer.reduce(int(to_scale_down[i]) * self.roots_of_unity_inv[i] * poly_degree_inv)
                  for i in range(num_coeffs)]  # CHANGED
        return result

# FFTContext and related complex/float code are unchanged

class FFTContext:
    ...
    # No change needed for complex-valued FFT; leave as-is

# (remainder of FFT code unchanged)


class FFTContext:
    """An instance of Fast Fourier Transform (FFT) parameters.

    The FFTContext keeps track of the length of the vector and precomputations
    to perform FFT.

    Attributes:
        fft_length (int): Length of the FFT vector. This must be twice the polynomial degree.
        roots_of_unity (list): The ith member of the list is w^i, where w
            is a root of unity.
        rot_group (list): Used for EMB only. Value at index i is 5i (mod fft_length)
            for 0 <= i < fft_length / 4.
        reversed_bits (list): The ith member of the list is the bits of i
            reversed, used in the iterative implementation of FFT.
    """
    def __init__(self, fft_length):
        """Inits FFTContext with a length for the FFT vector.

        Args:
            fft_length (int): Length of the FFT vector.
        """
        self.fft_length = fft_length
        self.precompute_fft()

    def precompute_fft(self):
        """Performs precomputations for the FFT.

        Precomputes all powers of roots of unity for the FFT and powers of inverse
        roots of unity for the inverse FFT.
        """
        self.roots_of_unity = [0] * self.fft_length
        self.roots_of_unity_inv = [0] * self.fft_length
        for i in range(self.fft_length):
            angle = 2 * pi * i / self.fft_length
            self.roots_of_unity[i] = complex(cos(angle), sin(angle))
            self.roots_of_unity_inv[i] = complex(cos(-angle), sin(-angle))

        # Compute precomputed array of reversed bits for iterated FFT.
        num_slots = self.fft_length // 4
        self.reversed_bits = [0] * num_slots
        width = int(log(num_slots, 2))
        for i in range(num_slots):
            self.reversed_bits[i] = reverse_bits(i, width) % num_slots

        # Compute rotation group for EMB with powers of 5.
        self.rot_group = [1] * num_slots
        for i in range(1, num_slots):
            self.rot_group[i] = (5 * self.rot_group[i - 1]) % self.fft_length

    def fft(self, coeffs, rou):
        """Runs FFT on the given coefficients.

        Runs iterated FFT with the given coefficients and roots of unity. See
        paper for pseudocode.

        Args:
            coeffs (list): List of coefficients to transform. Must be the
                length of the polynomial degree.
            rou (list): Powers of roots of unity to be used for transformation.
                For inverse NTT, this is the powers of the inverse root of unity.

        Returns:
            List of transformed coefficients.
        """
        num_coeffs = len(coeffs)
        assert len(rou) >= num_coeffs, \
            "Length of the roots of unity is too small. Length is " + str(len(rou))

        result = bit_reverse_vec(coeffs)

        log_num_coeffs = int(log(num_coeffs, 2))

        for logm in range(1, log_num_coeffs + 1):
            for j in range(0, num_coeffs, (1 << logm)):
                for i in range(1 << (logm - 1)):
                    index_even = j + i
                    index_odd = j + i + (1 << (logm - 1))

                    rou_idx = (i * self.fft_length) >> logm
                    omega_factor = rou[rou_idx] * result[index_odd]

                    butterfly_plus = result[index_even] + omega_factor
                    butterfly_minus = result[index_even] - omega_factor

                    result[index_even] = butterfly_plus
                    result[index_odd] = butterfly_minus

        return result

    def fft_fwd(self, coeffs):
        """Runs forward FFT on the given values.

        Runs forward FFT with the given values and parameters in the context.

        Args:
            coeffs (list): List of complex numbers to transform.

        Returns:
            List of transformed coefficients.
        """
        return self.fft(coeffs, rou=self.roots_of_unity)

    def fft_inv(self, coeffs):
        """Runs inverse FFT on the given values.

        Runs inverse FFT with the given values and parameters in the context.

        Args:
            coeffs (list): List of complex numbers to transform.

        Returns:
            List of transformed coefficients.
        """
        num_coeffs = len(coeffs)
        result = self.fft(coeffs, rou=self.roots_of_unity_inv)

        for i in range(num_coeffs):
            result[i] /= num_coeffs

        return result

    def check_embedding_input(self, values):
        """Checks that the length of the input vector to embedding is the correct size.

        Throws an error if the length of the input vector to embedding is not 1/4 the size
        of the FFT vector.

        Args:
            values (list): Input vector of complex numbers.
        """
        assert len(values) <= self.fft_length / 4, "Input vector must have length at most " \
            + str(self.fft_length / 4) + " < " + str(len(values)) + " = len(values)"

    def embedding(self, coeffs):
        """Computes a variant of the canonical embedding on the given coefficients.

        Computes the canonical embedding which consists of evaluating a given polynomial at roots of unity
        that are indexed 1 (mod 4), w, w^5, w^9, ...
        The evaluations are returned in the order: w, w^5, w^(5^2), ...

        Args:
            coeffs (list): List of complex numbers to transform.

        Returns:
            List of transformed coefficients.
        """
        self.check_embedding_input(coeffs)
        num_coeffs = len(coeffs)
        result = bit_reverse_vec(coeffs)
        log_num_coeffs = int(log(num_coeffs, 2))

        for logm in range(1, log_num_coeffs + 1):
            idx_mod = 1 << (logm + 2)
            gap = self.fft_length // idx_mod
            for j in range(0, num_coeffs, (1 << logm)):
                for i in range(1 << (logm - 1)):
                    index_even = j + i
                    index_odd = j + i + (1 << (logm - 1))

                    rou_idx = (self.rot_group[i] % idx_mod) * gap
                    omega_factor = self.roots_of_unity[rou_idx] * result[index_odd]

                    butterfly_plus = result[index_even] + omega_factor
                    butterfly_minus = result[index_even] - omega_factor

                    result[index_even] = butterfly_plus
                    result[index_odd] = butterfly_minus

        return result

    def embedding_inv(self, coeffs):
        """Computes the inverse variant of the canonical embedding.

        Args:
            values (list): List of complex numbers to transform.

        Returns:
            List of transformed coefficients.
        """
        self.check_embedding_input(coeffs)
        num_coeffs = len(coeffs)
        result = coeffs.copy()
        log_num_coeffs = int(log(num_coeffs, 2))

        for logm in range(log_num_coeffs, 0, -1):
            idx_mod = 1 << (logm + 2)
            gap = self.fft_length // idx_mod
            for j in range(0, num_coeffs, 1 << logm):
                for i in range(1 << (logm - 1)):
                    index_even = j + i
                    index_odd = j + i + (1 << (logm - 1))

                    rou_idx = (self.rot_group[i] % idx_mod) * gap

                    butterfly_plus = result[index_even] + result[index_odd]
                    butterfly_minus = result[index_even] - result[index_odd]
                    butterfly_minus *= self.roots_of_unity_inv[rou_idx]

                    result[index_even] = butterfly_plus
                    result[index_odd] = butterfly_minus

        to_scale_down = bit_reverse_vec(result)

        for i in range(num_coeffs):
            to_scale_down[i] /= num_coeffs

        return to_scale_down
