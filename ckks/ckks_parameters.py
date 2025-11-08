"""A module to keep track of parameters for the CKKS scheme."""
import math

class CKKSParameters:
    """An instance of parameters for the CKKS scheme.
    Attributes:
        poly_degree (int): Degree d of polynomial that determines the quotient ring R.
        ciph_modulus (int): Coefficient modulus of ciphertexts.
        big_modulus (int): Large modulus used for bootstrapping.
        scaling_factor (float): Scaling factor to multiply by.
        hamming_weight (int): Hamming weight parameter for sampling secret key.
        num_taylor_iterations (int): Number of iterations to perform for Taylor series in bootstrapping.
        prime_size (int): Minimum number of bits in primes for RNS representation.
    """
    def __init__(self, poly_degree, ciph_modulus, big_modulus, scaling_factor, taylor_iterations=6, prime_size=59):
        """Inits Parameters with the given parameters.
        Args:
            poly_degree (int): Degree d of polynomial of ring R.
            ciph_modulus (int): Coefficient modulus of ciphertexts.
            big_modulus (int): Large modulus used for bootstrapping.
            scaling_factor (float): Scaling factor to multiply by.
            taylor_iterations (int): Number of iterations to perform for Taylor series in bootstrapping.
            prime_size (int): Minimum number of bits in primes for RNS representation. (Not used if CRT not needed)
        """
        self.poly_degree = poly_degree
        self.ciph_modulus = ciph_modulus
        self.big_modulus = big_modulus
        self.scaling_factor = scaling_factor
        self.num_taylor_iterations = taylor_iterations
        self.hamming_weight = poly_degree // 4
        self.prime_size = prime_size

    def print_parameters(self):
        """Prints parameters."""
        print("Encryption parameters")
        print("\t Polynomial degree: %d" % (self.poly_degree))
        print("\t Ciphertext modulus size: %d bits" % (int(math.log(self.ciph_modulus, 2))))
        print("\t Big ciphertext modulus size: %d bits" % (int(math.log(self.big_modulus, 2))))
        print("\t Scaling factor size: %d bits" % (int(math.log(self.scaling_factor, 2))))
        print("\t Number of Taylor iterations: %d" % (self.num_taylor_iterations))
        print("\t RNS: No")
