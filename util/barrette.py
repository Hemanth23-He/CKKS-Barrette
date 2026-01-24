class BarrettReducer:
    def __init__(self, modulus):
        self.modulus = modulus
        self.k = modulus.bit_length() * 2
        self.mu = (1 << self.k) // modulus

    def reduce(self, x):
        q = ((x * self.mu) >> self.k)
        r = x - q * self.modulus
        if r >= self.modulus:
            r -= self.modulus
        if r < 0:
            r += self.modulus
        return r


    def inverse(self, x):
    a, b = x % self.modulus, self.modulus
    u0, u1 = 1, 0
    while b:
        q = a // b
        a, b = b, a - q * b
        u0, u1 = u1, u0 - q * u1
    return u0% self.modulus
