from Crypto.Util.number import bytes_to_long, getPrime
import os
import sys
from secret import FLAG

p = getPrime(128)
q = getPrime(128)
n = p * q
phi = (p - 1) * (q - 1)
e = 65537 # public key
d = pow(e, -1, phi) # private key

m = bytes_to_long(FLAG + os.urandom(32 - len(FLAG)))
assert m < n
enc = pow(m, e, n)

print('n=', n)
print('e=', e)
print('encrypted msg=', enc)

while True:
    inp = input('What is the flag? ').strip().encode()
    if inp == FLAG:
        print('Good!')
        exit(0)
    else:
        print('No no no!')
