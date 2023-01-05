from Crypto.Util.number import bytes_to_long, getPrime
import os
from secret import FLAG

p = getPrime(46)
q = getPrime(46)
n = p * q
phi = (p - 1) * (q - 1)

e = 65537 # public key
d = pow(e, -1, phi) # private key

m_len = 10
assert len(FLAG) <= m_len
m = bytes_to_long(FLAG + os.urandom(m_len - len(FLAG)))
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
