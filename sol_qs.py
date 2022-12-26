from pwn import *
from Crypto.Util.number import long_to_bytes
import subprocess

r = remote('libos446.duckdns.org', 10001)
n = r.recvline().decode()
print(n, end='')
n = int(n[n.find('= ')+2:])
e = r.recvline().decode()
print(e, end='')
e = int(e[e.find('= ')+2:])
enc = r.recvline().decode()
print(enc, end='')
enc = int(enc[enc.find('= ')+2:])

subprocess.run(['make', '-C', 'qs'])

io = process(['./qs/quadratic', f'{n}', '12'], stderr=open('/dev/null', 'w+b'))
pq = io.recvline().decode().split()
p = int(pq[0])
q = int(pq[1])
assert p * q == n, "Dixon Computation Error"

phi = (p - 1) * (q - 1)
d = pow(e, -1, phi)
m = pow(enc, d, n)
m = long_to_bytes(m).decode('ascii', errors='ignore')
m = m[:m.find('}')+1]
msg = r.recvuntil(b'?').decode()
print(msg)
print(f'The flag is:{m}')
r.sendline(m.encode('ascii'))
msg = r.recvline().decode()
print(msg)

