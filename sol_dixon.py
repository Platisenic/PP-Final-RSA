from pwn import *
from Crypto.Util.number import long_to_bytes
import subprocess
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--method', type=str, default='dixon_omp', help='dixon, dixon_omp, dixon_mpi')
    parser.add_argument('--threads', type=int, default=12, help='number of threads')
    parser.add_argument('--host', type=str, default='libos446.duckdns.org', help='attacked server')
    parser.add_argument('--port', type=int, default=10001, help='attacked server')
    args = parser.parse_args()

    r = remote(args.host, args.port)
    n = r.recvline().decode()
    print(n, end='')
    n = int(n[n.find('= ')+2:])
    e = r.recvline().decode()
    print(e, end='')
    e = int(e[e.find('= ')+2:])
    enc = r.recvline().decode()
    print(enc, end='')
    enc = int(enc[enc.find('= ')+2:])

    subprocess.run(['make', '-C', 'dixon'])

    io = process([f'./dixon/{args.method}', str(n), str(args.threads)])
    pq = io.recvline().decode().split()
    p = int(pq[0])
    q = int(pq[1])
    assert p * q == n, 'Dixon Computation Error'

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

