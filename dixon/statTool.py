#!/usr/bin/env python3

import sys
import subprocess

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} N", file=sys.stderr)
        print(f"N: Big number of RSA", file=sys.stderr)
        sys.exit(1)
    N = sys.argv[1]
    threads = [1, 2, 4, 6, 8, 12]
    pattern = "time : "
    repeat = 3
    subprocess.run("make")
    for thread in threads:
        totalTime = 0
        for _ in range(repeat):
            out = subprocess.run(["./dixon_omp", N, f"{thread}"], capture_output=True)
            out = out.stdout.decode()
            out = out[out.find(pattern)+len(pattern):out.rfind("\n")]
            totalTime += float(out)
        avgTime = totalTime / repeat
        print(f'N: {N} thread: {thread:2} time: {avgTime:.5f}')

