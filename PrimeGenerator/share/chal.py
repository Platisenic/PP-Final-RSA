from Crypto.Util.number import getPrime

print('Welcome to libos Random Prime Generator')

while True:
    inp = input('\nHow many bits do you need? ').strip()
    try:
        inp = int(inp)//2
    except ValueError:
        print('invalid number :(')
    else:
        if inp > 1000:
            print('too many bits :(')
            continue
        if inp <= 1:
            print('invalid number :(')
            continue

        p = getPrime(inp)
        q = getPrime(inp)
        pq = p * q
        print('p =', p)
        print('q =', q)
        print('pq=', pq)

