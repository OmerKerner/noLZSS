import time, sys, statistics as stats
import noLZSS

def run_once(data: bytes):
    t0 = time.perf_counter()
    noLZSS.factorize(data)
    return time.perf_counter() - t0

def main():
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'rb') as f:
            data = f.read()
    else:
        data = (b"abracadabra$") * 1000
    if not data.endswith(b"$"):
        data += b"$"
    for _ in range(3):
        run_once(data)
    times = [run_once(data) for _ in range(5)]
    mean = stats.mean(times)
    print(f"Bytes: {len(data)} Factors/s: {len(data)/mean:,.0f} Time/iter: {mean*1e3:.2f} ms")

if __name__ == '__main__':
    main()
