import sys
import re
import csv

proc = {
    'TPs': int,
    'FPs': int,
    'FNs': int,
    'P': float,
    'R': float,
}

prefix = sys.argv[1] if len(sys.argv) > 1 else ""

res = {}

# truth = re.fullmatch(r"Truth {([0-9]+)}", sys.stdin.readline().strip())
# assert(truth is not None)
# res[prefix + "truth"] = int(truth.group(1))

# called = re.fullmatch(r"Called {([0-9]+)}", sys.stdin.readline().strip())
# assert(called is not None)
# res[prefix + "called"] = int(called.group(1))

assert(sys.stdin.readline().startswith("Truth "))
assert(sys.stdin.readline().startswith("Called"))

reader = csv.DictReader(sys.stdin)
for row in reader:
    for key in row:
        res[prefix + key] = proc[key](row[key])

writer = csv.DictWriter(sys.stdout, fieldnames=res.keys())
writer.writeheader()
writer.writerow(res)
