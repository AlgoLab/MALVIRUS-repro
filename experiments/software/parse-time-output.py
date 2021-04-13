import csv
import sys

interesting = {
#    "Command being timed": "command",
    "User time (seconds)": "user",
    "System time (seconds)": "system",
    "Percent of CPU this job got": "perc_cpu",
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": "wall",
    "Maximum resident set size (kbytes)": "max_rss",
    "Exit status": "exit",
}

def parse_wall(x):
    x = x.split(":")
    res = 0
    for y in x:
        res = 60 * res + float(y)
    return res

proc = {
#   "command": str,
    "user": float,
    "system": float,
    "perc_cpu": lambda x: int(x.rstrip('%')),
    "wall": parse_wall,
    "max_rss": int,
    "exit": int,
}

prefix = sys.argv[1] if len(sys.argv) > 1 else ""

res = {}
for line in sys.stdin:
    line = line.strip()
    (key, value) = line.split(': ', 1)
    if key in interesting.keys():
        res[prefix + interesting[key]] = proc[interesting[key]](value)

writer = csv.DictWriter(sys.stdout, fieldnames=res.keys())
writer.writeheader()
writer.writerow(res)
