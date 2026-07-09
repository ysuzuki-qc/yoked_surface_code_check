
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

flist = glob.glob("./data/gap_stat_*.txt")
for fname in flist:
    name = fname.replace("\\","/").split("/")[-1].replace(".txt", "").replace("gap_stat_", "")
    value = {}
    with open(fname, "r") as f:
        for line in f:
            db, count = line.strip().split(" ")
            value[int(db)] = int(count)
    x = []
    y0 = []
    y1 = []
    y2 = []
    for db in range(0, max(value.keys())+1):
        success = value.get(db, 0)
        fail = value.get(-db, 0)
        if success + fail == 0:
            logical_error_rate = 0.0
        else:
            logical_error_rate = fail / (success + fail)
        x.append(db)
        y0.append(logical_error_rate)
        def calc_logical_error_rate(db):
            relative_error_rate = 10**(-db/10)
            return relative_error_rate/(1+relative_error_rate)
        y1.append(calc_logical_error_rate(db))
        y2.append(calc_logical_error_rate(db*0.9))
    plt.figure(figsize=(10, 10))
    plt.plot(x,y0,label="Logical Error Rate")
    plt.plot(x,y1,label="10^(-dB/10)")
    plt.plot(x,y2,label="10^(-0.9*dB/10)")
    plt.xlim(0,80)
    plt.ylim(1e-8,1e0)
    plt.xlabel("Complementary gap [dB]")
    plt.ylabel("Logical Error Rate")
    plt.legend()
    plt.yscale("log")
    plt.savefig(f"./figure/calibration_{name}.png")
