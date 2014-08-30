#!/usr/bin/env python2
import Tkinter as tk
import tkFont
import tkFileDialog
import tkMessageBox
import pcap
import dpkt
import numpy as np
import matplotlib
import time

matplotlib.use("TKAgg")


import matplotlib.pyplot as plt
from classes import StaticClasses


# Make a dict of statistical observations for a list/array
def mmm(ds):
    if len(ds) < 1:
        return {"len": 0}
    return {"len": len(ds),
            "min": min(ds),
            "mean": np.mean(ds),
            "std": np.std(ds),
            "max": max(ds)}


# Windowing functions for FFT transformation
dict_win = {
    "none": lambda x: np.tile(1.0, x),
    "hanning": np.hanning,
    "hamming": np.hamming,
    "blackman": lambda x: np.blackman(x, 14),
    "kaiser": lambda x: np.kaiser(x, 14),
    "bartlett": np.bartlett}


# Make an FFT based on input data and an optional window method
def make_fft(sli, window="none"):
    if window and dict_win[window]:
        win = dict_win[window](len(sli))
        filtered = np.zeros_like(sli)
        for x in range(len(sli)):
            filtered[x] = win[x] * sli[x]
    else:
        filtered = sli

    f = np.fft.fft(filtered)
    data = abs(f[:len(f)/2])

    if min(filtered) == max(filtered):
        # No data(signal) to analyze, zero as output
        return np.zeros_like(data)

    mx = float(max(data))
    return map(lambda x: x/mx, data)

    return data


# Make an autocorrelation using fftconvolve
def make_acor(sli):
    from scipy.signal import fftconvolve

    # TREND REMOVAL/DIFF, currently not fully implemented
    dataset = sli
    # dataset = map(lambda x,y: x-y, sli[1:], sli[:-1])
    # m = np.mean(sli)
    # dataset = map(lambda x: x - m, sli)

    aaa = fftconvolve(dataset, dataset[::-1], mode='full')
    bbb = aaa[len(aaa)/2:]
    ccc = map(lambda x: x / bbb[0], bbb)
    return ccc[1:]


# Generic function to generate waterfall plot,
# uses make_it as the 'make a line' function
def make_thing_water(tiel, ds, width, stride, per_item, times, make_it):
    xr = min(len(ds), width)
    yr = int(1 + (len(ds) / float(stride)))

    line = make_it(ds[0:width])
    xr = len(line)

    img = np.zeros(shape=(yr, xr), dtype=float)
    print(np.shape(img))

    for (i, val) in enumerate(line):
        img[0, i] = val

    for yl in range(1, yr - 1):
        st = stride * yl
        line = make_it(ds[st:st + width])
        for (i, val) in enumerate(line):
            img[yl, i] = val

    plt.figure(tiel)

    plt.ion()
    plt.clf()

    # plt.hist(img.flatten(), bins=1500)
    plt.imshow(img)
    plt.xlabel("")
    plt.ylabel("Offset in time %g:%g" % (times[0], times[1]))

    # plt.legend(scatterpoints=1, framealpha=0.5)
    # plt.grid(True)
    plt.title(tiel)
    plt.draw()


# Class to isolate some processing from the UI aspects
class Thedata:

    def __init__(self):
        self.clear()
        self.results = Results()

        self.clearcache()
        self.dpmap = {}

        # Rudimentary 'on demand' computation/cache
        self.how = {
            "classes": self.mkclasses,
            "classes_sorted": self.mksortclasses,
            "classstats_fft": self.mkfft,
            "classstats_cor": self.mkacor,
            "stats_perclass": self.mkstats,
            "stats_global": self.mkglob,
        }

    def sel_classes(self):
        return map(lambda x: self.dpmap[self.var["class_select"].get(x)],
                   self.var["class_select"].curselection())

    def clear(self):
        print("Clearing dataset")
        self.dataset = []

    # Open a PCAP file, adding its packets to the 'database'
    def open(self, fname):
        print("Opening %s" % fname)
        p = pcap.pcap(fname)
        bpf = self.var["BPF"].get()
        if len(bpf) > 0:
            print("Using BPF filter '%s'" % bpf)
            p.setfilter(bpf)

        def cb(timestamp, pack):
            pkt = pack[:]
            self.dataset.append((timestamp, pkt))
        p.loop(0, cb)
        print("Total %d packets in dataset" % len(self.dataset))

    # Make FFTs for all selected classes
    def mkfft(self):
        classes = self.get_data("classes")
        whitelist = self.sel_classes()
        wind = map(lambda x: self.var["Window"].get(x),
                   self.var["Window"].curselection()
                   )[0]

        slicesize = self.var["SlicePackets"].get()
        slicestart = self.var["SliceStart"].get()
        for (k, c) in classes.items():
            if k not in whitelist:
                continue
            times = {}
            for a in c:
                if a[0] in times:
                    times[a[0]] += 1
                else:
                    times[a[0]] = 1
            timeline = sorted(times.items())[slicestart:slicesize+slicestart]

            if timeline[-1][0] == timeline[0][0]:
                self.results.log("Timerange of length 1",
                                 level="event",
                                 more="Class: %s\nPackets: %d" % (k, len(c)))
                continue

            timeline_slotted = np.zeros(self.var["Resolution"].get(),
                                        dtype=int)
            # timeline_slotted = np.zeros(self.resolution)

            offset = timeline[0][0]
            per_item = (timeline[-1][0] - offset) / \
                       (self.var["Resolution"].get() - 1)

            for a in timeline:
                i = int(np.floor(((a[0] - offset) / per_item)))
                timeline_slotted[i] += 1

            ff = make_fft(timeline_slotted, window=wind)

            magn = [(x * x + y * y) ** 0.5 for (x, y) in zip(ff.real, ff.imag)]

            fftdata = mmm(magn[1:])
            fftdata["data"] = magn
            fftdata["freqs"] = np.fft.fftfreq(
                len(timeline_slotted),
                d=1.0 / per_item
                )[1:1 + len(magn)]
            # rate / fftlen

            if "classstats_fft" not in self.cache.keys():
                cs = self.set_data("classstats_fft", {})
            else:
                cs = self.get_data("classstats_fft")
            cs[k] = fftdata

    def mkacorwater(self):
        self.mkthingwater(make_acor)

    def mkfftwater(self):
        wind = map(lambda x: self.var["Window"].get(x),
                   self.var["Window"].curselection())[0]
        print(wind)

        fn = lambda x: make_fft(x, window=wind)
        self.mkthingwater(fn)

    # Make a waterfall plot of all selected classes
    # using 'make_it' as the 'make a line' function
    def mkthingwater(self, make_it):
        classes = self.get_data("classes")
        whitelist = self.sel_classes()
        for (k, c) in classes.items():
            slicest = self.var["SliceStart"].get()
            slicelen = self.var["SlicePackets"].get()
            width = self.var["WaterWidth"].get()
            stride = self.var["WaterStride"].get()
            if k not in whitelist:
                continue
            times = {}
            for a in c:
                if a[0] in times:
                    times[a[0]] += 1
                else:
                    times[a[0]] = 1
            timeline = sorted(times.items())[slicest: slicest+slicelen]

            thebegin = timeline[0][0]
            theend = timeline[-1][0]

            if theend == thebegin:
                self.results.log("Timerange of length 1",
                                 level="event",
                                 more="Class: %s\nPackets: %d" % (k, len(c)))
                continue

            timeline_slotted = np.zeros(self.var["Resolution"].get(),
                                        dtype=int)

            timespan = theend - thebegin
            per_item = timespan / (self.var["Resolution"].get() - 1)

            # width = int(width / per_item)
            stride = int(width * stride)
            print(width, stride)

            for a in timeline:
                i = int(np.floor((a[0] - thebegin) / per_item))
                timeline_slotted[i] += 1

            ti = "%s: %d:%d [%s.%s] (%d seconds)" % (
                k,
                slicest,
                slicest + len(timeline),
                time.ctime(thebegin),
                time.ctime(theend),
                int(theend-thebegin))
            print(ti)

            make_thing_water(
                ti,
                timeline_slotted,
                width,
                stride,
                per_item,
                (thebegin, theend),
                make_it)

    # Make autocorrelations for all selected classes
    def mkacor(self):
        classes = self.get_data("classes")
        whitelist = self.sel_classes()
        slicesize = self.var["SlicePackets"].get()
        slicestart = self.var["SliceStart"].get()
        for (k, c) in classes.items():
            if k not in whitelist:
                continue
            times = {}
            for a in c:
                if a[0] in times:
                    times[a[0]] += 1
                else:
                    times[a[0]] = 1
            timeline = sorted(times.items())[slicestart:slicesize+slicestart]

            thebegin = timeline[0][0]
            theend = timeline[-1][0]
            if theend == thebegin:
                self.results.log("Timerange of length 1",
                                 level="event",
                                 more="Class: %s\nPackets: %d" % (k, len(c)))
                continue

            timeline_slotted = np.zeros(self.var["Resolution"].get(),
                                        dtype=int)

            timespan = theend - thebegin
            per_item = timespan / (self.var["Resolution"].get() - 1)

            for a in timeline:
                i = int(np.floor((a[0] - thebegin) / per_item))
                timeline_slotted[i] += 1

            cor = make_acor(timeline_slotted)

            cordata = mmm(cor)
            cordata["data"] = cor
            cordata["per_item"] = per_item

            if "classstats_cor" not in self.cache.keys():
                cs = self.set_data("classstats_cor", {})
            else:
                cs = self.get_data("classstats_cor")
            cs[k] = cordata

    # Make all classes/process packets to classify
    def mkclasses(self):
        print("Making classes")

        # Clear old'ns
        self.var["class_select"].delete(0, tk.END)

        StaticClasses.clear()
        StaticClasses.setformat(self.var["Filter"].get())
        classes = {}
        for (ts, p) in self.dataset:
            c = StaticClasses.findclass(p)
            ent = (ts, p, c)
            if c in classes:
                classes[c].append(ent)
            else:
                classes[c] = [ent]
        print("Got %d packets over %d classes" % (
            len(self.dataset),
            len(classes.keys())))

        self.set_data("classes", classes)
        self.classlist()

    # Make the sorted liss of classes (for use in UI)
    def mksortclasses(self):
        delim = self.var["Delim"].get()
        num = self.var["Sortkey"].get()

        if delim[0] == "$":
            names = delim[1:]
            where, val = names.split(":")
            cvs = self.get_data(where)

            def sp(classna):
                return str(cvs[classna][val])

            classes = self.get_data('classes')
            classes_s = sorted(classes.items(), key=sp)
            self.set_data("classes_sorted", classes_s)
            return

        sp_def = lambda x: x[0].split(delim)[num]
        sp = {
            "NUM": lambda x: len(x[1]),
            "BYTES": lambda x: sum(map(lambda y: len(y[1]), classes[x[0]])),
            }.get(delim, sp_def)
        print(sp)

        classes = self.get_data('classes')
        classes_s = sorted(classes.items(), key=sp)
        self.set_data("classes_sorted", classes_s)

    # Update the classlist for the UI
    def classlist(self):
        self.mksortclasses()  # Resort
        classes = self.get_data('classes_sorted')
        self.var["class_select"].delete(0, tk.END)

        delim = self.var["Delim"].get()
        num = self.var["Sortkey"].get()
        if delim[0] == "$":
            names = delim[1:]
            where, val = names.split(":")
            cvs = self.get_data(where)
            mkidd = lambda classna, t: "%s: %d::  %s" % (
                classna,
                len(t),
                str(cvs[classna][val]))
        else:
            mkidd = lambda k, v: "%s_%d" % (k, len(v))

        for (k, v) in classes:
            idd = mkidd(k, v)
            self.dpmap[idd] = k

            self.var["class_select"].insert(tk.END, idd)
        self.var["class_select"].selection_set(0, tk.END)

    def mkstats(self):
        print("Making stats")
        self.stats_perclass = {}
        classes = self.get_data("classes")
        whitelist = self.sel_classes()

        if "stats_perclass" not in self.cache:
            spc = self.set_data("stats_perclass", {})
        else:
            spc = self.get_data("stats_perclass")

        for (k, c) in classes.items():
            if k not in whitelist:
                continue

            cs = {}
            cs["packetcount"] = len(c)

            cs["size_min"] = min(map(lambda x: len(x[1]), c))
            cs["size_max"] = max(map(lambda x: len(x[1]), c))
            cs["size_std"] = np.std(map(lambda x: len(x[1]), c))
            cs["size_mean"] = np.mean(map(lambda x: len(x[1]), c))

            cs["time_min"] = min(map(lambda x: (x[0]), c))
            cs["time_max"] = max(map(lambda x: (x[0]), c))
            cs["time_std"] = np.std(map(lambda x: (x[0]), c))
            cs["time_mean"] = np.mean(map(lambda x: (x[0]), c))
            cs["time_span"] = max(
                map(lambda x: (x[0]), c)
                ) + min(
                map(lambda x: (x[0]), c)
                )

            spc[k] = cs

    def mkglob(self):
        c = self.dataset
        cs = {}

        cs["size_min"] = min(map(lambda x: len(x[1]), c))
        cs["size_max"] = max(map(lambda x: len(x[1]), c))
        cs["size_std"] = np.std(map(lambda x: len(x[1]), c))
        cs["size_mean"] = np.mean(map(lambda x: len(x[1]), c))

        cs["time_min"] = min(map(lambda x: (x[0]), c))
        cs["time_max"] = max(map(lambda x: (x[0]), c))
        cs["time_std"] = np.std(map(lambda x: (x[0]), c))
        cs["time_mean"] = np.mean(map(lambda x: (x[0]), c))
        cs["time_span"] = max(
            map(lambda x: (x[0]), c)
            ) + min(
            map(lambda x: (x[0]), c)
            )

        self.set_data("stats_global", cs)

    def storecs(self, fn):
        # Class sizes
        f = open(fn, "wb")
        ii = []
        iii = []
        classes = self.get_data("classes_sorted")
        whitelist = self.sel_classes()
        for (k, v) in classes:
            if k not in whitelist:
                continue
            f.write("%s, %g" % (k, len(v)))
        f.close()

    def plot_data(self, key, title, perx):
        whitelist = self.sel_classes()
        csf = self.get_data(key)
        for (k, i) in csf.items():
            if k not in whitelist:
                continue
            data = i["data"]
            plt.figure(title)
            plt.ion()
            plt.clf()
            xd = []
            if "freqs" in i.keys():
                xd = i["freqs"]
                print(xd[:50], xd[:-50])
            else:
                xd = map(lambda x: x*i["per_item"], xrange(len(data)))
            plt.scatter(xd, data)
            # plt.plot(xd, data)
            plt.xlabel(perx)
            plt.ylabel("Magnitude")

            plt.legend(scatterpoints=1, framealpha=0.5)
            plt.grid(True)
            plt.title("%s: %s" % (title, k))
            plt.xscale(self.var["Xscale"].get())
            plt.yscale(self.var["Yscale"].get())
            plt.draw()
            print("Plotting %s,%s,%s" % (key, title, k))
        print("Plotted %s,%s" % (key, title))

    def plotfft(self):
        self.plot_data("classstats_fft", "FFT", "Hamsters")

    def plotacor(self):
        self.plot_data("classstats_cor", "Autocorrelation", "Seconds")

    def plotcs(self):
        ii = []
        iii = []
        classes = self.get_data("classes_sorted")
        whitelist = self.sel_classes()
        for (k, v) in classes:
            if k not in whitelist:
                continue
            ii.append(k)
            iii.append(len(v))
        plt.figure("Classsizes")
        plt.ion()
        plt.clf()
        try:
            v = map(float, ii)
            plt.bar(map(lambda x: x - 0.4, v), iii)
            plt.xlabel("Value")
        except ValueError as xx:
            print(xx)
            plt.bar(np.arange(-0.4, len(iii) + 0.5, 1.0)[:len(iii)], iii)
            plt.xlabel("Index")

        plt.legend(scatterpoints=1, framealpha=0.5)
        plt.grid(True)
        plt.ylabel("Count")
        plt.title("Class sizes")
        plt.xscale(self.var["Xscale"].get())
        plt.yscale(self.var["Yscale"].get())
        plt.draw()

    def storepackets(self, fn):
        f = open(fn, "wb")
        w = dpkt.pcap.Writer(f)

        whitelist = self.sel_classes()
        for (k, ds) in self.get_data("classes"):
            if k not in whitelist:
                continue
            for (ts, p) in ds:
                w.writepkt(p, ts=ts)
        f.close()

    def plotstats(self):
        print("Plotting")
        n = [0]

        def p(data, ti):
            if len(data) < 1:
                return
            plt.figure(ti)
            plt.ion()
            plt.clf()
            plt.hist(data, bins=150)
            plt.legend(scatterpoints=1, framealpha=0.5)
            plt.grid(True)
            plt.xlabel("Bin")
            plt.ylabel("Count")
            plt.title(ti)
            plt.xscale("linear")
            plt.yscale("linear")
            plt.draw()

            n[0] += 1
        a = {"smean": [], "smin": [], "smax": [], "sstd": [], "svar": [],
             "tmean": [], "tmin": [], "tmax": [], "tstd": [], "tvar": [],
             "tspan": [], "packetcount": []}

        spc = self.get_data("stats_perclass")
        whitelist = self.sel_classes()
        for (c, l) in spc.items():
            if c not in whitelist:
                continue
            a["packetcount"].append(l["packetcount"])
            print(l)
            a["smin"].append(l["size_min"])
            a["smax"].append(l["size_max"])
            a["sstd"].append(l["size_std"])
            a["smean"].append(l["size_mean"])

            a["tmin"].append(l["time_min"])
            a["tmax"].append(l["time_max"])
            a["tstd"].append(l["time_std"])
            a["tmean"].append(l["time_mean"])
            a["tspan"].append(l["time_span"])
        for (k, v) in a.items():
            p(v, k)

    def storestats(self, fn):
        print("Storing stats")
        f = open(fn, "w")
        whitelist = self.sel_classes()
        f.write("Global, packets, %d\n" % len(self.dataset))
        f.write("Global, Format, %s\n" % self.var["Filter"].get())

        classes = self.get_data("classes")
        f.write("Global, Classes, %d\n" % len(classes.keys()))

        spc = self.get_data("stats_perclass")
        for (k, v) in spc.items():
            if k not in whitelist:
                continue
            for (n, g) in v.items():
                f.write("%s, %s, %g\n" % (k, n, g))

        gs = self.get_data("stats_global")
        for (n, g) in gs.items():
            f.write("global, %s, %g\n" % (n, g))

        csf = self.get_data("classstats_fft")
        for (k, i) in csf.items():
            if k not in whitelist:
                continue
            f.write("FFT: %s:%s\n" % (k, i.__repr__()))

        csc = self.get_data("classstats_cor")
        for (k, i) in csc.items():
            if k not in whitelist:
                continue
            f.write("Cor: %s:%s\n" % (k, i.__repr__()))
        f.close()

    def storeperclass(self, fn):
        print("Storing per class")
        classes = self.get_data("classes")
        whitelist = self.sel_classes()
        for (k, v) in classes.items():
            if k not in whitelist:
                continue
            print(k)
            f = open("%s_%s" % (fn, k), "wb")

            w = dpkt.pcap.Writer(f)
            for (ts, p, c) in v:
                w.writepkt(p, ts=ts)
            f.close()

    def set_data(self, k, v):
        print("CSET %s")
        self.cache[k] = v
        return v

    def get_data(self, what):
        print("Getting data %s" % what)
        if what in self.cache.keys():
            rv = self.cache[what]
            print("Hit %s" % what)
            return rv
        met = self.how[what]
        print("Making data %s" % str(met))
        print("Miss '%s': by %s" % (what, str(met)))
        print(":"+"\n:".join(self.cache.keys()))
        met()
        rv = self.cache[what]
        # print("Thing: '%s'" % rv)
        return rv

    def clearcache(self):
        self.cache = {}


# Class to use as a sort of logger, not fully completed
class Results:

    def __init__(self, fname=None):
        if fname:
            self.fd = open(fname, "wb")
        else:
            self.fd = False
        self.hits = {"info": 0, "note": 0, "event": 0, "alert": 0, "error": 0}

    def log(self, msg, level="info", more=""):
        # LOG
        self.hits[level] += 1
        t = "[%s]%s" % (level, msg)
        m = "[%s]%s (%s)\n" % (level, msg, more)
        print(m)
        if self.fd:
            self.fd.write(m)
        if level == "alert":
            m = tkMessageBox.showerror(t, m)
        # m = tkMessageBox(m, option=tkMessageBox.OK, icon=tkMessageBox.ERROR)


# Class responsible for the UI, interaction and 'bootup'
class TheUI(tk.Frame):
    stick = tk.N + tk.S + tk.W + tk.E

    def setDefaults(self):
        # self.option_add("*button*font", self.font)
        self.option_add("*font", self.font)

    def clear(self):
        self.pd = Thedata()
        self.pd.var = self.var

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        # Change this to fix Font
        # (development was on a high pixel density screen)
        self.font = tkFont.Font(family='sans', size=-35, weight='normal')

        self.setDefaults()

        self.grid(sticky=TheUI.stick)
        self.createWidgets()

        # Install a killswitch (not perfect, only on the main window)
        self.bind_all("<Key-Print>", self.doquit)

        self.clear()
        self.curfiles = []

    def __handleprocess(self):
        self.pd.mkclasses()
        self.pd.mkstats()
        self.pd.analyze()

    def __handleplots(self):
        self.pd.plotstats()

    def __handleplotcs(self):
        self.pd.plotcs()

    def createMenu(self):
        self.menu = tk.Menu(self.top, tearoff=0)
        self.top["menu"] = self.menu

        self.submenu_file = tk.Menu(self.menu)
        self.menu.add_cascade(menu=self.submenu_file, label="File")
        self.submenu_file.add_command(label="Open", command=self.__openhandler)
        self.submenu_file.add_command(
            label="Open+", command=self.__openphandler)
        self.submenu_file.add_command(
            label="Reopen!", command=self.__reopenhandler)
        self.submenu_file.add_command(
            label="Export blob", command=self.__exportb)
        self.submenu_file.add_command(
            label="Export classes", command=self.__exportc)
        self.submenu_file.add_command(
            label="Export classsizes", command=self.__exportz)
        self.submenu_file.add_command(
            label="Export stats", command=self.__exports)
        self.submenu_file.add_command(
            label="Quit", command=self.doquit)

        self.submenu_about = tk.Menu(self.menu)
        self.menu.add_cascade(menu=self.submenu_about, label="About")
        self.submenu_about.add_command(
            label="About", command=self.__abouthandler)

    def doquit(self, *args, **kargs):
        import sys
        sys.exit(0)

    def __abouthandler(self):
        print("You clicked a button, how about that?")

    def __openhandler(self):
        self.clear()
        self.curfiles = []
        f = self.__promptfile()
        self.curfiles.append(f)
        self.pd.open(f)

    def __reopenhandler(self):
        print("Reopen files")
        self.clear()
        for a in self.curfiles:
            print("--%s" % a)
            self.pd.open(a)
        print("Done")

    def __openphandler(self):
        f = self.__promptfile()
        self.curfiles.append(f)
        self.pd.open(f)

    def __promptfile(self):
        o = tkFileDialog.askopenfilename(defaultextension=".pcap", filetypes=[
            ("PCAP", "*.pcap*"),
            ("tcpdump", "*.tcpdump*"),
            ("All", "*"), ],
            initialdir="~/DATA/", title="Open packet data")
        return o

    def __exportz(self):
        fn = tkFileDialog.asksaveasfilename(defaultextension=".pcap")
        self.pd.storecs(fn)

    def __exportb(self):
        fn = tkFileDialog.asksaveasfilename(defaultextension=".pcap")
        self.pd.storepackets(fn)

    def __exportc(self):
        fn = tkFileDialog.asksaveasfilename(defaultextension=".pcap")
        self.pd.storeperclass(fn)

    def __exports(self):
        fn = tkFileDialog.asksaveasfilename(defaultextension=".txt")
        self.pd.storestats(fn)

    def createWidgets(self):
        self.top = self.winfo_toplevel()
        self.top.columnconfigure(0, weight=1)
        self.top.rowconfigure(0, weight=1)

        self.rowconfigure(0, weight=0)  # Menu
        self.rowconfigure(1, weight=0)  # Labels

        # self.rowconfigure(2, weight=1)
        # self.columnconfigure(0, weight=1)

        self.createMenu()

        # Col 0: buttons
        # Col 1: Labels
        # Col 2: entry boxes
        # For button
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=0)
        self.columnconfigure(2, weight=2)
        self.buttons = []
        self.var = {}

        def mkbutton(ti, fn):
            b = tk.Button(self, text=ti, command=fn)
            # b.grid(row=len(self.buttons), column=0, sticky=tk.W)
            b.grid(row=1 + len(self.buttons), column=0, sticky=TheUI.stick)
            self.buttons.append(b)

        lb = tk.Label(self, text="Settings:")
        lb.grid(row=1, column=2)

        def mksetting(ti, v, de):
            rn = 2 + len(self.var.items())
            val = v()
            val.set(de)
            self.var[ti] = val

            l = tk.Label(self, text="%s:" % ti)
            l.grid(row=rn, column=1, sticky=TheUI.stick)

            if v == tk.StringVar:
                e = tk.Entry(self, textvariable=val)
            else:
                e = tk.Entry(self, textvariable=val)
            e.grid(row=rn, column=2, sticky=TheUI.stick)
        self.listboxes = []

        def mklistbox(varname, varcontent=[], high=10, sel=-1, mode="many"):
            # dict_mode = {"many": tk.MULTIPLE, "one": tk.SINGLE}
            dict_mode = {"many": tk.EXTENDED, "one": tk.SINGLE}
            container = tk.Frame(self)

            c = 3 + len(self.listboxes)
            self.columnconfigure(c, weight=1)  # Content
            container.grid(row=2, column=c, sticky=TheUI.stick)
            self.rowconfigure(2, weight=2)
            self.columnconfigure(c, weight=2)

            container.rowconfigure(0, weight=2)
            container.columnconfigure(0, weight=2)

            lbl = tk.Label(self, text=varname+":")
            lbl.grid(row=1, column=c, sticky=tk.N+tk.W+tk.E+tk.S)

            # s1 = tk.Scrollbar(container, orient=tk.HORIZONTAL)
            # s1.grid(row=2, column=0, sticky=tk.E+tk.W)
            #
            # s2 = tk.Scrollbar(container, orient=tk.VERTICAL)
            # s2.grid(row=1, column=1, sticky=tk.N+tk.S)

            var = tk.StringVar()
            lb = tk.Listbox(
                container,
                exportselection=0,
                activestyle="underline",
                # xscrollcommand=s1.set,
                # yscrollcommand=s2.set,
                selectmode=dict_mode[mode],
                height=high,
                listvariable=var,
                # state=tk.NORMAL
                )
            # print(var)
            lb.grid(row=0,
                    column=0,
                    rowspan=high,
                    columnspan=1,
                    sticky=TheUI.stick)
            self.listboxes.append(lb)
            self.var[varname] = lb

            for a in varcontent:
                lb.insert(tk.END, a)

            if sel == -1:
                return
            if sel == "all":
                lb.selection_set(0, tk.END)
            else:
                lb.selection_set(sel)

        mkbutton("(re)Do Cache", self.__clearcache)
        mkbutton("Classlist", lambda: self.pd.classlist())
        # mkbutton("Proc bulk", self.__handleprocess)
        mkbutton("Plot classes", self.__handleplotcs)
        mkbutton("Plot stats", self.__handleplots)
        mkbutton("Plot FFT", lambda: self.pd.plotfft())
        mkbutton("Plot acor", lambda: self.pd.plotacor())
        mkbutton("Plot awater", lambda: self.pd.mkacorwater())
        mkbutton("Plot fwater", lambda: self.pd.mkfftwater())

        mksetting("Filter", tk.StringVar, "prot_tcplowport_udplowport")
        mksetting("BPF", tk.StringVar, "")
        mksetting("Xscale", tk.StringVar, "linear")
        mksetting("Yscale", tk.StringVar, "linear")
        mksetting("Resolution", tk.IntVar, 100000)
        mksetting("SlicePackets", tk.IntVar, 10000)
        mksetting("SliceStart", tk.IntVar, 0)
        mksetting("Delim", tk.StringVar, "NUM")
        mksetting("Sortkey", tk.IntVar, -1)
        mksetting("WaterWidth", tk.IntVar, 1000)
        mksetting("WaterStride", tk.DoubleVar, 1.0)

        mklistbox(
            "class_select",
            varcontent=map(str, range(100)),
            sel="all",
            mode="many")
        mklistbox("Window", varcontent=dict_win.keys(), mode="one", sel=1)

    def __clearcache(self):
        self.pd.clearcache()
        print(len(self.pd.get_data("classes")))


def main():
    app = TheUI()
    app.master.title("UI")
    app.mainloop()


if __name__ == "__main__":
    main()
