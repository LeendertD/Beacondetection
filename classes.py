import dpkt
import zlib


def __init__(self):
    pass


def to_ip(r):
    oc = map(str, map(ord, r))
    return ".".join(oc)


def srcdst(p):
    a = [to_ip(p.dst), to_ip(p.src)]
    a.sort()
    return ":".join(a)


# Use zlib compression to estimate entropy
def entropy_zlib(ip, raw, decs):
    data = ip.data
    mlen, sc = map(float, decs.split(","))
    if len(data) < mlen:
        return 0.0
    olen = 1.0 * len(raw)
    clen = len(zlib.compress(raw))
    rat = clen / olen
    rat *= sc
    return int(rat)


# Check the length of an IP packet againts a range
def iplen(ip, raw, args):
    mi, ma = args.split(",")
    return mi <= ip.len < ma


# Use an expert rule to determine class 'title'/'name'/'identifier'
class StaticClasses:
    class_history = {}
    packets = 0

    fmt = [lambda x, y, z: "NO_FORMAT"]

    @classmethod
    def clear(self):
        StaticClasses.class_history = {}
        StaticClasses.packets = 0
        StaticClasses.fmt = [lambda x, y, z: "NO_FORMAT"]

    @classmethod
    def setformat(self, fmt):
        trans = {
            "": lambda x, y: "blank",
            "src": lambda ip, raw: to_ip(ip.src),
            "dst": lambda ip, raw: to_ip(ip.dst),
            "srcdst": lambda ip, raw: srcdst(ip),
            "len": iplen,
            "datalen": lambda ip, c: len(ip.data),
            "prot=": lambda ip, raw, arg: ip.p == arg,
            "prot": lambda ip, raw: ip.p,
            "ttl": lambda ip, raw: ip.ttl,
            "tos": lambda ip, raw: ip.tos,
            "icmp": lambda ip, raw: isinstance(ip.data, dpkt.icmp.ICMP),
            "icmp.type": lambda ip, raw:
                isinstance(ip.data, dpkt.icmp.ICMP) and ip.data.type,
            "icmp.code": lambda ip, raw:
                isinstance(ip.data, dpkt.icmp.ICMP) and ip.data.code,

            "tcpsport": lambda ip, raw:
                isinstance(ip.data, dpkt.tcp.TCP) and ip.data.sport,
            "tcpdport": lambda ip, raw:
                isinstance(ip.data, dpkt.tcp.TCP) and ip.data.dport,
            "tcplowport": lambda ip, raw:
                isinstance(ip.data, dpkt.tcp.TCP) and
                min(ip.data.dport, ip.data.sport),
            "tcphighport": lambda ip, raw:
                isinstance(ip.data, dpkt.tcp.TCP) and
                max(ip.data.dport, ip.data.sport),

            "udpsport": lambda ip, raw:
                isinstance(ip.data, dpkt.udp.UDP) and ip.data.sport,
            "udpdport": lambda ip, raw:
                isinstance(ip.data, dpkt.udp.UDP) and ip.data.dport,
            "udplowport": lambda ip, raw:
                isinstance(ip.data, dpkt.udp.UDP) and
                min(ip.data.dport, ip.data.sport),
            "udphighport": lambda ip, raw:
                isinstance(ip.data, dpkt.udp.UDP) and
                max(ip.data.dport, ip.data.sport),
            "entropy.zlib": entropy_zlib,
            }
        self.fmt = []
        for a in fmt.split("_"):
            if ":" in a:
                cmd, args = a.split(":", 1)
                fn = lambda x, y: trans[cmd](x, y, args)
            else:
                fn = trans[a]

            self.fmt.append(fn)

    @classmethod
    def findclass(self, raw):
        e = dpkt.ethernet.Ethernet(raw)
        if isinstance(e.data, dpkt.ip.IP):
            ip = e.data

            items = map(lambda x: str(x(ip, raw)), self.fmt)
            clid = "_".join(items)
        else:
            clid = "Unk"

        if clid in self.class_history:
            self.class_history[clid] += 1
        else:
            self.class_history[clid] = 1
        self.packets += 1
        return clid
