# Synthetic placements with an independent oracle for the alnqc suite.
import math
import os
import tempfile

L1, L2 = 50, 40
NAMES = ["t0", "t1", "t2", "chr1", "u0", "u1"]
LENGTHS = [1000, 1000, 500, 10000, 800, 800]
TXGENES = [
    {"tx": "t0", "gene": "g0", "geneName": "G0", "biotype": "protein_coding"},
    {"tx": "t1", "gene": "g0", "geneName": "G0", "biotype": "protein_coding"},
    {"tx": "t2", "gene": "g1", "geneName": "G1", "biotype": "rRNA"},
]


def _hit(tx, pos, fwd=True, frag_len=200, mates=0, weight=1.0):
    return {"tx": tx, "pos": pos, "fwd": fwd, "fragLen": frag_len, "mates": mates,
            "score": 100, "weight": weight, "format": -1}


def fragments():
    """(kind, FragHits) in order."""
    out = []
    out += [{"hits": [], "decoy": False}] * 10
    out += [{"hits": [], "decoy": True}] * 5
    out += [{"hits": [_hit(0, 10 * i)], "decoy": False} for i in range(20)]
    out += [{"hits": [_hit(0, 300, fwd=False), _hit(1, 700, fwd=False)], "decoy": False} for _ in range(12)]
    out += [{"hits": [_hit(0, 50), _hit(2, 100, weight=0.5)], "decoy": False} for _ in range(4)]
    out += [{"hits": [_hit(2, 100 + i, frag_len=0, mates=2)], "decoy": False} for i in range(6)]
    out += [{"hits": [_hit(2, 480, frag_len=0, mates=1)], "decoy": False}]
    # references no gene names: several of them are not one gene
    out += [{"hits": [_hit(4, 100), _hit(5, 100)], "decoy": False} for _ in range(3)]
    out += [{"hits": [_hit(4, 300)], "decoy": False} for _ in range(2)]
    return out


def ref_names():
    return NAMES


def ref_lengths():
    return LENGTHS


def tx_genes():
    return TXGENES


def placements():
    return fragments()


def write_pairs():
    d = tempfile.mkdtemp(prefix="alnqc-test-")
    p = os.path.join(d, "r")
    n = len(fragments())
    for k, l in ((1, L1), (2, L2)):
        with open("%s_%d.fq" % (p, k), "w") as f:
            for i in range(n):
                f.write("@r%d\n%s\n+\n%s\n" % (i, "A" * l, "I" * l))
    return p


def _points(length):
    return [(length - 1) * k // 100 for k in range(1, 101)]


def expected():
    gene = {"t0": "g0", "t1": "g0", "t2": "g1"}
    btype = {"t0": "protein_coding", "t1": "protein_coding", "t2": "rRNA"}
    e = {"fragments": 0, "unmapped": 0, "decoy": 0, "oneTx": 0, "oneGene": 0, "multiGene": 0}
    bio = {}
    cov = {}
    nfr = {}
    for f in fragments():
        e["fragments"] += 1
        hs = f["hits"]
        if not hs:
            e["decoy" if f["decoy"] else "unmapped"] += 1
            continue
        txs = {h["tx"] for h in hs}
        genes = {gene.get(NAMES[t]) for t in txs}
        one_gene = len(txs) == 1 or (len(genes) == 1 and None not in genes)
        if len(txs) == 1:
            e["oneTx"] += 1
        elif one_gene:
            e["oneGene"] += 1
        else:
            e["multiGene"] += 1
        tot = sum(h["weight"] for h in hs)
        for h in hs:
            b = btype.get(NAMES[h["tx"]])
            if b is not None:
                bio[b] = bio.get(b, 0.0) + h["weight"] / tot
        if not one_gene:
            continue
        for t in txs:
            nfr[t] = nfr.get(t, 0) + 1
        for h in hs:
            t, n = h["tx"], LENGTHS[h["tx"]]
            covered = set()
            if h["mates"] == 0 and h["fragLen"] > 0:
                left, right = (L1, L2) if h["fwd"] else (L2, L1)
                covered |= set(range(h["pos"], h["pos"] + left))
                covered |= set(range(h["pos"] + h["fragLen"] - right, h["pos"] + h["fragLen"]))
            else:
                l = L2 if h["mates"] == 2 else L1
                covered |= set(range(h["pos"], h["pos"] + l))
            c = cov.setdefault(t, [0] * 100)
            for k, p in enumerate(_points(n)):
                if p in covered:
                    c[k] += 1
    return e, bio, cov, nfr


def _tin(c):
    tot = sum(c)
    h = -sum(x / tot * math.log(x / tot) for x in c if x > 0)
    return 100 * math.exp(h) / 100


def check(q, refs, body, tins):
    """'' when the summary, gene body coverage and TINs match the oracle."""
    e, bio, cov, nfr = expected()
    for k, v in e.items():
        if q[k] != v:
            return "%s: got %r, want %r" % (k, q[k], v)
    got_bio = {n: c for n, c in zip(refs["biotypeNames"], q["biotypeCounts"])}
    if set(got_bio) != set(bio) or any(abs(got_bio[k] - bio[k]) > 1e-9 for k in bio):
        return "biotypes: got %r, want %r" % (got_bio, bio)
    got = {p["tx"]: list(p["points"]) for p in q["txPoints"]}
    if got != cov:
        diff = [t for t in set(got) | set(cov) if got.get(t) != cov.get(t)]
        return "points differ on %r" % diff
    if [p["fragments"] for p in q["txPoints"]] != [nfr[t] for t in sorted(nfr)]:
        return "fragment counts differ"
    want_body = [sum(cov[t][k] for t in cov) for k in range(100)]
    if list(body) != want_body:
        return "gene body differs"
    want_tin = [(t, _tin(cov[t])) for t in sorted(cov) if nfr[t] >= 10]
    if len(tins) != len(want_tin) or any(a[0] != b[0] or abs(a[1] - b[1]) > 1e-9 for a, b in zip(tins, want_tin)):
        return "TIN: got %r, want %r" % (tins, want_tin)
    return ""


def check_reports(body_text, tin_text, bio_text, rrna_text, origin_text):
    lines = body_text.strip().split("\n")
    if not lines[0].startswith("Percentile\t1\t2") or not lines[1].startswith("S\t") or len(lines[1].split("\t")) != 101:
        return "geneBodyCoverage text"
    t = tin_text.strip().split("\n")
    if t[0] != "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)" or not t[1].startswith("S\t"):
        return "TIN text"
    rows = [l.split("\t") for l in bio_text.split("\n") if l and not l.startswith("#")]
    if sorted(r[0] for r in rows) != ["protein_coding", "rRNA"]:
        return "biotype text %r" % rows
    _, bio, _, _ = expected()
    want = 100 * bio["rRNA"] / sum(bio.values())
    last = rrna_text.strip().split("\n")[-1].split("\t")
    if last[0] != "'S'" or abs(float(last[1]) - want) > 1e-9:
        return "rRNA text %r" % last
    o = origin_text.strip().split("\n")[-1].split("\t")
    if o[0] != "S" or abs(sum(float(x) for x in o[1:]) - 100) > 1e-9:
        return "origin text %r" % o
    return ""
