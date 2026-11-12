"""Figures of the report (CERN ROOT) and the ROOT file of the sweep."""
import math

import numpy as np
import ROOT

from common import FIG, RES, describe, plain, system, tlabel

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Categorical slots in fixed order (validated palette, light mode).
HEX = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
       "#6250d6", "#e34948"]
COL = [ROOT.TColor.GetColor(h) for h in HEX]
INK = ROOT.TColor.GetColor("#0b0b0b")
INK2 = ROOT.TColor.GetColor("#52514e")
GRID = ROOT.TColor.GetColor("#e4e3df")
BAND = ROOT.TColor.GetColor("#dcdad3")


def style():
    s = ROOT.gStyle
    s.SetOptStat(0)
    s.SetOptTitle(0)
    s.SetTextFont(42)
    for ax in "xyz":
        s.SetLabelFont(42, ax)
        s.SetTitleFont(42, ax)
        s.SetLabelSize(0.035, ax)
        s.SetTitleSize(0.040, ax)
        s.SetLabelColor(INK2, ax)
        s.SetTitleColor(INK, ax)
        s.SetAxisColor(INK2, ax)
    s.SetLegendFont(42)
    s.SetLegendBorderSize(0)
    s.SetLegendFillColor(0)
    s.SetLegendTextSize(0.035)
    if hasattr(s, "SetLegendFillStyle"):
        s.SetLegendFillStyle(0)  # transparent: grid lines stay visible
    s.SetFrameLineColor(INK2)
    s.SetGridColor(GRID)
    s.SetGridStyle(1)
    s.SetPadTickX(0)
    s.SetPadTickY(0)
    s.SetEndErrorSize(0)
    s.SetCanvasColor(0)
    s.SetPadColor(0)
    s.SetFrameFillColor(0)


style()
KEEP = []  # PyROOT: keep drawn objects alive until the canvas is saved


def keep(*o):
    KEEP.extend(o)
    return o[0] if len(o) == 1 else o


def save(c, name):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        c.SaveAs(str(FIG / f"{name}.{ext}"))
    KEEP.clear()


def text(x, y, s, size=0.035, color=INK, align=11, ndc=True):
    t = ROOT.TLatex(x, y, s)
    if ndc:
        t.SetNDC()
    t.SetTextSize(size)
    t.SetTextColor(color)
    t.SetTextAlign(align)
    t.SetTextFont(42)
    t.Draw()
    return keep(t)


def fig_ranking(rows):
    """Dot plot: the theory floor (Paule-Mandel) of each configuration with its
    68 % Q-profile interval, both systems (a diagnostic; the grades are the
    goodness-of-fit tests)."""
    n = len(rows)
    c = ROOT.TCanvas("rank", "", 1100, 160 + 22 * n)
    c.SetLeftMargin(0.30)
    c.SetRightMargin(0.04)
    c.SetTopMargin(60.0 / (160 + 22 * n))
    c.SetBottomMargin(100.0 / (160 + 22 * n))
    c.SetLogx()
    c.SetGridx()
    xmax = max(max(r["ccbar_floor"], r["bbar_floor"]) for r in rows) * 1.6
    h = keep(ROOT.TH2F("hr", "", 10, 1.0, xmax, n, 0, n))
    h.GetXaxis().SetTitle("theory floor #sigma_{th} (MeV)  #leftarrow better")
    h.GetXaxis().SetMoreLogLabels()
    h.GetXaxis().SetNoExponent()
    for k, r in enumerate(rows):
        lbl = describe(r["tag"]) + ("  (degenerate)" if r["degenerate"] else "")
        h.GetYaxis().SetBinLabel(n - k, lbl)
    # pixel-sized fonts (precision 3): the canvas height grows with n
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(15)
    h.GetYaxis().SetTickLength(0)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(20)
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(22)
    h.GetXaxis().SetTitleOffset(1.6)
    h.Draw()
    for j, (sname, lab) in enumerate((("ccbar", "c#bar{c}"),
                                      ("bbar", "b#bar{b}"))):
        g = keep(ROOT.TGraphAsymmErrors())
        for k, r in enumerate(rows):
            f = max(r[f"{sname}_floor"], 1.01)
            g.SetPoint(k, f, n - k - 0.5 + (0.12 if j else -0.12))
            lo, hi = r.get(f"{sname}_floor_lo", math.nan), r.get(f"{sname}_floor_hi", math.nan)
            if math.isfinite(lo) and math.isfinite(hi):
                g.SetPointError(k, max(f - max(lo, 1.0), 0), max(hi - f, 0), 0, 0)
        g.SetMarkerStyle(20 if j == 0 else 21)
        g.SetMarkerSize(1.1)
        g.SetMarkerColor(COL[j])
        g.SetLineColor(COL[j])
        g.Draw("PZ same")
    leg = keep(ROOT.TLegend(0.72, 1 - 55.0 / (160 + 22 * n), 0.96, 1 - 8.0 / (160 + 22 * n)))
    leg.SetTextFont(43)
    leg.SetTextSize(20)
    leg.SetNColumns(2)
    leg.AddEntry(KEEP[-3], "c#bar{c}", "p")
    leg.AddEntry(KEEP[-2], "b#bar{b}", "p")
    leg.Draw()
    save(c, "ranking_floor")


def fig_residuals(runs, tags, sname, name):
    """Model - PDG for every state, a few models side by side."""
    s0 = system(runs[tags[0]], sname)
    labels = [x["label"] for x in s0["states"]]
    n = len(labels)
    c = ROOT.TCanvas(name, "", 1200, 650)
    c.SetBottomMargin(0.24)
    c.SetLeftMargin(0.09)
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.12)
    c.SetGridy()
    lo, hi = -130, 150
    h = keep(ROOT.TH1F(f"h{name}", "", n, 0, n))
    for k, lbl in enumerate(labels):
        h.GetXaxis().SetBinLabel(k + 1, tlabel(lbl))
    h.GetXaxis().LabelsOption("v")
    h.GetXaxis().SetLabelSize(0.040)
    h.SetMinimum(lo)
    h.SetMaximum(hi)
    h.GetYaxis().SetTitle("M_{model} #minus M_{PDG} (MeV)")
    h.Draw("AXIS")
    nfit = sum(1 for x in s0["states"] if x["in_fit"])
    box = keep(ROOT.TBox(nfit, lo, n, hi))
    box.SetFillColorAlpha(BAND, 0.35)
    box.Draw()
    text(nfit + 0.1, hi - 12, "above threshold: predictions", 0.032, INK2,
         11, ndc=False)
    text(0.1, hi - 12, "fitted", 0.032, INK2, 11, ndc=False)
    zero = keep(ROOT.TLine(0, 0, n, 0))
    zero.SetLineColor(INK2)
    zero.Draw()
    leg = keep(ROOT.TLegend(0.09, 0.89, 0.98, 0.99))
    leg.SetNColumns(len(tags))
    for j, t in enumerate(tags):
        s = system(runs[t], sname)
        g = keep(ROOT.TGraphErrors())
        for k, x in enumerate(s["states"]):
            dx = (j - (len(tags) - 1) / 2) * 0.14
            v = 1e3 * (x["model_GeV"] - x["pdg_GeV"])
            g.SetPoint(k, k + 0.5 + dx, max(lo + 2, min(hi - 2, v)))
            g.SetPointError(k, 0, 1e3 * math.hypot(x["sigma_par_GeV"],
                                                    s["sigma_th_GeV"]))
        g.SetMarkerStyle(20)
        g.SetMarkerSize(1.0)
        g.SetMarkerColor(COL[j])
        g.SetLineColor(COL[j])
        g.SetLineWidth(2)
        g.Draw("PZ same")
        leg.AddEntry(g, f"{describe(t)} (#sigma_{{th}} = {1e3 * s['sigma_th_GeV']:.1f})", "p")
    leg.Draw()
    h.Draw("AXIS same")
    save(c, name)


def fig_spectrum(run, sname, name, title):
    """Level diagram by J^PC column: PDG (black line, grey band +- sigma_exp)
    and model (coloured line, box +- sqrt(sigma_th^2 + sigma_par^2 +
    sigma_num^2))."""
    s = system(run, sname)
    cols = ["1S0", "3S1", "1P1", "3P0", "3P1", "3P2", "3D1", "3D2"]
    jpc = {"1S0": "0^{#minus+}", "3S1": "1^{#minus#minus}",
           "1P1": "1^{+#minus}", "3P0": "0^{++}", "3P1": "1^{++}",
           "3P2": "2^{++}", "3D1": "1^{#minus#minus} (D)",
           "3D2": "2^{#minus#minus}"}
    used = [cc for cc in cols if any(x["channel"] == cc for x in s["states"])]
    ms = [x["pdg_GeV"] for x in s["states"]] + [x["model_GeV"] for x in s["states"]]
    lo, hi = min(ms) - 0.08, max(ms) + 0.08
    c = ROOT.TCanvas(name, "", 1000, 700)
    c.SetLeftMargin(0.10)
    c.SetRightMargin(0.03)
    c.SetTopMargin(0.08)
    h = keep(ROOT.TH1F(f"h{name}", "", len(used), 0, len(used)))
    for k, cc in enumerate(used):
        h.GetXaxis().SetBinLabel(k + 1, jpc[cc])
    h.GetXaxis().SetLabelSize(0.045)
    h.SetMinimum(lo)
    h.SetMaximum(hi)
    h.GetYaxis().SetTitle("mass (GeV)")
    h.Draw("AXIS")
    thr = {"ccbar": (3.7297, "D#bar{D} threshold"),
           "bbar": (10.5590, "B#bar{B} threshold")}[sname]
    ln = keep(ROOT.TLine(0, thr[0], len(used), thr[0]))
    ln.SetLineStyle(2)
    ln.SetLineColor(INK2)
    ln.Draw()
    text(0.05, thr[0] + 0.01, thr[1], 0.030, INK2, 11, ndc=False)
    sth = s.get("sigma_th_GeV", 0.0)
    for x in s["states"]:
        k = used.index(x["channel"])
        e = x.get("sigma_exp_GeV") or 0.0
        if e > 0.0005:
            band = keep(ROOT.TBox(k + 0.12, x["pdg_GeV"] - e, k + 0.48, x["pdg_GeV"] + e))
            band.SetFillColorAlpha(INK2, 0.35)
            band.Draw()
        pdg = keep(ROOT.TLine(k + 0.12, x["pdg_GeV"], k + 0.48, x["pdg_GeV"]))
        pdg.SetLineColor(INK)
        pdg.SetLineWidth(3)
        pdg.Draw()
        m = x["model_GeV"]
        if m is None:
            continue
        err = math.sqrt((sth * x.get("floor_w", 1.0)) ** 2 + (x.get("sigma_par_GeV") or 0) ** 2
                        + (x.get("sigma_num_GeV") or 0) ** 2)
        err = max(err, 0.0015)
        mod = keep(ROOT.TBox(k + 0.52, m - err, k + 0.88, m + err))
        mod.SetFillColorAlpha(COL[0] if x["in_fit"] else COL[1], 0.45)
        mod.Draw()
        mid = keep(ROOT.TLine(k + 0.52, m, k + 0.88, m))
        mid.SetLineColor(COL[0] if x["in_fit"] else COL[1])
        mid.SetLineWidth(2)
        mid.Draw()
    leg = keep(ROOT.TLegend(0.36, 0.80, 0.97, 0.92))
    l1 = keep(ROOT.TLine())
    l1.SetLineWidth(3)
    l1.SetLineColor(INK)
    b1 = keep(ROOT.TBox())
    b1.SetFillColor(COL[0])
    b2 = keep(ROOT.TBox())
    b2.SetFillColor(COL[1])
    leg.AddEntry(l1, "PDG (band: #pm#sigma_{exp})", "l")
    leg.AddEntry(b1, "model, fitted (box: #pm#sqrt{#sigma_{th}^{2}+#sigma_{par}^{2}+#sigma_{num}^{2}})", "f")
    leg.AddEntry(b2, "model, prediction (same)", "f")
    leg.Draw()
    text(0.10, 0.94, title, 0.035)
    save(c, name)


def fig_alpha_s():
    """alpha_s(Q)/alpha_s(0) of every running form and its Gaussian-sum fit."""
    ROOT.gInterpreter.Declare("""
    #include <cmath>
    double qg_frozen(double Q2) { return 4*M_PI/(25./3.*std::log((Q2+4*0.25)/0.09)); }
    double qg_apt(double Q2) { double x=Q2/0.09; if (std::fabs(x-1)<1e-6) return 4*M_PI/(25./3.)*0.5;
                                return 4*M_PI/(25./3.)*(1/std::log(x)+1/(1-x)); }
    double qg_rich(double Q2) { return 4*M_PI/(9.*std::log1p(Q2/0.16)); }
    double qg_gi(double Q2) { double a[3]={0.25,0.15,0.20}, g[3]={0.5,std::sqrt(10.)/2,std::sqrt(1000.)/2};
                               double s=0; for(int k=0;k<3;++k) s+=a[k]*std::exp(-Q2/(4*g[k]*g[k])); return s; }
    """)
    c = ROOT.TCanvas("as", "", 1000, 650)
    c.SetLogx()
    c.SetGridy()
    c.SetLeftMargin(0.10)
    c.SetRightMargin(0.03)
    h = keep(ROOT.TH1F("has", "", 10, 0.02, 50))
    h.SetMinimum(0)
    h.SetMaximum(1.6)
    h.GetXaxis().SetTitle("Q (GeV)")
    h.GetYaxis().SetTitle("#alpha_{s}(Q^{2})")
    h.GetXaxis().SetMoreLogLabels()
    h.GetXaxis().SetNoExponent()
    h.Draw("AXIS")
    forms = [("GI Gaussian sum (C/D, A/B gi)", ROOT.qg_gi, 0),
             ("frozen, #Lambda = 0.30, M_{g} = 0.50 GeV", ROOT.qg_frozen, 1),
             ("APT, #Lambda = 0.30 GeV", ROOT.qg_apt, 2),
             ("Richardson, #Lambda = 0.40 GeV (confining, not fitted)",
              ROOT.qg_rich, 3)]
    leg = keep(ROOT.TLegend(0.12, 0.13, 0.52, 0.40))
    leg.SetFillStyle(0)
    for lab, fn, j in forms:
        g = keep(ROOT.TGraph())
        m = 0
        for k in range(300):
            q = 0.02 * (50 / 0.02) ** (k / 299)
            if fn(q * q) < 1.6:
                g.SetPoint(m, q, fn(q * q))
                m += 1
        g.SetLineColor(COL[j])
        g.SetLineWidth(2)
        g.Draw("L same")
        leg.AddEntry(g, lab, "l")
    leg.Draw()
    text(0.96, 0.86, "A/B fit #alpha_{s}(0) with the shape fixed;", 0.028, INK2, 31)
    text(0.96, 0.81, "C/D 'free' fits the three GI amplitudes", 0.028, INK2, 31)
    save(c, "alpha_s_forms")


def fig_loo_vs_floor(rows):
    c = ROOT.TCanvas("loo", "", 1000, 700)
    c.SetLogx()
    c.SetLogy()
    c.SetGrid()
    c.SetLeftMargin(0.11)
    c.SetRightMargin(0.03)
    h = keep(ROOT.TH2F("hloo", "", 10, 1, 500, 10, 1, 500))
    h.GetXaxis().SetTitle("theory floor #sigma_{th} (MeV)")
    h.GetYaxis().SetTitle("exact leave-one-out RMS (MeV)")
    h.GetXaxis().SetMoreLogLabels()
    h.GetYaxis().SetMoreLogLabels()
    h.GetXaxis().SetNoExponent()
    h.GetYaxis().SetNoExponent()
    h.Draw()
    diag = keep(ROOT.TLine(1, 1, 500, 500))
    diag.SetLineColor(INK2)
    diag.SetLineStyle(2)
    diag.Draw()
    leg = keep(ROOT.TLegend(0.13, 0.75, 0.50, 0.90))
    for j, (sname, lab) in enumerate((("ccbar", "c#bar{c}"),
                                      ("bbar", "b#bar{b}"))):
        g = keep(ROOT.TGraphAsymmErrors())
        k = 0
        for r in rows:
            x, y = r[f"{sname}_floor"], r[f"{sname}_loo_rms"]
            if math.isfinite(y) and x > 0:
                g.SetPoint(k, max(x, 1.05), max(y, 1.05))
                lo, hi = r.get(f"{sname}_floor_lo", math.nan), r.get(f"{sname}_floor_hi", math.nan)
                ey = y / math.sqrt(2 * r[f"{sname}_N"])  # sampling error of an RMS
                g.SetPointError(k, x - lo if math.isfinite(lo) else 0,
                                hi - x if math.isfinite(hi) else 0, min(ey, y - 1.0), ey)
                k += 1
        g.SetMarkerStyle(20 if j == 0 else 21)
        g.SetMarkerColor(COL[j])
        g.SetLineColor(COL[j])
        g.SetMarkerSize(1.1)
        g.Draw("PZ same")
        leg.AddEntry(g, lab, "p")
    leg.Draw()
    text(0.55, 0.18, "dashed: LOO RMS = #sigma_{th}", 0.03, INK2)
    save(c, "loo_vs_floor")


def fig_widths(runs, tags):
    """Model / PDG for every annihilation width with a PDG value: LO (open)
    and NLO with the scale band mu = M/4 ... M (filled, error bar)."""
    keys = []
    for t in tags:
        for sname in ("ccbar", "bbar"):
            for w in system(runs[t], sname).get("widths", []):
                k = (w["label"], w["process"])
                if w.get("pdg_keV") and k not in keys:
                    keys.append(k)
    if not keys:
        return
    c = ROOT.TCanvas("w", "", 1200, 650)
    c.SetLeftMargin(0.08)
    c.SetRightMargin(0.02)
    c.SetBottomMargin(0.22)
    c.SetLogy()
    c.SetGridy()
    h = keep(ROOT.TH1F("hw", "", len(keys), 0, len(keys)))
    for k, (lab, proc) in enumerate(keys):
        h.GetXaxis().SetBinLabel(k + 1, tlabel(lab) + (" #gamma#gamma" if proc == "gamma gamma"
                                                        else " ee"))
    h.GetXaxis().SetLabelSize(0.04)
    h.GetXaxis().LabelsOption("v")
    h.SetMinimum(0.05)
    h.SetMaximum(20)
    h.GetYaxis().SetTitle("#Gamma(model) / #Gamma(PDG)")
    h.GetYaxis().SetMoreLogLabels()
    h.GetYaxis().SetNoExponent()
    h.Draw("AXIS")
    one = keep(ROOT.TLine(0, 1, len(keys), 1))
    one.SetLineColor(INK2)
    one.SetLineStyle(2)
    one.Draw()
    leg = keep(ROOT.TLegend(0.09, 0.80, 0.97, 0.97))
    leg.SetNColumns(2)
    for j, t in enumerate(tags):
        g0 = keep(ROOT.TGraph())
        g1 = keep(ROOT.TGraphAsymmErrors())
        for sname in ("ccbar", "bbar"):
            for w in system(runs[t], sname).get("widths", []):
                k = (w["label"], w["process"])
                if k not in keys or not w.get("pdg_keV"):
                    continue
                x = keys.index(k) + 0.5 + (j - (len(tags) - 1) / 2) * 0.15
                pdg = w["pdg_keV"]
                g0.SetPoint(g0.GetN(), x, w["lo_keV"] / pdg)
                n = g1.GetN()
                g1.SetPoint(n, x, w["nlo_keV"] / pdg)
                g1.SetPointError(n, 0, 0, (w["nlo_keV"] - w["nlo_band_lo_keV"]) / pdg,
                                 (w["nlo_band_hi_keV"] - w["nlo_keV"]) / pdg)
        g0.SetMarkerStyle(24)
        g1.SetMarkerStyle(20)
        for g in (g0, g1):
            g.SetMarkerColor(COL[j])
            g.SetLineColor(COL[j])
            g.SetMarkerSize(1.1)
        g0.Draw("P same")
        g1.Draw("PZ same")
        leg.AddEntry(g1, f"{describe(t)}: open LO, filled NLO (band #mu = M/4 ... M)", "p")
    leg.Draw()
    save(c, "widths")


def fig_pulls(run, name):
    c = ROOT.TCanvas(name, "", 900, 600)
    c.SetLeftMargin(0.11)
    c.SetRightMargin(0.03)
    h = keep(ROOT.TH1F(f"h{name}", "", 16, -4, 4))
    hp = keep(ROOT.TH1F(f"hp{name}", "", 16, -4, 4))
    for sname in ("ccbar", "bbar"):
        s = system(run, sname)
        for x in s["states"]:
            if x["pull"] is None:
                continue
            (h if x["in_fit"] else hp).Fill(max(-3.99, min(3.99, x["pull"])))
    h.SetFillColor(COL[0])
    h.SetLineColor(0)
    hp.SetFillColor(COL[1])
    hp.SetLineColor(0)
    st = keep(ROOT.THStack(f"s{name}", ""))
    st.Add(h)
    st.Add(hp)
    st.Draw("hist")
    st.GetXaxis().SetTitle("pull (leave-one-out for fitted states)")
    st.GetYaxis().SetTitle("states")
    f = keep(ROOT.TF1(f"f{name}", "[0]*exp(-0.5*x*x)", -4, 4))
    f.SetParameter(0, (h.Integral() + hp.Integral()) * 0.5 / math.sqrt(2 * math.pi))
    f.SetLineColor(INK2)
    f.SetLineStyle(2)
    f.Draw("same")
    leg = keep(ROOT.TLegend(0.65, 0.75, 0.97, 0.92))
    leg.AddEntry(h, "fitted (LOO pull)", "f")
    leg.AddEntry(hp, "predicted", "f")
    leg.AddEntry(f, "N(0,1)", "l")
    leg.Draw()
    save(c, name)


def fig_params(runs, tags):
    """Fitted heavy-quark mass, coupling and string tension across models."""
    pars = [("m1", "m_{c} (GeV)", ("ccbar",)), ("m1", "m_{b} (GeV)", ("bbar",)),
            ("alpha_s", "#alpha_{s} (or #alpha_{s}(0))", ("ccbar", "bbar")),
            ("b", "b (GeV^{2})", ("ccbar", "bbar"))]
    c = ROOT.TCanvas("pars", "", 1400, 160 + 28 * len(tags))
    c.Divide(len(pars), 1, 0.002, 0.0)
    for ip, (pn, title, systems) in enumerate(pars):
        pad = c.cd(ip + 1)
        pad.SetLeftMargin(0.50 if ip == 0 else 0.05)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.16)
        pad.SetTopMargin(0.12)
        pad.SetGridx()
        vals = []
        for t in tags:
            for sname in systems:
                an = system(runs[t], sname)["analysis"]
                if pn in an["names"]:
                    k = an["names"].index(pn)
                    vals.append((t, sname, an["value"][k], an["sigma"][k]))
        if not vals:
            continue
        lo = min(v - 2 * e for _, _, v, e in vals)
        hi = max(v + 2 * e for _, _, v, e in vals)
        n = len(tags)
        h = keep(ROOT.TH2F(f"hp{pn}{ip}", "", 10, lo, hi, n, 0, n))
        for k, t in enumerate(tags):
            h.GetYaxis().SetBinLabel(n - k, describe(t) if ip == 0 else "")
        h.GetYaxis().SetLabelSize(0.085 if ip == 0 else 0)
        h.GetYaxis().SetTickLength(0)
        h.GetXaxis().SetTitle(title)
        h.GetXaxis().SetTitleSize(0.06)
        h.GetXaxis().SetLabelSize(0.05)
        h.GetXaxis().SetNdivisions(505)
        h.Draw()
        for j, sname in enumerate(("ccbar", "bbar")):
            if sname not in systems:
                continue
            g = keep(ROOT.TGraphErrors())
            m = 0
            for t, sn, v, e in vals:
                if sn == sname:
                    g.SetPoint(m, v, n - tags.index(t) - 0.5 + (0.12 if j else -0.12))
                    g.SetPointError(m, e, 0)
                    m += 1
            g.SetMarkerStyle(20 if j == 0 else 21)
            g.SetMarkerColor(COL[j])
            g.SetLineColor(COL[j])
            g.SetLineWidth(2)
            g.Draw("PZ same")
            if ip == 2:
                if j == 0:
                    leg = keep(ROOT.TLegend(0.30, 0.885, 0.98, 0.995))
                    leg.SetNColumns(2)
                    leg.SetTextSize(0.06)
                leg.AddEntry(g, "c#bar{c}" if j == 0 else "b#bar{b}", "p")
                if j == 1:
                    leg.Draw()
    c.cd(0)
    save(c, "parameters")


def fig_learning_curves(cv_runs, tags):
    c = ROOT.TCanvas("lc", "", 1300, 560)
    c.Divide(2, 1, 0.004, 0.0)
    for ip, sname in enumerate(("ccbar", "bbar")):
        pad = c.cd(ip + 1)
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.16)
        pad.SetBottomMargin(0.13)
        pad.SetLogy()
        pad.SetGridy()
        curves = []
        for t in tags:
            lc = system(cv_runs[t], sname)["analysis"].get("learning_curve") if t in cv_runs else None
            if lc and lc["train_size"]:
                curves.append((t, lc))
        if not curves:
            continue
        xs = [m for _, lc in curves for m in lc["train_size"]]
        h = keep(ROOT.TH1F(f"hlc{ip}", "", 10, min(xs) - 0.5, max(xs) + 0.5))
        h.SetMinimum(0.8)
        h.SetMaximum(100)
        h.GetXaxis().SetTitle("states in the training set")
        h.GetYaxis().SetTitle("RMS of model #minus PDG (MeV)")
        h.GetYaxis().SetMoreLogLabels()
        h.GetYaxis().SetNoExponent()
        h.Draw("AXIS")
        text(0.13, 0.95, {"ccbar": "c#bar{c}", "bbar": "b#bar{b}"}[sname], 0.05, INK, 13)
        leg = keep(ROOT.TLegend(0.13, 0.85, 0.97, 0.995))
        leg.SetNColumns(2)
        leg.SetTextSize(0.033)
        for j, (t, lc) in enumerate(curves):
            for key, style, lab in (("train_rms_GeV", 2, "train"),
                                    ("test_median_GeV", 1, "test (median)")):
                g = keep(ROOT.TGraph())
                for k, (m, v) in enumerate(zip(lc["train_size"], lc[key])):
                    if v is not None and math.isfinite(v):
                        g.SetPoint(g.GetN(), m, 1e3 * v)
                g.SetLineColor(COL[j])
                g.SetMarkerColor(COL[j])
                g.SetLineStyle(style)
                g.SetLineWidth(2)
                g.SetMarkerStyle(20 if style == 1 else 24)
                g.Draw("LP same")
                if style == 1:
                    leg.AddEntry(g, describe(t), "lp")
        leg.Draw()
    c.cd(0)
    save(c, "learning_curves")


def fig_bgrid(brows):
    tags = sorted({r["tag"] for r in brows})
    if not tags:
        return
    r1s = sorted({r["r1"] for r in brows})
    rms_ = sorted({r["rmax"] for r in brows})
    c = ROOT.TCanvas("bg", "", 1200, 1000)
    c.Divide(2, len(tags), 0.004, 0.004)
    ROOT.gStyle.SetPalette(ROOT.kBird)
    for it, tag in enumerate(tags):
        for js, sname in enumerate(("ccbar", "bbar")):
            pad = c.cd(2 * it + js + 1)
            pad.SetLogz()
            pad.SetRightMargin(0.16)
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.14)
            pad.SetTopMargin(0.10)
            h = keep(ROOT.TH2F(f"hbg{it}{js}", "", len(r1s), 0, len(r1s), len(rms_), 0,
                               len(rms_)))
            for i, a in enumerate(r1s):
                h.GetXaxis().SetBinLabel(i + 1, f"{a:g}")
            for j, b in enumerate(rms_):
                h.GetYaxis().SetBinLabel(j + 1, f"{b:g}")
            for r in brows:
                if r["tag"] == tag and r["system"] == sname and math.isfinite(r["trace_excess_MeV"]):
                    h.SetBinContent(r1s.index(r["r1"]) + 1, rms_.index(r["rmax"]) + 1,
                                    max(r["trace_excess_MeV"], 1e-3))
            h.SetMinimum(1e-3)
            h.GetXaxis().SetTitle("r_{1} (GeV^{-1})")
            h.GetYaxis().SetTitle("r_{max} (GeV^{-1})")
            h.GetZaxis().SetTitle("level sum above the best (MeV)")
            h.GetXaxis().SetTitleOffset(1.1)
            ROOT.gStyle.SetPaintTextFormat(".2g")
            h.SetMarkerSize(1.3)
            h.Draw("COLZ TEXT")
            text(0.14, 0.93, f"{describe(tag)},  {'c#bar{c}' if sname == 'ccbar' else 'b#bar{b}'}",
                 0.05, INK, 11)
    c.cd(0)
    save(c, "basis_grid")


def fig_minimizers(mrows):
    c = ROOT.TCanvas("mz", "", 1000, 750)
    c.SetLogx()
    c.SetLogy()
    c.SetGrid()
    c.SetLeftMargin(0.11)
    c.SetRightMargin(0.03)
    h = keep(ROOT.TH2F("hmz", "", 10, 1, 1000, 10, 1, 1000))
    h.GetXaxis().SetTitle("final #sigma_{th}: best minimum any method found (MeV)")
    h.GetYaxis().SetTitle("#sigma_{th} reached by one method (MeV)")
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetMoreLogLabels()
        ax.SetNoExponent()
    h.Draw()
    d = keep(ROOT.TLine(1, 1, 1000, 1000))
    d.SetLineColor(INK2)
    d.SetLineStyle(2)
    d.Draw()
    leg = keep(ROOT.TLegend(0.50, 0.24, 0.95, 0.40))
    for j, (key, lab, mk) in enumerate((("lm_single_floor", "Levenberg-Marquardt, default start", 20),
                                        ("lm_multistart_floor", "Levenberg-Marquardt, best of 9 starts", 21),
                                        ("minuit_floor", "Minuit2 Migrad, default start", 25))):
        g = keep(ROOT.TGraph())
        k = 0
        for r in mrows:
            x, y = r["final_floor"], r[key]
            if math.isfinite(x) and math.isfinite(y):
                g.SetPoint(k, max(x, 1.05), max(y, 1.05))
                k += 1
        g.SetMarkerStyle(mk)
        g.SetMarkerColor(COL[j])
        g.SetMarkerSize(1.2)
        g.Draw("P same")
        leg.AddEntry(g, lab, "p")
    leg.Draw()
    text(0.50, 0.18, "above the diagonal: stuck in a local minimum", 0.03, INK2)
    save(c, "minimizers")


def fig_error_methods(erows):
    sel = [r for r in erows if math.isfinite(r["bootstrap"]) and r["linear"] > 0]
    if not sel:
        return
    n = len(sel)
    H = 160 + 20 * n
    c = ROOT.TCanvas("em", "", 1000, H)
    c.SetLeftMargin(0.36)
    c.SetRightMargin(0.03)
    c.SetTopMargin(60.0 / H)
    c.SetBottomMargin(100.0 / H)
    c.SetGridx()
    h = keep(ROOT.TH2F("hem", "", 10, 0, 3, n, 0, n))
    for k, r in enumerate(sel):
        flag = (f"  (MCMC #hat{{R}} = {r['mcmc_rhat']:.1f})"
                if math.isfinite(r["mcmc_rhat"]) and r["mcmc_rhat"] > 1.1 else "")
        sl = {"ccbar": "c#bar{c}", "bbar": "b#bar{b}"}[r["system"]]
        pl = {"m1": "m_{Q}", "alpha_s": "#alpha_{s}", "sigma": "#sigma",
              "sigma0": "#sigma_{0}"}.get(r["param"], r["param"])
        h.GetYaxis().SetBinLabel(n - k, f"{plain(r['tag'])}  {sl}  {pl}{flag}")
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(14)
    h.GetYaxis().SetTickLength(0)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(20)
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(22)
    h.GetXaxis().SetTitleOffset(1.5)
    h.GetXaxis().SetTitle("error / linear (J^{T}J)^{-1} error")
    h.Draw()
    one = keep(ROOT.TLine(1, 0, 1, n))
    one.SetLineColor(INK2)
    one.SetLineStyle(2)
    one.Draw()
    leg = keep(ROOT.TLegend(0.36, 1 - 55.0 / H, 0.97, 1 - 8.0 / H))
    leg.SetTextFont(43)
    leg.SetTextSize(20)
    leg.SetNColumns(3)
    for j, (key, lab, mk) in enumerate((("minos", "Minos (half width)", 22),
                                        ("bootstrap", "bootstrap std", 20),
                                        ("mcmc", "MCMC std", 21))):
        g = keep(ROOT.TGraph())
        m = 0
        for k, r in enumerate(sel):
            v = (0.5 * (r["minos_hi"] - r["minos_lo"]) if key == "minos"
                 else r[key])
            if math.isfinite(v):
                g.SetPoint(m, min(2.95, v / r["linear"]), n - k - 0.5)
                m += 1
        g.SetMarkerStyle(mk)
        g.SetMarkerColor(COL[j])
        g.SetMarkerSize(1.1)
        g.Draw("P same")
        leg.AddEntry(g, lab, "p")
    leg.Draw()
    save(c, "error_methods")


def fig_tau(names, M):
    n = len(names)
    c = ROOT.TCanvas("tau", "", 1000, 900)
    c.SetLeftMargin(0.27)
    c.SetBottomMargin(0.27)
    c.SetRightMargin(0.13)
    c.SetTopMargin(0.04)
    h = keep(ROOT.TH2F("htau", "", n, 0, n, n, 0, n))
    for i in range(n):
        h.GetXaxis().SetBinLabel(i + 1, names[i])
        h.GetYaxis().SetBinLabel(i + 1, names[i])
        for j in range(n):
            h.SetBinContent(i + 1, j + 1, M[i, j])
    h.SetMinimum(-1)
    h.SetMaximum(1)
    h.GetXaxis().LabelsOption("v")
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetLabelFont(43)
        ax.SetLabelSize(17)
    h.GetZaxis().SetLabelSize(0.03)
    # diverging palette through a neutral midpoint (tau = 0)
    stops = np.array([0.0, 0.5, 1.0])
    reds = np.array([0.92, 0.95, 0.16])
    greens = np.array([0.41, 0.95, 0.47])
    blues = np.array([0.20, 0.95, 0.84])
    ROOT.TColor.CreateGradientColorTable(3, stops, reds, greens, blues, 99)
    h.SetMarkerSize(1.1)
    ROOT.gStyle.SetPaintTextFormat(".2f")
    h.Draw("COLZ TEXT")
    save(c, "criteria_tau")
    ROOT.gStyle.SetPalette(ROOT.kBird)


def fig_weights(wrows, keys):
    """Weight that each criterion gives each physical model."""
    tags = [r["tag"] for r in wrows]
    n = len(tags)
    H = 160 + 30 * n
    c = ROOT.TCanvas("wts", "", 1000, H)
    c.SetLeftMargin(0.26)
    c.SetRightMargin(0.03)
    c.SetTopMargin(70.0 / H)
    c.SetBottomMargin(90.0 / H)
    c.SetGridx()
    h = keep(ROOT.TH2F("hwts", "", 10, 0, 0.6, n, 0, n))
    for k, t in enumerate(tags):
        h.GetYaxis().SetBinLabel(n - k, describe(t))
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(16)
    h.GetYaxis().SetTickLength(0)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(18)
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(20)
    h.GetXaxis().SetTitleOffset(1.6)
    h.GetXaxis().SetTitle("model weight")
    h.Draw()
    leg = keep(ROOT.TLegend(0.26, 1 - 62.0 / H, 0.97, 1 - 6.0 / H))
    leg.SetNColumns(3)
    leg.SetTextFont(43)
    leg.SetTextSize(16)
    marks = [20, 21, 22, 23, 33, 34]
    for j, (key, lab) in enumerate(keys):
        g = keep(ROOT.TGraph())
        m = 0
        for k, r in enumerate(wrows):
            v = r.get(key, float("nan"))
            if isinstance(v, float) and math.isfinite(v):
                g.SetPoint(m, min(v, 0.595), n - k - 0.5 + 0.12 * (j - (len(keys) - 1) / 2))
                m += 1
        g.SetMarkerStyle(marks[j % len(marks)])
        g.SetMarkerColor(COL[j % len(COL)])
        g.SetMarkerSize(1.2)
        g.Draw("P same")
        leg.AddEntry(g, lab, "p")
    leg.Draw()
    save(c, "model_weights")


def fig_averaged(arows, keys):
    """Averaged prediction minus PDG for every unfitted state, per weighting."""
    n = len(arows)
    c = ROOT.TCanvas("avg", "", 1200, 600)
    c.SetBottomMargin(0.24)
    c.SetLeftMargin(0.08)
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.10)
    c.SetGridy()
    lo = min(r[f"{k}_diff"] - r[f"{k}_unc"] for r in arows for k in keys) - 10
    hi = max(r[f"{k}_diff"] + r[f"{k}_unc"] for r in arows for k in keys) + 10
    h = keep(ROOT.TH1F("havg", "", n, 0, n))
    for i, r in enumerate(arows):
        h.GetXaxis().SetBinLabel(i + 1, tlabel(r["state"]))
    h.SetMinimum(lo)
    h.SetMaximum(hi)
    h.GetXaxis().LabelsOption("v")
    h.GetXaxis().SetLabelSize(0.045)
    h.GetYaxis().SetTitle("averaged model #minus PDG (MeV)")
    h.GetYaxis().SetTitleOffset(0.8)
    h.Draw("AXIS")
    zero = keep(ROOT.TLine(0, 0, n, 0))
    zero.SetLineColor(INK2)
    zero.Draw()
    leg = keep(ROOT.TLegend(0.08, 0.91, 0.98, 0.99))
    leg.SetNColumns(len(keys))
    for j, k in enumerate(keys):
        g = keep(ROOT.TGraphErrors())
        for i, r in enumerate(arows):
            g.SetPoint(i, i + 0.5 + 0.15 * (j - (len(keys) - 1) / 2), r[f"{k}_diff"])
            g.SetPointError(i, 0, r[f"{k}_unc"])
        g.SetMarkerStyle(20 + j)
        g.SetMarkerColor(COL[j])
        g.SetLineColor(COL[j])
        g.SetLineWidth(2)
        g.Draw("PZ same")
        leg.AddEntry(g, {"AICc": "AICc weights", "pBMA": "pseudo-BMA (LOO)",
                         "stack": "stacking (LOO)", "equal": "equal weights"}.get(k, k), "p")
    leg.Draw()
    save(c, "averaged_predictions")


def fig_selection(rows):
    """The four goodness-of-fit p-values of every model (masses and decays,
    both systems) against alpha = 0.05; filled: passes the physics checks,
    open: fails one. Sorted as the selection table."""
    keys = [("p_mass_ccbar", "masses c#bar{c}"), ("p_mass_bbar", "masses b#bar{b}"),
            ("p_decay_ccbar", "decays c#bar{c}"), ("p_decay_bbar", "decays b#bar{b}")]
    n = len(rows)
    H = 150 + 24 * n
    c = ROOT.TCanvas("sel", "", 1100, H)
    c.SetLeftMargin(0.30)
    c.SetRightMargin(0.04)
    c.SetTopMargin(60.0 / H)
    c.SetBottomMargin(90.0 / H)
    c.SetLogx()
    c.SetGridx()
    h = keep(ROOT.TH2F("hsel", "", 10, 1e-6, 1.5, n, 0, n))
    h.GetXaxis().SetTitle("p-value (Pearson #chi^{2})")
    for k, r in enumerate(rows):
        h.GetYaxis().SetBinLabel(n - k, describe(r["tag"]) + ("  #checkmark" if r["verdict"] == "pass" else ""))
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetLabelFont(43)
        ax.SetTitleFont(43)
    h.GetYaxis().SetLabelSize(15)
    h.GetYaxis().SetTickLength(0)
    h.GetXaxis().SetLabelSize(18)
    h.GetXaxis().SetTitleSize(20)
    h.GetXaxis().SetTitleOffset(1.6)
    h.Draw()
    line = keep(ROOT.TLine(0.05, 0, 0.05, n))
    line.SetLineColor(INK2)
    line.SetLineStyle(2)
    line.Draw()
    leg = keep(ROOT.TLegend(0.31, 1 - 55.0 / H, 0.96, 1 - 6.0 / H))
    leg.SetNColumns(4)
    leg.SetTextFont(43)
    leg.SetTextSize(15)
    for j, (key, lab) in enumerate(keys):
        for physics_ok in (True, False):
            g = keep(ROOT.TGraph())
            for k, r in enumerate(rows):
                gate2 = not any(not f.startswith(("mass_", "decay_")) for f in r["fails"])
                if key in r and math.isfinite(r[key]) and gate2 == physics_ok:
                    g.SetPoint(g.GetN(), max(r[key], 1.1e-6), n - k - 0.5 + (j - 1.5) * 0.15)
            g.SetMarkerStyle((20, 21, 22, 23)[j] if physics_ok else (24, 25, 26, 32)[j])
            g.SetMarkerColor(COL[j])
            g.SetMarkerSize(1.0)
            if g.GetN():
                g.Draw("P same")
            if physics_ok:
                leg.AddEntry(g, lab, "p")
    leg.Draw()
    save(c, "selection")


def fig_floor_intervals(grows):
    """Forest plot: the theory floor of each model of the physical set with the
    68 % (thick) and 95 % (thin) intervals of the inverted chi^2 test, one
    panel per system. Overlapping intervals: the data do not separate the
    models by accuracy."""
    tags = list(dict.fromkeys(r["tag"] for r in grows))
    n = len(tags)
    H = 150 + 24 * n
    c = ROOT.TCanvas("fint", "", 1300, H)
    c.Divide(2, 1, 0.002, 0.0)
    for ip, (sname, title) in enumerate((("ccbar", "c#bar{c}"), ("bbar", "b#bar{b}"))):
        pad = c.cd(ip + 1)
        pad.SetLogx()
        pad.SetGridx()
        pad.SetLeftMargin(0.42 if ip == 0 else 0.04)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(50.0 / H)
        pad.SetBottomMargin(90.0 / H)
        rs = {r["tag"]: r for r in grows if r["system"] == sname}
        lo = min(r["ci95_lo"] for r in rs.values()) * 0.8
        hi = max(r["ci95_hi"] for r in rs.values()) * 1.25
        h = keep(ROOT.TH2F(f"hfi{ip}", "", 10, max(lo, 0.3), hi, n, 0, n))
        h.GetXaxis().SetTitle(f"{title}: #sigma_{{th}} (MeV), 68 % and 95 % CL")
        h.GetXaxis().SetMoreLogLabels()
        h.GetXaxis().SetNoExponent()
        for k, t in enumerate(tags):
            h.GetYaxis().SetBinLabel(n - k, describe(t) if ip == 0 else "")
        for ax in (h.GetXaxis(), h.GetYaxis()):
            ax.SetLabelFont(43)
            ax.SetTitleFont(43)
        h.GetYaxis().SetLabelSize(15)
        h.GetYaxis().SetTickLength(0)
        h.GetXaxis().SetLabelSize(18)
        h.GetXaxis().SetTitleSize(20)
        h.GetXaxis().SetTitleOffset(1.8)
        h.Draw()
        for key, width in (("95", 1), ("68", 5)):
            g = keep(ROOT.TGraphAsymmErrors())
            for k, t in enumerate(tags):
                r = rs.get(t)
                if not r:
                    continue
                j = g.GetN()
                g.SetPoint(j, r["floor"], n - k - 0.5)
                g.SetPointError(j, r["floor"] - r[f"ci{key}_lo"],
                                r[f"ci{key}_hi"] - r["floor"], 0, 0)
            g.SetLineColor(COL[ip])
            g.SetLineWidth(width)
            g.SetMarkerStyle(20)
            g.SetMarkerSize(0.9)
            g.SetMarkerColor(INK)
            g.Draw("PZ same")
    c.cd(0)
    save(c, "floor_intervals")


def fig_wavefunctions(run, title):
    """Reduced radial wave functions u(r) = r R(r) of the S and P states,
    rebuilt from the logged expansion (nu, N_nL, c) alone, with their 68 %
    parametric band (40 draws from the fit's covariance): one panel per
    system, r in fm."""
    import validity
    c = ROOT.TCanvas("wf", "", 1300, 560)
    c.Divide(2, 1, 0.004, 0.0)
    hbarc = 0.1973269804
    for ip, sname in enumerate(("ccbar", "bbar")):
        s = system(run, sname)
        if s is None or "basis" not in s:
            continue
        pad = c.cd(ip + 1)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.10)
        pad.SetBottomMargin(0.13)
        nu = np.array(s["basis"]["nu_GeV2"])
        norms = {int(k): np.array(v) for k, v in s["basis"]["norm_by_L"].items()}
        pick = [x for x in s["states"] if x.get("wavefunction") and x["channel"] in
                ("3S1", "3P2") and x["in_fit"] or x["label"] in ("psi(4040)", "Upsilon(4S)")]
        rmax = 2.0 if sname == "ccbar" else 1.4
        h = keep(ROOT.TH1F(f"hwf{ip}", "", 10, 0, rmax))
        h.SetMinimum(-0.9)
        h.SetMaximum(1.6)
        h.GetXaxis().SetTitle("r (fm)")
        h.GetYaxis().SetTitle("u(r) = r R(r)  (GeV^{1/2})")
        h.Draw("AXIS")
        zero = keep(ROOT.TLine(0, 0, rmax, 0))
        zero.SetLineColor(INK2)
        zero.Draw()
        leg = keep(ROOT.TLegend(0.55, 0.62, 0.97, 0.89))
        leg.SetHeader(("charmonium" if sname == "ccbar" else "bottomonium") + ", " + title)
        r = np.linspace(1e-3, rmax / hbarc, 600)
        for j, x in enumerate(pick[:6]):
            dom = max(x["wavefunction"]["components"], key=lambda q: q["weight"])
            u = r * validity.radial(nu, norms[dom["L"]], dom["c"], dom["L"], r)
            sgn = np.sign(u[np.argmax(np.abs(u) > 1e-3 * np.abs(u).max())])  # u > 0 at small r
            u *= sgn
            b = x.get("wavefunction_band")
            if b:  # 16-84 % of u(r) over parameter draws
                rb = np.array(b["r_fm"])
                ok = rb <= rmax
                lo, hi = sgn * np.array(b["u_q16"])[ok], sgn * np.array(b["u_q84"])[ok]
                lo, hi = np.minimum(lo, hi), np.maximum(lo, hi)
                xs = np.concatenate([rb[ok], rb[ok][::-1]])
                ys = np.concatenate([lo, hi[::-1]])
                band = keep(ROOT.TGraph(len(xs), xs.astype(float), ys.astype(float)))
                band.SetFillColorAlpha(COL[j % len(COL)], 0.30)
                band.SetLineWidth(0)
                band.Draw("F same")
            g = keep(ROOT.TGraph(len(r), (r * hbarc).astype(float), u.astype(float)))
            g.SetLineColor(COL[j % len(COL)])
            g.SetLineWidth(2)
            g.Draw("L same")
            leg.AddEntry(g, tlabel(x["label"]), "l")
        leg.Draw()
    c.cd(0)
    save(c, "wavefunctions")


def root_file(runs, rows):
    f = ROOT.TFile(str(RES / "qgem_sweep.root"), "RECREATE")
    ft = ROOT.TTree("fits", "one entry per configuration and system")
    st = ROOT.TTree("states", "one entry per configuration and PDG state")
    import array
    tag = ROOT.std.string()
    sysn = ROOT.std.string()
    lab = ROOT.std.string()
    chn = ROOT.std.string()
    d = {k: array.array("d", [0.0]) for k in
         ("floor", "chi2", "cond", "loo_rms", "pred_rms", "pdg", "model",
          "sig_par", "sig_num", "pull")}
    i = {k: array.array("i", [0]) for k in ("N", "p", "in_fit", "n")}
    ft.Branch("tag", tag)
    ft.Branch("system", sysn)
    for k in ("floor", "chi2", "cond", "loo_rms", "pred_rms"):
        ft.Branch(k, d[k], f"{k}/D")
    for k in ("N", "p"):
        ft.Branch(k, i[k], f"{k}/I")
    st.Branch("tag", tag)
    st.Branch("system", sysn)
    st.Branch("label", lab)
    st.Branch("channel", chn)
    for k in ("pdg", "model", "sig_par", "sig_num", "pull"):
        st.Branch(k, d[k], f"{k}/D")
    for k in ("in_fit", "n"):
        st.Branch(k, i[k], f"{k}/I")
    for r in rows:
        run = runs[r["tag"]]
        for sname in ("ccbar", "bbar", "bc"):
            s = system(run, sname)
            if s is None:
                continue
            tag.assign(r["tag"])
            sysn.assign(sname)
            if sname != "bc":
                d["floor"][0] = r[f"{sname}_floor"]
                d["chi2"][0] = r[f"{sname}_chi2"]
                d["cond"][0] = r[f"{sname}_cond"]
                d["loo_rms"][0] = r[f"{sname}_loo_rms"]
                d["pred_rms"][0] = r[f"{sname}_pred_rms"]
                i["N"][0] = r[f"{sname}_N"]
                i["p"][0] = r[f"{sname}_p"]
                ft.Fill()
            for x in s["states"]:
                lab.assign(x["label"])
                chn.assign(x["channel"])
                d["pdg"][0] = x["pdg_GeV"]
                nan = float("nan")
                d["model"][0] = nan if x["model_GeV"] is None else x["model_GeV"]
                d["sig_par"][0] = x["sigma_par_GeV"] if x["sigma_par_GeV"] is not None else nan
                d["sig_num"][0] = x["sigma_num_GeV"] if x["sigma_num_GeV"] is not None else nan
                d["pull"][0] = x["pull"] if x["pull"] is not None else 0.0
                i["in_fit"][0] = x["in_fit"]
                i["n"][0] = x["n"]
                st.Fill()
    ft.Write()
    st.Write()
    f.Close()


def fig_limits(rows, tags, r_break=None):
    """Out-of-sample |model - PDG| against the size and the velocity of each
    state: where the potential description stops working. Filled: fitted
    (exact leave-one-out miss); open: predicted."""
    c = ROOT.TCanvas("lim", "", 1300, 560)
    c.Divide(2, 1, 0.004, 0.0)
    for ip, (key, title, lo, hi) in enumerate(
            (("r_rms_fm", "r_{rms} (fm)", 0.0, 1.8),
             ("v2", "#LTv^{2}#GT = #LTp^{2}#GT/m_{Q}^{2}", 0.0, 0.4))):
        pad = c.cd(ip + 1)
        pad.SetLogy()
        pad.SetGridy()
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.14)
        pad.SetBottomMargin(0.13)
        h = keep(ROOT.TH1F(f"hlim{ip}", "", 10, lo, hi))
        h.SetMinimum(0.05)
        h.SetMaximum(400)
        h.GetXaxis().SetTitle(title)
        h.GetYaxis().SetTitle("out-of-sample |M_{model} #minus M_{PDG}| (MeV)")
        h.GetYaxis().SetMoreLogLabels()
        h.GetYaxis().SetNoExponent()
        h.Draw("AXIS")
        leg = keep(ROOT.TLegend(0.13, 0.87, 0.97, 0.995))
        leg.SetNColumns(2)
        leg.SetTextSize(0.032)
        for j, t in enumerate(tags):
            for fitted, mk in ((1, 20), (0, 24)):
                g = keep(ROOT.TGraphAsymmErrors())
                for r in rows:
                    if r["tag"] == t and r["in_fit"] == fitted:
                        y = max(abs(r["resid_oos_MeV"]), 0.06)
                        n = g.GetN()
                        g.SetPoint(n, r[key], y)
                        e = r.get("err_MeV", 0.0)
                        g.SetPointError(n, 0, 0, min(e, y - 0.05), e)
                g.SetMarkerStyle(mk)
                g.SetMarkerColor(COL[j])
                g.SetLineColor(COL[j])
                g.SetMarkerSize(1.1)
                if g.GetN():
                    g.Draw("PZ same")
                if fitted:
                    leg.AddEntry(g, describe(t), "p")
        if key == "r_rms_fm" and r_break:
            ln = keep(ROOT.TLine(r_break, 0.05, r_break, 400))
            ln.SetLineStyle(2)
            ln.SetLineColor(INK2)
            ln.Draw()
            text(r_break + 0.02, 0.1, "string breaking (lattice)", 0.032, INK2, 11,
                 ndc=False)
        leg.Draw()
    c.cd(0)
    save(c, "model_limits")
