import ROOT

def projslice(heff,poscorrmin, poscorrmax, xslice,yslice,slwidth):
         
    can6 = ROOT.TCanvas("c6","c6",100,10,800,600);
    can6.cd()
    heff.GetXaxis().SetRangeUser(poscorrmin,poscorrmax)
    heff.GetYaxis().SetRangeUser(yslice,yslice+slwidth)
    h1 = heff.ProjectionX()
    h1.Scale(5/slwidth)
    h1.SetStats(0)
    h1.Draw("E")
    can6.Update()
    can7 = ROOT.TCanvas("c7","c7",100,10,800,600);
    can7.cd()
    heff.GetYaxis().SetRangeUser(poscorrmin,poscorrmax)
    heff.GetXaxis().SetRangeUser(xslice,xslice+slwidth)
    h1 = heff.ProjectionY()
    h1.Scale(5/slwidth)
    h1.SetStats(0)
    h1.Draw("E")
    can7.Update()
    raw_input('Waiting for you to finish editing')

def proj():
    effcorr= ROOT.TFile.Open("efficiencyQuick_mapcorr_Irr_50.root")
    canvas = effcorr.Get("c3")
    heff = canvas.GetPrimitive("DUTPlane0TrackResidualHitFine")
    
    can = ROOT.TCanvas("can","can",900,600)
    heff.Draw("colz")
    can.Update()
    can.WaitPrimitive()
    
    poscorrmin = -30
    poscorrmax = 85
    
    xslice = float(raw_input("Start point of xslice: "))
    yslice = float(raw_input("Start point of yslice: "))
    slwidth = float(raw_input("Width of your slice: "))
    
    projslice(heff,poscorrmin, poscorrmax, xslice,yslice,slwidth)
    
proj()
