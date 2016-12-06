import ROOT

#-----------------------------------------------------------------------------
# Load necessary shared libraries
#
def setPlotDefaults(root, options = None):

    #root.gROOT.SetStyle('Plain')

    root.gStyle.SetFillColor(10)
    root.gStyle.SetFrameFillColor(10)
    root.gStyle.SetCanvasColor(10)
    root.gStyle.SetPadColor(10)
    root.gStyle.SetTitleFillColor(0)
    root.gStyle.SetStatColor(10)
    
    root.gStyle.SetCanvasBorderMode(0)
    root.gStyle.SetFrameBorderMode(0) 
    root.gStyle.SetPadBorderMode(0)   
    root.gStyle.SetDrawBorder(0)      
    root.gStyle.SetTitleBorderSize(0)
    
    root.gStyle.SetFuncWidth(2)
    root.gStyle.SetHistLineWidth(2)
    root.gStyle.SetFuncColor(2)
    
    root.gStyle.SetPadTopMargin(0.08)
    root.gStyle.SetPadBottomMargin(0.16)
    root.gStyle.SetPadLeftMargin(0.16)
    root.gStyle.SetPadRightMargin(0.12)
  
    # set axis ticks on top and right
    root.gStyle.SetPadTickX(1)         
    root.gStyle.SetPadTickY(1)         
  
    # Set the background color to white
    root.gStyle.SetFillColor(10)           
    root.gStyle.SetFrameFillColor(10)      
    root.gStyle.SetCanvasColor(10)         
    root.gStyle.SetPadColor(10)            
    root.gStyle.SetTitleFillColor(0)       
    root.gStyle.SetStatColor(10)           
  
  
    # Turn off all borders
    root.gStyle.SetCanvasBorderMode(0)
    root.gStyle.SetFrameBorderMode(0) 
    root.gStyle.SetPadBorderMode(0)   
    root.gStyle.SetDrawBorder(0)      
    root.gStyle.SetTitleBorderSize(0) 
  
    # Set the size of the default canvas
    root.gStyle.SetCanvasDefH(400)          
    root.gStyle.SetCanvasDefW(650)          
    #gStyle->SetCanvasDefX(10)
    #gStyle->SetCanvasDefY(10)   
  
    # Set fonts
    font = 42
    #root.gStyle.SetLabelFont(font,'xyz')
    #root.gStyle.SetStatFont(font)       
    #root.gStyle.SetTitleFont(font)      
    #root.gStyle.SetTitleFont(font,'xyz')
    #root.gStyle.SetTextFont(font)       
    #root.gStyle.SetTitleX(0.3)        
    #root.gStyle.SetTitleW(0.4)        
  
   # Set Line Widths
   #gStyle->SetFrameLineWidth(0)
   #root.gStyle.SetFuncWidth(2)
   #root.gStyle.SetHistLineWidth(2)
   #root.gStyle.SetFuncColor(2)
   #
   # Set tick marks and turn off grids
    root.gStyle.SetNdivisions(505,'xyz')
   #
   # Set Data/Stat/... and other options
   #root.gStyle.SetOptDate(0)
   #root.gStyle.SetOptDate(0)   
   #root.gStyle.SetDateX(0.1)
   #root.gStyle.SetDateY(0.1)
   #gStyle->SetOptFile(0)
   ##root.gStyle.SetOptStat(1110)
    root.gStyle.SetOptStat(0)
    #root.gStyle.SetOptFit(111)
    root.gStyle.SetStatFormat('4.3f')
    root.gStyle.SetFitFormat('4.3f')
   #gStyle->SetStatTextColor(1)
   #gStyle->SetStatColor(1)
   #gStyle->SetOptFit(1)
   #gStyle->SetStatH(0.20)
   #gStyle->SetStatStyle(0)
   #gStyle->SetStatW(0.30)
   #gStyle -SetStatLineColor(0)
   #root.gStyle.SetStatX(0.919)
   #root.gStyle.SetStatY(0.919)
   #root.gStyle.SetOptTitle(0)
   #gStyle->SetStatStyle(0000)    # transparent mode of Stats PaveLabel
   #root.gStyle.SetStatBorderSize(0)
   #
    #root.gStyle.SetLabelSize(0.065,'xyz')
    #gStyle -> SetLabelOffset(0.005,'xyz')
    #root.gStyle.SetTitleY(.5)
    root.gStyle.SetTitleOffset(1.0,'xz')
    root.gStyle.SetTitleOffset(1.1,'y')
    root.gStyle.SetTitleSize(0.065, 'xyz')
    root.gStyle.SetLabelSize(0.065, 'xyz')
    #root.gStyle.SetTextAlign(22)
    root.gStyle.SetTextSize(0.1)
   #
   ##root.gStyle.SetPaperSize(root.TStyle.kA4)  
    root.gStyle.SetPalette(1)
   #
   ##root.gStyle.SetHistMinimumZero(True)
   
    root.gROOT.ForceStyle()
#-----------------------------------------
def Format(h):

    h.SetLineColor(1)
    h.SetMarkerColor(1)
    h.SetTitle('Projection of 2D Efficiency')
#-----------------------------------------
def Style():
    ROOT.gROOT.LoadMacro('/Users/schae/testarea/CAFAna/HWWMVACode/atlasstyle-00-03-05/AtlasStyle.C')                   
    ROOT.gROOT.LoadMacro('/Users/schae/testarea/CAFAna/HWWMVACode/atlasstyle-00-03-05/AtlasUtils.C')
    ROOT.SetAtlasStyle()

setPlotDefaults(ROOT)
name='effx'
xmin=-5.0
xmax=45.0

f=ROOT.TFile.Open(name+'.root')
can = f.Get('c1_n2')

can.Draw()
color=2
h=None
leg=ROOT.TLegend(0.2,0.2,0.6,0.4)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.065)
stat=None
syststat=None
for j in can.GetListOfPrimitives():
    print j.GetName()
    if j.InheritsFrom("TH1D"):
        j.SetStats(0)
        print 'Hist:',j.GetName()
        h=j.Clone()
    if j.InheritsFrom("TH1F"):
        j.SetStats(0)
        print 'HistF:',j.GetName()
    if j.InheritsFrom("TH1"):
        j.SetStats(0)
        print 'His:',j.GetName()
    if j.InheritsFrom("TGraphAsymmErrors"):
        #j.SetStats(0)
        print 'err: ',j.GetName()
        j.SetLineColor(2)
    if j.InheritsFrom("TGraph"):
        #j.SetStats(0)
        j.GetXaxis().SetRangeUser(xmin,xmax)
        j.GetYaxis().SetRangeUser(0.75,1.01)        
        print 'reg',j.GetName()
        j.SetLineColor(color)
        j.SetMarkerColor(color)        
        if color==1:
            h=j.Clone()
            stat=j.Clone()
            leg.AddEntry(j,'Syst. Unc.')
        else:
            j.SetFillStyle(3002)
            j.SetFillColor(1)
            j.SetMarkerColor(1)            
            syststat=j.Clone()
            syststat.SetFillStyle(3002)
            syststat.SetFillColor(1)
            j.SetLineWidth(0)
            leg.AddEntry(j,'Stat.+Syst. Unc.')            
        color-=1
        
            
line = ROOT.TLine(xmin,1.0,xmax,1.0)
line.SetLineWidth(3)
line.SetLineColor(1)
line.Draw()


lineelec = ROOT.TLine(xmin,0.95,xmin,1.01)
lineelec.SetLineWidth(4)
lineelec.SetLineColor(2)
lineelec.Draw()

lineelec2 = ROOT.TLine(xmax,0.95,xmax,1.01)
lineelec2.SetLineWidth(4)
lineelec2.SetLineColor(2)
lineelec2.Draw()
l = ROOT.TLatex(0.12, 0.95, 'Pixel Center')
l.SetNDC()
l.SetTextFont(72)
l.SetTextSize(0.045)
l.SetTextColor(2)
l.Draw()

l2 = ROOT.TLatex(0.8, 0.95, 'Pixel Center')
l2.SetNDC()
l2.SetTextFont(72)
l2.SetTextSize(0.045)
l2.SetTextColor(2)
l2.Draw()

linecent = ROOT.TLine((xmax-xmin)/2.0+xmin,0.95,(xmax-xmin)/2.0+xmin,1.01)
linecent.SetLineWidth(4)
linecent.SetLineColor(4)
linecent.SetLineStyle(2)
linecent.Draw()
lc = ROOT.TLatex(0.46, 0.95, 'Pixel Edge')
lc.SetNDC()
lc.SetTextFont(72)
lc.SetTextSize(0.045)
lc.SetTextColor(4)
lc.Draw()

syststat.Draw('e2 same')
stat.Draw(' same')
syststat
leg.Draw()

can.Update()
can.WaitPrimitive()

can.SaveAs(name+'_format'+'.eps')
can.SaveAs(name+'_format'+'.pdf')
can.SaveAs(name+'_format'+'.C')
can.SaveAs(name+'_format'+'.root')
f.Close()
