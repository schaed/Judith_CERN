import ROOT
import math
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
   #root.gStyle.SetDateX(0.1)
   #root.gStyle.SetDateY(0.1)
   #gStyle->SetOptFile(0)
   ##root.gStyle.SetOptStat(1110)
    root.gStyle.SetOptStat(1111)
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
#def Style():
    #ROOT.gROOT.LoadMacro('/Users/schae/testarea/CAFAna/HWWMVACode/atlasstyle-00-03-05/AtlasStyle.C')                   
    #ROOT.gROOT.LoadMacro('/Users/schae/testarea/CAFAna/HWWMVACode/atlasstyle-00-03-05/AtlasUtils.C')
    #ROOT.SetAtlasStyle()

def ScaleX(x,shift):
    v = x + shift 
    return v

def ScaleY(y,shift):
    v = y + shift 
    return v

def ScaleAxis(a, Scale,shift):

  if (not a): 
      return
  if a.GetXbins().GetSize():
      X = ROOT.TArrayD(a.GetXbins())
      for i in range(0, X.GetSize()): 
          X[i] = Scale(X[i],shift)
      a.Set((X.GetSize() - 1), X.GetArray()) 
  else: 
      a.Set( a.GetNbins(),
              Scale(a.GetXmin(),shift), 
              Scale(a.GetXmax() ,shift))
  return


def ScaleXaxis(h, Scale,shift):
  if (not h): 
      return
  ScaleAxis(h.GetXaxis(), Scale,shift)
  return


def ScaleYaxis(h,Scale,shift):
  if (not h): 
      return
  ScaleAxis(h.GetYaxis(), Scale,shift)
  return

def corrEff(ring,ringup,ringun,hmap,hmapErr,shiftmapx,shiftmapy,rangemapx,rangemapy):
    for i in range (0,rangemapx):
      for j in range (0,rangemapy):
        
        if ring.GetBinContent(i+1+shiftmapx,j+1+shiftmapy) != 0:
           deo = ring.GetBinError(i+1+shiftmapx,j+1+shiftmapy)/ring.GetBinContent(i+1+shiftmapx,j+1+shiftmapy)
        else:
            deo = 0
	if hmap.GetBinContent(i+1,j+1) <0.5:
	  result = 0
	  resultup = 0
	  resultun = 0
          dem = 0
	else:
          dem = hmapErr.GetBinContent(i+1,j+1)/hmap.GetBinContent(i+1,j+1)
	  result = ring.GetBinContent(i+1+shiftmapx,j+1+shiftmapy)/(hmap.GetBinContent(i+1,j+1))
	  resultup = ring.GetBinContent(i+1+shiftmapx,j+1+shiftmapy)/(hmap.GetBinContent(i+1,j+1)*(1+hmapErr.GetBinContent(i+1,j+1)))
	  resultun = ring.GetBinContent(i+1+shiftmapx,j+1+shiftmapy)/(hmap.GetBinContent(i+1,j+1)*(1-hmapErr.GetBinContent(i+1,j+1)))
	ring.SetBinContent(i+1+shiftmapx,j+1+shiftmapy,result)
        ring.SetBinError(i+1+shiftmapx,j+1+shiftmapy,math.sqrt(pow(deo,2) + pow(dem,2)))
	ringup.SetBinContent(i+1+shiftmapx,j+1+shiftmapy,resultup)
	ringun.SetBinContent(i+1+shiftmapx,j+1+shiftmapy,resultun)

def corrEff_red(ring,hmap,shiftmapx,shiftmapy,rangemapx,rangemapy):
    for i in range (0,rangemapx):
      for j in range (0,rangemapy):
	if hmap.GetBinContent(i+1,j+1) <0.1:
	  result = 0
	else:
	  result = ring.GetBinContent(i+1+shiftmapx,j+1+shiftmapy)/(hmap.GetBinContent(i+1,j+1))
	ring.SetBinContent(i+1+shiftmapx,j+1+shiftmapy,result)
        
    
def Fit(file_name='',map1='',_suffix=''):

    f = ROOT.TFile.Open(file_name)
    map1 = ROOT.TFile.Open(map1)
    twoD = f.Get('Efficiency/DUTPlane0TrackResidualHitFine')
    twoDall = f.Get('Efficiency/DUTPlane0TrackResidualFine') 
    m1 = map1.Get('cMap')
    hmap = m1.GetPrimitive("map")
    m2 = map1.Get('cMapErr')
    hmapErr = m2.GetPrimitive("mapErr")
    
    can = ROOT.TCanvas("c2","c2",100,10,800,600);
    can.cd()
    twoD_prev = twoD.Clone()
    twoD_prev.RebinX(5)
    twoD_prev.RebinY(5)
    twoD_prev.Draw('colz')
    can.Update()
    t1 = int(raw_input('input the x shift. typically the x-mean: '))
    t2 = int(raw_input('input the y shift. typically the y-mean: '))
    pixwidth = int(raw_input('input the pixel width: '))
    
    xshift = t1
    yshift = t2

    ring = twoD.Clone()
    ringall= twoDall.Clone()

    for i in range (0,800):
      for j in range (0,800):
	result = twoD.GetBinContent(int(i+(xshift)),int(j+(yshift)))
	result1 = twoDall.GetBinContent(int(i+(xshift)),int(j+(yshift)))
	ring.SetBinContent(i,j,result)
	ringall.SetBinContent(i,j,result1)
        

    # ring.Draw("colz")
    # can.Update()
    # can.WaitPrimitive()
    # ringall.Draw("colz")
    # can.Update()
    # can.WaitPrimitive()
    #ScaleXaxis(twoD, ScaleX, xshift)
    #ScaleXaxis(twoDall, ScaleX, xshift)
    #ScaleYaxis(twoD, ScaleY, yshift)
    #ScaleYaxis(twoDall, ScaleY, yshift)
    
    ring1= ring.Clone()
    ringup = ring.Clone() 
    ringun = ring.Clone()
    
    posmin = -90
    posmax = 130
    poscorrmin = -40
    poscorrmax = 80

    rangemapx = hmap.GetNbinsX()
    rangemapy = hmap.GetNbinsY()
    shiftmapx = int(abs(twoD.GetXaxis().GetXmin()))-pixwidth
    shiftmapy = int(abs(twoD.GetYaxis().GetXmin()))-pixwidth
    
    hmap.Sumw2()
    corrEff(ring,ringup,ringun,hmap,hmapErr,shiftmapx,shiftmapy,rangemapx,rangemapy)
    
    print "Signal Hits: ",ring.Integral(),' ',ring.Integral(0,10001,0,10001)
    print "Telescope Hits: ",twoDall.Integral()
    
    rebinF=5
    ring.RebinX(rebinF)
    ring.RebinY(rebinF)
    ring1.RebinX(rebinF)
    ring1.RebinY(rebinF)
    ringup.RebinX(rebinF)
    ringup.RebinY(rebinF)
    ringun.RebinX(rebinF)
    ringun.RebinY(rebinF)
    ringall.RebinX(rebinF)
    ringall.RebinY(rebinF)
    
    # can.cd()
    # ring1.Draw("colz")
    # can.Update()
    # can.WaitPrimitive()
    # raw_input("Wait")
    
    ring.Sumw2()
    ring1.Sumw2()
    ringall.Sumw2()
    
    ring.Divide(ring,ringall,1.0,1.0,"B")
    ring1.Divide(ring1,ringall,1.0,1.0,"B")
    ringup.Divide(ringall)
    ringun.Divide(ringall)
    ringpr = ring.Clone()
    ringpr1 = ring1.Clone()
    ringprerr = ringup.Clone()
    ringprall = ringall.Clone()

    can.cd()
    ring1.GetXaxis().SetRangeUser(posmin,posmax)
    ring1.GetYaxis().SetRangeUser(posmin,posmax)
    ring1.GetZaxis().SetRangeUser(0.00,1.00)
    ring1.SetTitle('')
    ring1.SetStats(0)
    ring1.Draw('colz')
    can.Update()
    
    poscorrmin = -30
    poscorrmax = 85
    
    can3 = ROOT.TCanvas("c3","c3",100,10,800,600);
    can3.cd()
    ring.GetXaxis().SetRangeUser(poscorrmin,poscorrmax)
    ring.GetYaxis().SetRangeUser(poscorrmin,poscorrmax)
    ring.GetZaxis().SetRangeUser(0.0,1.0)
    ring.SetTitle('')
    ring.SetStats(0)
    ring.Draw('colz')
    can3.Update()


    xslice = float(raw_input("Start point of xslice: "))
    yslice = float(raw_input("Start point of yslice: "))
    slwidth = float(raw_input("Width of your slice: "))

    raw_input('Waiting for you to finish editing')
    can.SaveAs('efficiencyQuick_'+_suffix+'.root')
    can3.SaveAs('efficiencyQuick_mapcorr_'+_suffix+'.root')

    can4 = ROOT.TCanvas("c4","c4",100,10,800,600);
    can4.cd()
    ringpr1.GetXaxis().SetRangeUser(posmin,posmax)
    ringpr1.GetYaxis().SetRangeUser(yslice,yslice+slwidth)
    ringprall.GetXaxis().SetRangeUser(posmin,posmax)
    ringprall.GetYaxis().SetRangeUser(yslice,yslice+slwidth)
    h1x = ringpr1.ProjectionX()
    h2x = ringprall.ProjectionX()
    h1x.Scale(5/slwidth)
    ex = []
    for i in range(h1x.GetNbinsX()):
       eps = h1x.GetBinContent(i+1)
       N = h2x.GetBinContent(i+1)
       if N !=0:
           error = math.sqrt(eps*abs(1-eps)/N)
       else:
           error = 0
       ex.append(error)
       h1x.SetBinError(i+1,error)
    h1x.SetStats(0)
    h1x.Draw("E")
    can4.Update()
    can5 = ROOT.TCanvas("c5","c5",100,10,800,600);
    can5.cd()
    ringpr1.GetYaxis().SetRangeUser(posmin,posmax)
    ringpr1.GetXaxis().SetRangeUser(xslice,xslice+slwidth)
    ringprall.GetYaxis().SetRangeUser(posmin,posmax)
    ringprall.GetXaxis().SetRangeUser(xslice,xslice+slwidth)
    h1y = ringpr1.ProjectionY()
    h2y = ringprall.ProjectionY()
    h1y.Scale(5/slwidth)
    ey = []
    for i in range(h1y.GetNbinsX()):
       eps = h1y.GetBinContent(i+1)
       N = h2y.GetBinContent(i+1)
       if N !=0:
           error = math.sqrt(eps*abs(1-eps)/N)
       else:
           error = 0
       ey.append(error)
       h1y.SetBinError(i+1,error)
    h1y.SetStats(0)
    h1y.Draw("E")
    can5.Update()
    
    raw_input('Waiting for you to finish editing')
    can4.SaveAs('efficiencyQuick_projectionX_'+_suffix+'.root') 
    can5.SaveAs('efficiencyQuick_projectionY_'+_suffix+'.root')

    can6 = ROOT.TCanvas("c6","c6",100,10,800,600);
    can6.cd()
    ringpr.GetXaxis().SetRangeUser(poscorrmin,poscorrmax)
    ringpr.GetYaxis().SetRangeUser(yslice,yslice+slwidth)
    h1 = ringpr.ProjectionX()
    h1.Scale(5/slwidth)
    for i in range(h1.GetNbinsX()):
        h1.SetBinError(i+1,ex[i+6])
    h1.SetStats(0)
    h1.Draw("E")
    can6.Update()
    can7 = ROOT.TCanvas("c7","c7",100,10,800,600);
    can7.cd()
    ringpr.GetYaxis().SetRangeUser(poscorrmin,poscorrmax)
    ringpr.GetXaxis().SetRangeUser(xslice,xslice+slwidth)
    h1 = ringpr.ProjectionY()
    h1.Scale(5/slwidth)
    for i in range(h1.GetNbinsX()):
        h1.SetBinError(i+1,ey[i+6])
    h1.SetStats(0)
    h1.Draw("E")
    can7.Update()
    
    raw_input('Waiting for you to finish editing')
    can6.SaveAs('efficiencyQuick_projectionX_corr_'+_suffix+'.root')
    can7.SaveAs('efficiencyQuick_projectionY_corr_'+_suffix+'.root')



setPlotDefaults(ROOT)
mapspath = '../TelescopeResSmear/maps/'
Fit('../../TowerJazz/run_33118_tj_W3R15_M129_6V_3DRS_fft12_sync-analysis-result.root',mapspath+'map_clusterType_4_val_pixX_51.000000_pixY_50.000000_resX_8.500000_resY_10.500000_gridX_1.000000_gridY_1.000000_integratorType_0.root','Irr_50')
Fit('../../TowerJazz/run_33118_tj_W3R15_M129_6V_3DRS_fft12_sync-analysis-result.root',mapspath+'map_clusterType_4_val_pixX_51.000000_pixY_50.000000_resX_9.500000_resY_11.500000_gridX_1.000000_gridY_1.000000_integratorType_0.root','Irr_50_err')



