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
def Style():
    ROOT.gROOT.LoadMacro('/Users/schae/testarea/CAFAna/HWWMVACode/atlasstyle-00-03-05/AtlasStyle.C')                   
    ROOT.gROOT.LoadMacro('/Users/schae/testarea/CAFAna/HWWMVACode/atlasstyle-00-03-05/AtlasUtils.C')
    ROOT.SetAtlasStyle()

def Fit(file_name='',_suffix=''):
    rebin=2
    #f = ROOT.TFile.Open('./../Judith_original/TestData/dut-60V_runs-2-3-4-5-6-7-1-4-5-6_settings1_sync_analysis-result.root')
    #f = ROOT.TFile.Open('./../Judith_original/TestData/dut-120V_runs-23-24-25-26-27-28-29-30-1-2-3-4-5-6-7_settings1_sync_analysis-result-cutslope.root')
    #f = ROOT.TFile.Open('./../Judith_original/TestData/dut-90V_runs-9-11-1-2-3-5-6-8-9_settings_sync_analysis-result.root')
    #f = ROOT.TFile.Open('./../Judith_original/TestData/dut-120V_runs-23-24-25-26-27-28-29-30-1-2-3-4-5-6-7_settings1_sync_analysis-result.root')
    #f = ROOT.TFile.Open('./../Judith_original/TestData/dut-90V_runs-9-11-1-2-3-5-6-8-9_settings_sync_analysis-result-nocut-align.root')
    f = ROOT.TFile.Open(file_name)
    #dut-60V_runs-2-3-4-5-6-7-1-4-5-6_settings1_sync_analysis-result-cutt0-align.root
    #dut-120V_runs-23-24-25-26-27-28-29-30-1-2-3-4-5-6-7_settings1_sync_analysis-result-nocut-align.root
    #f = ROOT.TFile.Open('./../Judith_original/TestData/dut-120V_runs-23-24-25-26-27-28-29-30-1-2-3-4-5-6-7_settings1_sync_analysis-result.root')
    #twoD = f.Get('Efficiency/sensor0_TrackResEffFine')
    twoD = f.Get('Efficiency/DUTPlane0TrackResidualHitFine')
    twoDall = f.Get('Efficiency/DUTPlane0TrackResidualFine')        
    
    ring = twoD.Clone()

    can = ROOT.TCanvas("c2","c2",100,10,800,600);
    can.cd()
    twoD.Draw('colz')
    can.Update()
    can.WaitPrimitive() 
    t1 = raw_input('input the x shift. typically the x-mean: ')
    t2 = raw_input('input the y shift. typically the y-mean: ')
    
    xshift = float(t1)
    yshift = float(t2)
    bin_width=5
    #bin_width=15
    midX=twoD.GetNbinsX()/2+int(xshift)
    midY=twoD.GetNbinsY()/2+int(yshift)
    midX_noshift=twoD.GetNbinsX()/2
    midY_noshift=twoD.GetNbinsY()/2
    #for q in range(0,int(midX/float(bin_width))):
    ##for q in range(0,midX):
    #    #i=q*bin_width
    #    #sum up and average
    #    mysum = 0.0
    #    entries = 0.0
    #    for i in range(q*bin_width,(q+1)*bin_width):
    #        for j in range(-i,i+1):
    #            mysum+=twoD.GetBinContent(midX+i,midY+j)
    #            mysum+=twoD.GetBinContent(midX-i,midY+j)
    #            entries+=twoDall.GetBinContent(midX+i,midY+j)
    #            entries+=twoDall.GetBinContent(midX-i,midY+j)                
    #            if j!=i and j!=-1:
    #                mysum+=twoD.GetBinContent(midX+j,midY+i)
    #                mysum+=twoD.GetBinContent(midX+j,midY-i)
    #                entries+=twoDall.GetBinContent(midX+j,midY+i)
    #                entries+=twoDall.GetBinContent(midX+j,midY-i)
    #            #else:
    #            #    entries-=1
    #
    #    if bin_width!=1 and q==1:
    #        mysum+=twoD.GetBinContent(midX,midY)
    #        entries+=twoDall.GetBinContent(midX,midY)            
    #        #entries+=1
    #    if entries==0.0:
    #        entries+=1.0
    #    #print 'entries: ',entries
    #    for i in range(q*bin_width,(q+1)*bin_width):
    #        for j in range(-i,i+1):
    #            ring.SetBinContent(midX_noshift+i,midY_noshift+j,mysum/entries/0.97)
    #            ring.SetBinContent(midX_noshift-i,midY_noshift+j,mysum/entries/0.97)
    #            #if j!=i and j!=-1:
    #            ring.SetBinContent(midX_noshift+j,midY_noshift+i,mysum/entries/0.97)
    #            ring.SetBinContent(midX_noshift+j,midY_noshift-i,mysum/entries/0.97)
    #if bin_width!=1:
    #    ring.SetBinContent(midX_noshift,midY_noshift,ring.GetBinContent(midX_noshift+1,midY_noshift)) 
    #
    #for q in range(0,int(midX/float(bin_width))):
    #    if ring.GetBinContent(midX+q*bin_width+1,midY)>0.0:
    #        print 'Distance from center: ',q*bin_width+1,' Eff: ',ring.GetBinContent(midX+q*bin_width+1,midY)
    print "Signal Hits: ",ring.Integral(),' ',ring.Integral(0,10001,0,10001)
    print "Telescope Hits: ",twoDall.Integral()
    rebinF=1
    ring.RebinX(rebinF)
    ring.RebinY(rebinF)
    twoDall.RebinX(rebinF)
    twoDall.RebinY(rebinF)
        
    ring.GetXaxis().SetRangeUser(-200,100)
    ring.GetYaxis().SetRangeUser(-200,1000)
    ring.Draw('colz')
    can.Update()
    can.WaitPrimitive()
            
    twoDall.GetXaxis().SetRangeUser(-200,100)
    twoDall.GetYaxis().SetRangeUser(-200,1000)
    twoDall.Draw('colz')
    can.Update()
    can.WaitPrimitive()
    
    ring.Divide(twoDall)
    ring.GetXaxis().SetRangeUser(-200,100)
    ring.GetYaxis().SetRangeUser(-200,1000)
    ring.GetXaxis().SetRangeUser(-200,200)
    ring.GetYaxis().SetRangeUser(-200.0,100)
    ring.SetTitle('DUT Hit Efficiency')      
    ring.GetZaxis().SetRangeUser(0.0,1.0)
    ring.SetTitle('')
    ring.SetStats(0)
    ring.Draw('colz')
    can.Update()
    can.WaitPrimitive()
    raw_input('Waiting for you to finish editting')
    can.SaveAs('efficiencyQuick_'+_suffix+'.eps')
    can.SaveAs('efficiencyQuick_'+_suffix+'.pdf')
    can.SaveAs('efficiencyQuick_'+_suffix+'.C')
    can.SaveAs('efficiencyQuick_'+_suffix+'.root')

#Style()
setPlotDefaults(ROOT)
#Fit('60V')
#Fit('./../Judith_original/TestData/dut-60V_runs-2-3-4-5-6-7-1-4-5-6_settings1_sync_analysis-result-cutt0wider-align.root','cutt0_60V')
#Fit('./../Judith_original/TestData/dut-90V_runs-9-11-1-2-3-5-6-8-9_settings_sync_analysis-result-cutt0wider-align.root','cutt0_90V')
#Fit('./../Judith_original/TestData/dut-120V_runs-23-24-25-26-27-28-29-30-1-2-3-4-5-6-7_settings1_sync_analysis-result-cutt0wider-align.root','cutt0_120V')

#Fit('./../Judith_original/TestData/dut-60V_runs-2-3-4-5-6-7-1-4-5-6_settings1_sync_analysis-result-nocut-align.root','nocut_60V')
#Fit('./../Judith_original/TestData/dut-90V_runs-9-11-1-2-3-5-6-8-9_settings_sync_analysis-result-nocut-align.root','nocut_90V')
#Fit('./../Judith_original/TestData/dut-120V_runs-23-24-25-26-27-28-29-30-1-2-3-4-5-6-7_settings1_sync_analysis-result-nocut-align.root','nocut_120V')
#
#Fit('./../Judith_original/TestData/dut-60V_runs-2-3-4-5-6-7-1-4-5-6_settings1_sync_analysis-result-cutslope-align.root','cutslope_60V')
#Fit('./../Judith_original/TestData/dut-90V_runs-9-11-1-2-3-5-6-8-9_settings_sync_analysis-result-cutslope-align.root','cutslope_90V')
#Fit('./../Judith_original/TestData/dut-120V_runs-23-24-25-26-27-28-29-30-1-2-3-4-5-6-7_settings1_sync_analysis-result-cutslope-align.root','cutslope_120V')
#Fit('TowerJazz/dut_run804_sync-analysis-result.root','No Cut')
#Fit('TowerJazz/run_804_tj_W3R15_50um_6V_new_sync-analysis-result.root','No Cut')
#Fit('TowerJazz/run_868_tj_W3R13_50um_6V_1DRS_new_sync-analysis-result.root','No Cut')
#Fit('TowerJazz/run_804_tj_W3R15_50um_6V_fft4_sync-analysis-result.root','UnIrr')
#Fit('TowerJazz/run_32880_tj_W3R13_50um_6V_3DRS-analysis-result.root','UnIrr')
#Fit('TowerJazz/run_32877_tj_W3R13_50um_6V_3DRS_sync-analysis-result.root','UnIrr')
Fit('TowerJazz/run_32878_tj_W3R13_50um_6V_3DRS_sync-analysis-result.root','UnIrr')
#Fit('TowerJazz/run_868_tj_W3R13_50um_6V_1DRS_fft4_sync-analysis-result.root','Irr')
#Fit('TowerJazz/run804_new5_sync-analysis-result.root','No_Cut')
#Fit('TowerJazz/tmp_804_sync-analysis-result.root','No Cut')
#Fit('TowerJazz/test_lowpass_sync-analysis-result.root','No Cut')
