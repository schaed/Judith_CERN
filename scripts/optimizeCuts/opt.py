import ROOT

class Pixel:

    def __init__(self, num=0, cut=''):

        self.num=num
        self.hCharge = ROOT.TH1F('charge%s_%s' %(cut,num),'charge%s_%s' %(cut,num),200,0.0,0.2)
        self.hT0 = ROOT.TH1F('t0%s_%s' %(cut,num),'t0%s_%s'%(cut,num),1000,0.0,1000.0)
        #self.hTime = ROOT.TH1F('timing%s_%s' %(cut,num),'timing%s_%s' %(cut,num),20000,-10000.0,10000.0)
        self.hTime = ROOT.TH1F('timing%s_%s' %(cut,num),'timing%s_%s' %(cut,num),500,0.0,500.0)
        self.hMag = ROOT.TH1F('fftmag%s_%s' %(cut,num),'fftmag%s_%s' %(cut,num),500,0.0,15.0)
        self.hPhs = ROOT.TH1F('fftphs%s_%s' %(cut,num),'fftphs%s_%s' %(cut,num),64,-3.2,3.2) 

        self.hCharge.GetXaxis().SetTitle('Charge [V]')
        self.hCharge.GetYaxis().SetTitle('Events')
        self.hT0.GetXaxis().SetTitle('Time after Trigger [ns]')
        self.hT0.GetYaxis().SetTitle('Events')
        self.hTime.GetXaxis().SetTitle('Charge Collection Time [ns]')
        self.hTime.GetYaxis().SetTitle('Events')
        self.hMag.GetXaxis().SetTitle('FFT Magnitude of low freq. [ns]')
        self.hMag.GetYaxis().SetTitle('Events')
        self.hPhs.GetXaxis().SetTitle('FFT Phase of low freq. [rad]')
        self.hPhs.GetYaxis().SetTitle('Events')
    def Fill(self, ch, t0, timing, mag, phs):
        self.hCharge.Fill(ch)
        self.hT0.Fill(t0)
        self.hTime.Fill(timing)
        self.hMag.Fill(mag)
        self.hPhs.Fill(phs)
    def Write(self,fout):
        self.hCharge.SetDirectory(fout)
        self.hT0.SetDirectory(fout)
        self.hTime.SetDirectory(fout)
        self.hMag.SetDirectory(fout)
        self.hPhs.SetDirectory(fout)
        #fout.Write()
        #fout.Close()
        
#f = ROOT.TFile.Open('../../TowerJazz/run_904_tj_W3R13_50um_6V_1DRS.root')
f = ROOT.TFile.Open('../../TowerJazz/run_804_tj_W3R15_50um_6V.root')
fout = ROOT.TFile.Open('fout_804_01V_60Time.root','RECREATE')

hits = f.Get('Plane0/Hits')
ii=0
pixels=[]
pixels_15V=[]
pixels_15V_time=[]
pixels_15V_t0=[]
for h in hits:

    if ii%10000==0:
        print 'Event: ',ii
    ii+=1        
    #if h.NHits==4:
    #    if h.Value[3]<0.0:
    #        continue
    #print 'Event: ',ii
    if len(pixels)==0 and h.NHits>0:
        for p in range(0,h.NHits):
            pixels+=[Pixel(p)]
            pixels_15V+=[Pixel(p,'_Val15V')]
            pixels_15V_time+=[Pixel(p,'_Val15VTime')]
            pixels_15V_t0+=[Pixel(p,'_Val15VT0')]
    for p in range(0,h.NHits):
        #print 'pixel: ',p,' Charge: ',h.Value[p], h.T0[p], h.Timing[p]
        #print 'Charge: ',h.Value[p]
        # Cuts:
        #if h.Value[p]<0.05:
        #    continue
        #if not (h.Timing[p]<380.0 and h.Timing[p]>370.0):
        #if not (h.Timing[p]<220.0 and h.Timing[p]>180.0):
        #    continue
        #if not (h.Timing[p]<230.0 and h.Timing[p]>218.0):
        #    continue
        #print ' PASS Charge pixel: ',p,' Charge: ',h.Value[p], h.T0[p], h.Timing[p]
        pixels[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])
        if h.Value[p]>0.003:
            pixels_15V[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])
            if h.Timing[p]>15.0 and  h.Timing[p]<60.0:
                pixels_15V_time[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])
            if (h.T0[p]<280.0 and h.T0[p]>200.0) and (h.Timing[p]>15.0 and  h.Timing[p]<60.0): 
                pixels_15V_t0[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])

for p in pixels:
    p.Write(fout)
for p in pixels_15V:
    p.Write(fout)
for p in pixels_15V_time:
    p.Write(fout)
for p in pixels_15V_t0:
    p.Write(fout)    
fout.Write()
fout.Close()
print 'Done'
