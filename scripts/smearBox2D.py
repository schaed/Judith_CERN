import ROOT

can = ROOT.TCanvas()
f = ROOT.TFile.Open('quick.root','recreate')

hbox=ROOT.TH2F('box','',100,-100.0,100.0,100,-100.0,100.0)
hsmear=ROOT.TH2F('smear','',100,-100.0,100.0,100,-100.0,100.0)


reso = 20.1
reso = 14.1
r = ROOT.TRandom3();
r.SetSeed(5)
NEVT=1000000
w=1.0/float(NEVT)*100.0
for e in range(0,NEVT): #for i in range(25,75):
    a=r.Rndm()*100.0 - 50.0
    b=r.Rndm()*100.0 - 50.0    
    #if round(a)==0:
    #    print a, round(a)
    hbox.Fill(round(a),round(b),w)
    a+=r.Gaus(0.0,reso)
    b+=r.Gaus(0.0,reso)        
    hsmear.Fill(round(a),round(b),w)


    
print 'Reso: ',reso
print 'eff for central 20%: ',hsmear.Integral(40,60,40,60)/hbox.Integral(40,60,40,60)
print 'eff for central 80%: ',hsmear.Integral(29,71,29,71)/hbox.Integral(29,71,29,71)
    
can.cd()
#hbox.Draw()
#hsmear.SetLineColor(2)
#hsmear.SetMarkerColor(2)
hsmear.RebinX(5)
hsmear.RebinY(5)
hbox.RebinX(5)
hbox.RebinY(5)
hsmear.Divide(hbox)
hsmear.Draw('colz')
can.Update()
can.WaitPrimitive()

hsmear.SetDirectory(f)
hbox.SetDirectory(f)

f.Write()
f.Close()
