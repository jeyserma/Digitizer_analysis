
import sys, os, glob, shutil, json, math, re, random, array
import ROOT

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


class Analyzer():

    fIn = None  # pointer to ROOT file
    t = None    # raw data tree
    
    tag = ""
    savePath = ""
    dqmPath = "" 
    basePath = ""
    
    scanid = -1
    HVPoint = -1
    savedir = ""
    
    verbose = 0
    
    cfg = None
    
    stripArea = -1
    
    # threshold in integer number of noise sigma
    threshold = -1
    
    
    timeVector = array.array('d')
    muonTimeVector = array.array('d') # x-values for the muon time window
    
    # noise and muon time window definitions (in ns)
    noiseTimeWindowBegin = -1
    noiseTimeWindowEnd = -1
    
    muonTimeWindowBegin = -1
    muonTimeWindowEnd = -1
    
    noiseTimeWindowBeginIndex = -1
    noiseTimeWindowEndIndex = -1
    
    muonTimeWindowBeginIndex = -1
    muonTimeWindowEndIndex = -1
    

    
    # results from clusterization
    muonCLS = -1
    muonCLS_err = -1

    # results from efficiency
    efficiencyAbs = -1
    efficiencyAbs_err = -1
   
    
    
    # drawing options (see function __drawAux())
    c1 = None # default square canvas
    c2 = None # default rectangular canvas
    textCMS = "#bf{CMS} 904,#scale[0.75]{ #it{Preliminary}}" # CMS GIF++ TLatex (top left on canvas)
    textAux = None # auxiliary info (top right on canvas)


    def __init__(self, dir, savePath, scanid, HVPoint, scanType):

        self.scanid = scanid
        self.HVPoint = HVPoint
        self.basePath = dir
        self.savePath = savePath
        
        if not os.path.exists(self.savePath): os.makedirs(self.savePath)
   
    
        # get the raw data
        self.fIn = ROOT.TFile("%s/HV%d_DIGITIZER/HV%d_DIGITIZER.root" % (dir, HVPoint, HVPoint))
        self.t = self.fIn.Get("data")


        # default square canvas
        self.c1 = ROOT.TCanvas("c1", "c1", 800, 800)
        self.c1.SetLeftMargin(0.12)
        self.c1.SetRightMargin(0.05)
        self.c1.SetTopMargin(0.05)
        self.c1.SetBottomMargin(0.1)
        
        # default rectangular canvas
        self.c2 = ROOT.TCanvas("c2", "c2", 900, 1200)
        self.c2.SetLeftMargin(0.12)
        self.c2.SetRightMargin(0.13)
        self.c2.SetTopMargin(0.05)
        self.c2.SetBottomMargin(0.1)
        
    def setVerbose(self, verbose):
    
        self.verbose = verbose
        
        
    def loadConfig(self, cfg):
    
        self.stripArea = cfg["stripArea"]   
        self.nStrips = len(cfg["DIG_channels"])
        
        self.DIG_strips = cfg["DIG_strips"]
        self.DIG_strips_mask = cfg["DIG_strips_mask"]
        self.DIG_channels = cfg["DIG_channels"]
 
    def setNoiseTimeWindow(self, start, end):
    
        self.noiseTimeWindowBegin = start
        self.noiseTimeWindowEnd = end
        
    def setMuonTimeWindow(self, start, end):
    
        self.muonTimeWindowBegin = start
        self.muonTimeWindowEnd = end
               
    def setThreshold(self, thrs):
    
        self.threshold = thrs

 
    def calibrateTime(self, sampling):
    
        time_ = self.fIn.Get("time")

        for t in time_: 
        
            t_ = t*sampling
            self.timeVector.append(t_)
            
            if self.noiseTimeWindowBeginIndex == -1 and t_ > self.noiseTimeWindowBegin: self.noiseTimeWindowBeginIndex = t
            if self.noiseTimeWindowEndIndex == -1 and t_ > self.noiseTimeWindowEnd: self.noiseTimeWindowEndIndex = t
            
            if self.muonTimeWindowBeginIndex == -1 and t_ > self.muonTimeWindowBegin: self.muonTimeWindowBeginIndex = t
            if self.muonTimeWindowEndIndex == -1 and t_ > self.muonTimeWindowEnd: self.muonTimeWindowEndIndex = t
            
            if t_ > self.muonTimeWindowBegin and t_ < self.muonTimeWindowEnd: self.muonTimeVector.append(t_)
       

       
        
    # convert DAC unit to VPP/2
    def DACtoV(self, pulse):

        #ret = ROOT.vector('double')()
        ret = array.array('d')
        for i in pulse: ret.append(1000.*(-.5 + i/4096.)) # conversion from DAC to -VPP/2 --> VPP/2
        return ret
        
       

       
    # performs a basic analysis on a single pulse: 
    #  - calculate the noise stdev and offset in noise window
    #  - substract noise offset from pulse
    #  - returns TGraph of pulse with correct timing axis (ns)
    def basePulseAnalysis(self, pulse, muonTimeWindow = True):
    
        # convert raw pulse from DAC to mV
        pulse = self.DACtoV(pulse)

        
        # calculate noise standard deviation and mean (in noise window)
        pulse_noise = array.array('d') # new pulse only in noise window
        for i,p in enumerate(pulse): 
        
            t_ = self.timeVector[i]
            if t_ < self.noiseTimeWindowBegin: continue
            if t_ > self.noiseTimeWindowEnd: break
            pulse_noise.append(p)

        noise_stdv = ROOT.TMath.RMS(len(pulse_noise), pulse_noise)
        noise_mean = ROOT.TMath.Mean(len(pulse_noise), pulse_noise) 

        
        # create new pulse, mean substracted
        pulse_muon = array.array('d')
        time_muon = array.array('d')
        it = 0 # tgraph iterator
        for i,p in enumerate(pulse): 
        
            t_ = self.timeVector[i]
            if muonTimeWindow:
                if t_ > self.muonTimeWindowBegin and t_ < self.muonTimeWindowEnd:                           
                    pulse_muon.append(p-noise_mean)
                    time_muon.append(t_)
                    it += 1
            else:
                pulse_muon.append(p-noise_mean)
                time_muon.append(t_)
                it += 1

        return pulse_muon, time_muon, noise_stdv
        
    def analyze(self, printPulses = False):
 
        print "Run analysis"
        
        nHitsAbs = 0
        nTrig = 0
        
        dir = self.savePath + "eventDisplay"
        #if os.path.exists(dir): shutil.rmtree(dir)
        if not os.path.exists(dir): os.makedirs(dir)
    
        pulses = [] # holding pointers to TTree
        graphs = [] # holding TGraphs per pulse per event, needs to be emptied after each event!
        g_ampl = [] # holding histogram for each channel for amplitude
        g_qint = [] # holding histogram for each channel for charge
        
        g_ampl_tot = ROOT.TH1D("ampl_tot", "Amplitude distribution", 50, 0, 50)
        g_qint_tot = ROOT.TH1D("qint_tot", "Charge distribution", 150, 0, 1500)
        
        for i, ch in enumerate(self.DIG_channels):
            
            print "Load digitizer channel %d" % ch
            pulse = ROOT.vector('double')()
            pulses.append(pulse)
            self.t.SetBranchAddress("pulse_ch%d" % ch, pulses[i]) # , "pulse[1024]/F"
            
            h_ampl = ROOT.TH1D("ampl_ch%d" % ch, "Amplitude channel %d" % ch, 50, 0, 50)
            h_qint = ROOT.TH1D("qint_ch%d" % ch, "Charge channel %d" % ch, 150, 0, 1500)
            g_ampl.append(h_ampl)
            g_qint.append(h_qint)
            
            
        # construct all histograms
        h_clustersize = ROOT.TH1D("clustersize", "Cluster size", 1000, 0, 1000)
            
        
        # loop over all events
        for evNum in range(0, self.t.GetEntries()+1):
        
            self.t.GetEntry(evNum)
            nTrig += 1
            
            miny = 1e10
            maxy = -1e10
            
            thrs = [] # storage of all thresholds for each channel
            ampl = [] # storage of amplitudes (max amplitude) for each channel
            qint = [] # storage of charges (integral) for each channel
            

            
            # fill all TGraphs
            nHits_ = 0 # count # hits over trheshold
            for i, ch in enumerate(self.DIG_channels):
            
                pulse_muon, time_muon, noise_stdv = self.basePulseAnalysis(pulses[i])
                graph = ROOT.TGraph(len(time_muon), time_muon, pulse_muon)
                graph.SetName("g_%d" % ch)
                graphs.append(graph)
                
                thrs_ = -1.0*noise_stdv*self.threshold # compute threshold
                thrs.append(thrs_)
                
                
                if min(pulse_muon) < miny: miny = min(pulse_muon)
                if max(pulse_muon) > maxy: maxy = max(pulse_muon)
                
                isTriggered = (min(pulse_muon) < thrs_) # pulse is negative!
                if isTriggered: nHits_ +=1
                
                
                # amplitude and charge analysis (only if strip is triggered)
                if isTriggered:
                
                    g_ampl[i].Fill(-1.0*min(pulse_muon)) # convert to positive amplitude
                    g_qint[i].Fill(graph.Integral()) # integral of muon pulse (in muon window)
                    
                    g_ampl_tot.Fill(-1.0*min(pulse_muon))
                    g_qint_tot.Fill(graph.Integral())


            if miny < 0: miny_ = 1.1*miny
            else: miny_ = 0.9*miny
            
            if maxy > 0: maxy_ = 1.1*maxy
            else: maxy_ = 0.9*maxy

            if printPulses: self.plotEvent(graphs, dir, evNum, min(self.muonTimeVector), max(self.muonTimeVector), miny_, maxy_, thrs)
            
            if nHits_ > 0: 
            
                nHitsAbs += 1
                h_clustersize.Fill(nHits_) # simple clusterisation: just # strips triggered (no spatial/time constraint)

            
            # clean
            #for g in graphs: g.Delete()
            del graphs[:]
            
            
            
        self.efficiencyAbs = 1.0*nHitsAbs / (1.0*nTrig)
        self.efficiencyAbs_err = math.sqrt(self.efficiencyAbs*(1.0-self.efficiencyAbs)/nTrig)
        
        self.muonCLS = h_clustersize.GetMean()
        self.plotMuonCLS(h_clustersize)
        self.plotSingle(g_ampl_tot, 0, 50, "Amplitude (mV)", "amplitude")
        self.plotSingle(g_qint_tot, 0, 1500, "Charge (a.u.)", "qint")
        self.plotChannels(g_ampl, 0, 50, "Amplitude (mV)", "amplitude")
        self.plotChannels(g_qint, 0, 1500, "Charge (a.u.)", "qint")
        
        
  
    
    def plotSingle(self, h, min_, max_, title, output):
    
      
        self.c1.cd()
        self.c1.Clear()

        if h.Integral() > 1: h.Scale(1.0/h.Integral())
            
        h.Draw("HIST")
        h.SetLineColor(ROOT.kBlue)
        h.GetYaxis().SetRangeUser(0, 1.3*h.GetMaximum())
        h.GetXaxis().SetRangeUser(min_, max_)
        h.SetLineWidth(2)  
            
        h.Draw("HIST SAME")
        h.SetLineWidth(2)  
        h.SetLineColor(ROOT.kBlue)

        h.GetXaxis().SetTitle(title)
        h.GetXaxis().SetTitleOffset(1.2)
        h.GetXaxis().SetLabelOffset(0.005)

        h.GetYaxis().SetTitle("Events (normalized)")   
        h.GetYaxis().SetTitleOffset(1.8)
        h.GetYaxis().SetLabelOffset(0.005)
            

        self.__drawAux(self.c1)
        self.c1.RedrawAxis()
        self.c1.Modify()
        if self.verbose > 0:
            self.c1.SaveAs("%s%s.png" % (self.savePath, output))
            self.c1.SaveAs("%s%s.pdf" % (self.savePath, output))
    
    def plotChannels(self, hists, min_, max_, title, output):
    
        for i, ch in enumerate(self.DIG_channels):
            
            h = hists[i]
      
            self.c1.cd()
            self.c1.Clear()

            leg = ROOT.TLegend(.15, 0.75, .4, .93)
            leg.SetBorderSize(0)
            leg.SetTextSize(0.03)
            leg.SetFillStyle(0)
            
            if h.Integral() > 1: h.Scale(1.0/h.Integral())
            
            h.Draw("HIST")
            h.SetLineColor(ROOT.kBlue)
            h.GetYaxis().SetRangeUser(0, 1.3*h.GetMaximum())
            h.GetXaxis().SetRangeUser(min_, max_)
            h.SetLineWidth(2)  
            
            h.Draw("HIST SAME")
            h.SetLineWidth(2)  
            h.SetLineColor(ROOT.kBlue)

            h.GetXaxis().SetTitle(title)
            h.GetXaxis().SetTitleOffset(1.2)
            h.GetXaxis().SetLabelOffset(0.005)

            h.GetYaxis().SetTitle("Events (normalized)")   
            h.GetYaxis().SetTitleOffset(1.8)
            h.GetYaxis().SetLabelOffset(0.005)
            
            params = ROOT.TLatex()
            params.SetTextFont(42)
            params.SetTextSize(0.03)
            params.SetNDC()
            params.DrawLatex(0.16, 0.9, "Channel %d" % self.DIG_strips[i])
            

            self.__drawAux(self.c1)
            self.c1.RedrawAxis()
            self.c1.Modify()
            if self.verbose > 0:
                self.c1.SaveAs("%s%s_ch%d.png" % (self.savePath, output, self.DIG_strips[i])) 
                self.c1.SaveAs("%s%s_ch%d.pdf" % (self.savePath, output, self.DIG_strips[i])) 
                
        
    
    def plotMuonCLS(self, h_clustersize):
    
        maxCLS = 8
  
        self.c1.cd()
        self.c1.Clear()

        leg = ROOT.TLegend(.15, 0.75, .4, .93)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.03)
        leg.SetFillStyle(0)
        
        if h_clustersize.Integral() > 1: h_clustersize.Scale(1.0/h_clustersize.Integral())
        
        h_clustersize.Draw("HIST")
        h_clustersize.SetLineColor(ROOT.kBlue)
        h_clustersize.GetYaxis().SetRangeUser(0, 1.3*h_clustersize.GetMaximum())
        h_clustersize.GetXaxis().SetRangeUser(0, maxCLS)
        h_clustersize.SetLineWidth(2)  
        
        h_clustersize.Draw("HIST SAME")
        h_clustersize.SetLineWidth(2)  
        h_clustersize.SetLineColor(ROOT.kBlue)

        h_clustersize.GetXaxis().SetTitle("Cluster size")
        h_clustersize.GetXaxis().SetTitleOffset(1.2)
        h_clustersize.GetXaxis().SetLabelOffset(0.005)

        h_clustersize.GetYaxis().SetTitle("Events (normalized)")   
        h_clustersize.GetYaxis().SetTitleOffset(1.8)
        h_clustersize.GetYaxis().SetLabelOffset(0.005)
        
        params = ROOT.TLatex()
        params.SetTextFont(42)
        params.SetTextSize(0.03)
        params.SetNDC()
        params.DrawLatex(0.16, 0.9, "Mean muon cluster size (CLS): %.2f" % (self.muonCLS))
        

        self.__drawAux(self.c1)
        self.c1.RedrawAxis()
        self.c1.Modify()
        if self.verbose > 0:
            self.c1.SaveAs("%sCLS_muon.png" % (self.savePath)) 
            self.c1.SaveAs("%sCLS_muon.pdf" % (self.savePath))   
            
    
    
    def plotEvent(self, graphs, outdir, evNum, minx, maxx, miny, maxy, thrs = []):

        # plot events on a divided canvas
        c = ROOT.TCanvas("evPlot", "c", 1200, 900)
        c.Divide(1, 8, 1e-5, 1e-5)
            
        thrs_lines = [] # hold all line objects

        # loop over all strips
        for i, ch in enumerate(self.DIG_channels):
                
            c.cd(i+1)
            p = c.GetPad(i+1)
            p.SetGrid()
                
            p.SetTopMargin(0.25)
            p.SetBottomMargin(0.05)
            p.SetLeftMargin(0.07)
            p.SetRightMargin(0.05)

            g = graphs[i]
            

            g.SetMarkerStyle(20)
            g.SetMarkerSize(.4)
            g.SetLineWidth(1)
            g.SetLineColor(ROOT.kRed)
            g.SetMarkerColor(ROOT.kRed)

            g.GetYaxis().SetTitleFont(43)
            g.GetYaxis().SetTitleSize(12)
            g.GetYaxis().SetLabelFont(43)
            g.GetYaxis().SetLabelSize(20)
            g.GetYaxis().SetNdivisions(3)
            

            if minx != -1 and maxx != -1: g.GetXaxis().SetRangeUser(minx, maxx)
            if miny != -1 and maxy != -1: g.GetYaxis().SetRangeUser(miny, maxy)
                
            #g.GetXaxis().SetNdivisions(8)
            g.GetXaxis().SetLabelSize(0.15)
            g.GetXaxis().SetLabelOffset(-.9)
            if i%2  == 0: g.GetXaxis().SetLabelSize(0)

            


            g.SetLineWidth(2)
            g.Draw("AL") #AL AXIS X+
            graphs.append(g)


            right = ROOT.TLatex()
            right.SetNDC()
            right.SetTextFont(43)
            right.SetTextSize(20)
            right.SetTextAlign(13)
            right.DrawLatex(.97, .45,"#%d" % self.DIG_strips[i])
            
            if len(thrs) != 0:
            
                line = ROOT.TLine(minx, thrs[i], maxx, thrs[i]);
                line.SetLineColor(ROOT.kBlue)
                line.SetLineWidth(2)
                line.SetLineStyle(2)
                line.Draw()
                thrs_lines.append(line) # add to collector

            p.Update()
            p.Modify()
            c.Update()
                
                

        ### General text on canvas
        c.cd(0)

        

        # toptext
        right = ROOT.TLatex()
        right.SetNDC()
        right.SetTextFont(43)
        right.SetTextSize(20)
        right.SetTextAlign(23)
        right.DrawLatex(.5, .995, "Scan ID: %06d, %s, Event number: %06d  [x-axis: ns, y-axis: mV]" % (self.scanid, self.HVPoint, evNum))

        c.Modify()
            
        

        #c.SaveAs("%s/Scan%06d_HV%s_%d.pdf" % (outdir, self.scanid, self.HVPoint, evNum))
        c.SaveAs("%s/Scan%06d_HV%s_%d.png" % (outdir, self.scanid, self.HVPoint, evNum))               
    
    
        
    def DQM(self):
    
        print "Run DQM"
    
        pulses = [] # holding pointers to TTree
        graphs = [] # holding TGraphs per pulse per event, needs to be emptied after each event!
        
        dir = self.savePath + "DQM"
        #if os.path.exists(dir): shutil.rmtree(dir)
        if not os.path.exists(dir): os.makedirs(dir)
        
        
        for i, ch in enumerate(self.DIG_channels):
            
            print "Load digitizer channel %d" % ch
            pulse = ROOT.vector('double')()
            pulses.append(pulse)
            self.t.SetBranchAddress("pulse_ch%d" % ch, pulses[i]) # , "pulse[1024]/F"
            
        
        # loop over all events
        for evNum in range(0, self.t.GetEntries()+1):
        
            self.t.GetEntry(evNum)
            
            miny = 1e10
            maxy = -1e10
            
            # fill all TGraphs
            for i, ch in enumerate(self.DIG_channels):
                
                pulse_muon, time_muon, noise_stdv = self.basePulseAnalysis(pulses[i], False)
                graph = ROOT.TGraph(len(time_muon), time_muon, pulse_muon)
                graph.SetName("g_%d" % ch)
                graphs.append(graph)
            
                
                if min(pulse_muon) < miny: miny = min(pulse_muon)
                if max(pulse_muon) > maxy: maxy = max(pulse_muon)
            
                
            
            
            if miny < 0: miny_ = 1.1*miny
            else: miny_ = 0.9*miny
            
            if maxy > 0: maxy_ = 1.1*maxy
            else: maxy_ = 0.9*maxy
            
            self.plotEvent(graphs, dir, evNum, min(self.timeVector), max(self.timeVector), miny_, maxy_)

            del graphs[:]
        
        
        
    
        
    def __drawAux(self, c, aux = ""):
    
        textLeft = ROOT.TLatex()
        textLeft.SetTextFont(42)
        textLeft.SetTextSize(0.04)
        textLeft.SetNDC()
        textLeft.DrawLatex(c.GetLeftMargin(), 0.96, self.textCMS)
        
        textRight = ROOT.TLatex()
        textRight.SetNDC()
        textRight.SetTextFont(42)
        textRight.SetTextSize(0.04)
        textRight.SetTextAlign(31)
        if aux == "": textRight.DrawLatex(1.0-c.GetRightMargin(), 0.96, "S%d/HV%d" % (self.scanid, self.HVPoint))
        else: textRight.DrawLatex(1.0-c.GetRightMargin(), 0.96, "S%d/HV%d/%s" % (self.scanid, self.HVPoint, aux))

        
    def validateEvent(self):

        
        ## Quality flag validation (see Alexis: https://github.com/afagot/GIF_OfflineAnalysis/blob/master/src/utils.cc)
        qFlag = self.t.Quality_flag
        #print self.t.Quality_flag
        tmpflag = qFlag
        
        IsCorrupted = False
        nDigits = 0
        while tmpflag / int(math.pow(10, nDigits)) != 0: nDigits += 1;
        
        while not IsCorrupted and nDigits != 0:
        
            tdcflag = tmpflag / int(math.pow(10, nDigits-1))

            if tdcflag == 2: 
                IsCorrupted = True

            tmpflag = tmpflag % int(math.pow(10,nDigits-1))
            nDigits -= 1
        
        return not IsCorrupted
        

        ## PMT validation
        '''
        if len(self.TDC_channels_PMT) == 0: return True ## NO VALIDATION
        for ch in self.TDC_channels_PMT:
            if not ch in self.t.TDC_channel: return False
        
        return True ## default
        '''
        
    # Input: raw TDC channel/time vectors,
    # Output: converted TDC channels to strip numbers, within the optinally given time window
    def __groupAndOrder(self, TDC_CH, TDC_TS, windowStart = -1e9, windowEnd = 1e9):
    
        STRIP = []
        TS = []
        for i,ch in enumerate(TDC_CH):
            if not ch in self.TDC_channels: continue # only consider channels from chamber
            if TDC_TS[i] < windowStart: continue # min time window
            if TDC_TS[i] > windowEnd: continue # max time window
            #if TDC_TS[i] < self.timeWindowReject: continue # reject TDC first events
            #stripNo = cfg.TDC_strips[cfg.TDC_channels.index(ch)]
            stripNo = self.TDC_channels.index(ch)
            STRIP.append(stripNo)
            TS.append(TDC_TS[i])
        
        return STRIP, TS
 
            


   
    def write(self):
    
        print "Write output JSON file"
        
        out = {}
        
        param_input = {

            "threshold"                 : self.threshold,
            
            "noiseTimeWindowBegin"      : self.noiseTimeWindowBegin,
            "noiseTimeWindowEnd"        : self.noiseTimeWindowEnd,
            
            "muonTimeWindowBegin"       : self.muonTimeWindowBegin,
            "muonTimeWindowEnd"         : self.muonTimeWindowEnd,
            
        }
        
        param_output = {

    
            "muonCLS"                   : self.muonCLS,
            "muonCLS_err"               : self.muonCLS_err,
    
            "efficiencyAbs"             : self.efficiencyAbs,
            "efficiencyAbs_err"         : self.efficiencyAbs_err,


        }        
   
        data = {
        
            "input_parameters"          :  param_input, 
            "output_parameters"         :  param_output, 
        }
    
        with open("%soutput.json" % self.savePath, 'w') as fp: json.dump(data, fp, indent=4)
    


