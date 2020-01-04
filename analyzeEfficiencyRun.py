
import sys, os, glob, shutil, json, math, re, random
import ROOT
import analyzerDigitizer as an
import config 

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)',text) ]


def drawAux(c):
    
    textLeft = ROOT.TLatex()
    textLeft.SetTextFont(42)
    textLeft.SetTextSize(0.04)
    textLeft.SetNDC()
    textLeft.DrawLatex(c.GetLeftMargin(), 0.96, "#bf{CMS} 904,#scale[0.75]{ #it{Work in progress}}")
        
    textRight = ROOT.TLatex()
    textRight.SetNDC()
    textRight.SetTextFont(42)
    textRight.SetTextSize(0.04)
    textRight.SetTextAlign(31)
    textRight.DrawLatex(1.0-c.GetRightMargin(), 0.96, "S%d" % (scanid))
    
def setGraphStyle(g, xlabel, ylabel):

    g.SetLineWidth(2)
    g.SetLineColor(ROOT.kBlue)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(ROOT.kBlue)
    g.GetXaxis().SetTitle(xlabel)
    g.GetYaxis().SetTitle(ylabel)
    g.GetYaxis().SetTitleOffset(1.8)
    g.GetYaxis().SetLabelOffset(2.0*g.GetYaxis().GetLabelOffset())
    g.GetXaxis().SetTitleOffset(1.2)
    g.GetXaxis().SetLabelOffset(2.0*g.GetXaxis().GetLabelOffset())
        
    return g



if __name__ == "__main__":

    ## tag: all the plots and results will be saved in a directory with its tagname
    tag = "Digitizer_efficiency"
    
    ## config: specify the configuration containing the mapping and strip dimensions (see config.py)
    cfg = config.cfg_GRAPHITE
    
    ## dir: ROOT directory of all raw data files 
    dir = "."
    
    ## xMin, xMax: typical voltage range for the analysis 
    xMin, xMax = 6000, 8000
    

    
    ##############################################################################################
    outputdir = "%s/%s/" % (dir, tag)
    #if os.path.exists(outputdir): shutil.rmtree(outputdir) # delete output dir, if exists
    if not os.path.exists(outputdir):os.makedirs(outputdir) # make output dir
    
    
    HVeff = [] # storage of HV eff points
    out = {}


    
    # prepare TGraphs
    g_iMon_BOT = ROOT.TGraphErrors()
    g_iMon_TOP = ROOT.TGraphErrors()
    
    g_muon_cls = ROOT.TGraphErrors()
    
    # necessary to evaluate the error on the parameters at WP using linear interpolation
    g_muon_cls_err = ROOT.TGraphErrors()
    
    g_eff_muon = ROOT.TGraphErrors()
    g_eff_muon_err = ROOT.TGraphErrors()
    
    iMonMax = -999
    
    # get the scan ID from the ROOT file
    files = glob.glob("%s/*CAEN.root" % dir)
    if len(files) == 0: sys.exit("No ROOT files in directory") 
    scanid = int(re.findall(r'\d+', files[0])[0])
    
    # get all ROOT files in the dir
    files.sort(key=natural_keys) # sort on file name, i.e. according to HV points
    for i,CAENFile in enumerate(files):
    
        #if i != 7 : continue
        HVPoint = int(os.path.basename(CAENFile).split("_")[1][2:])
        print "Analyze HV point %d " % HVPoint
        
        saveDir = outputdir + "HV%d/" % HVPoint
        if not os.path.exists(saveDir): os.makedirs(saveDir)

        analyzer = an.Analyzer(dir, saveDir, scanid, HVPoint, "efficiency")
        
        analyzer.loadConfig(cfg) # cfg_GIFPP cfg_FEB12
        analyzer.setVerbose(1)
        
        analyzer.setNoiseTimeWindow(20, 60) # define start and end time of the noise window
        analyzer.setMuonTimeWindow(80, 160) # define start and end time of the muon window
        analyzer.setThreshold(6) # set threshold, expressed as integer times noise stdv 

        # calibrate time from bits to ns, depending on digitizer sampling rate in the ini file
        # 5 GHz/MS:     0.2048     
        # 2.5 GHz/MS:   0.4096      
        # 1 GHz/MS:     0.1024
        analyzer.calibrateTime(0.2048) 
        
        
        analyzer.analyze(False)
        analyzer.write()
        
        #analyzer.DQM() # plot all events over entire time range 

        
        
        

        # load CAEN for HV and currents
        CAEN = ROOT.TFile(CAENFile)
        current_TOP = CAEN.Get("Imon_%s" % cfg['topGapName']).GetMean()
        current_BOT = CAEN.Get("Imon_%s" % cfg['botGapName']).GetMean()
        HV = CAEN.Get("HVeff_%s" % cfg['topGapName']).GetMean()
        if HV < 20: HV = CAEN.Get("HVeff_%s" % cfg['botGapName']).GetMean()
        HVeff.append(HV)
        
        if HV > xMax: xMax = HV
        if HV < xMin: xMin = HV
        
        CAEN.Close()

        
        # load analyzer results
        with open("%s/output.json" % (saveDir)) as f_in: analyzerResults = json.load(f_in)
        analyzerResults = analyzerResults['output_parameters']
        

        if current_BOT > iMonMax: iMonMax = current_BOT
        if current_TOP > iMonMax: iMonMax = current_TOP
    
        g_iMon_BOT.SetPoint(i, HV, current_BOT)
        g_iMon_TOP.SetPoint(i, HV, current_TOP)

  
        
        g_muon_cls.SetPoint(i, HV, analyzerResults['muonCLS'])
        g_muon_cls.SetPointError(i, 0, analyzerResults['muonCLS_err'])
        #g_muon_cls_err.SetPoint(i, 0, analyzerResults['muonCLS_err'])
        
       
        
        g_eff_muon.SetPoint(i, HV, 100.*analyzerResults['efficiencyAbs'])
        g_eff_muon.SetPointError(i, 0, 100.*analyzerResults['efficiencyAbs_err'])

             
    
    # do plotting and fitting
    c = ROOT.TCanvas("c1", "c1", 800, 800)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    
    params = ROOT.TLatex()
    params.SetTextFont(42)
    params.SetTextSize(0.03)
    params.SetNDC()
    
    
    
    
    ############################
    # muon efficiency + fit
    ############################ 
    c.Clear()
    g_eff_muon = setGraphStyle(g_eff_muon, "HV_{eff} (V)", "Efficiency (%)")
    g_eff_muon.GetYaxis().SetRangeUser(0,100)
    g_eff_muon.GetXaxis().SetRangeUser(xMin, xMax)
    g_eff_muon.SetLineWidth(2)
    g_eff_muon.SetLineColor(ROOT.kBlue)
    g_eff_muon.SetMarkerStyle(20)
    g_eff_muon.SetMarkerColor(ROOT.kBlue)
    g_eff_muon.Draw("ALP")


    sigmoid = ROOT.TF1("sigmoid","[0]/(1+exp([1]*([2]-x)))", min(HVeff), max(HVeff));
    sigmoid.SetParName(0,"#epsilon_{max}");
    sigmoid.SetParName(1,"#lambda");
    sigmoid.SetParName(2,"HV_{50%}");
    sigmoid.SetParameter(0, 0.98); 
    sigmoid.SetParameter(1, 0.011); 
    sigmoid.SetParameter(2, 6700);   


    g_eff_muon.Fit("sigmoid")

    fitted = g_eff_muon.GetFunction("sigmoid")
    emax = fitted.GetParameter(0)
    lam = fitted.GetParameter(1)
    hv50 = fitted.GetParameter(2)
    
    emax_err = fitted.GetParError(0)
    lam_err = fitted.GetParError(1)
    hv50_err = fitted.GetParError(2)
    


    WP = (math.log(19)/lam + hv50 + 150)
    dLambdaInverse = lam_err / (lam*lam) # error on 1/lambda
    WP_err = math.sqrt((math.log(19)*dLambdaInverse)*(math.log(19)*dLambdaInverse) + hv50_err*hv50_err) # total error on WP
    out["workingPoint"] = WP
    out["workingPoint_err"] = WP_err
    


    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextAlign(13)
    latex.DrawLatex(0.5, 0.5, "#epsilon_{max} = %.1f %%" % (fitted.GetParameter(0)))
    latex.DrawLatex(0.5, 0.45, "#lambda = %.3f" % fitted.GetParameter(1))
    latex.DrawLatex(0.5, 0.4, "HV_{50%%} = %.1f V" % fitted.GetParameter(2))
    latex.DrawLatex(0.5, 0.35, "WP = %.1f V" % WP)
    latex.DrawLatex(0.5, 0.3, "eff(WP) = %.1f %%" % (fitted.Eval(WP)))
    out["eff"] = fitted.Eval(WP)
    out["effMax"] = fitted.GetParameter(0)
    
    drawAux(c)
    c.SaveAs("%s/muonEfficiency.png" % outputdir)
    c.SaveAs("%s/muonEfficiency.pdf" % outputdir)
    
    
  


    
    ############################
    # gap currents
    ############################
    c.Clear()
    g_iMon_BOT = setGraphStyle(g_iMon_BOT, "HV_{eff} (V)", "Current (#muA)")
    g_iMon_TOP = setGraphStyle(g_iMon_TOP, "HV_{eff} (V)", "Current (#muA)")
    
    g_iMon_BOT.Draw("ALP")
    g_iMon_BOT.GetYaxis().SetRangeUser(0, 1.1*iMonMax)
    
    g_iMon_TOP.SetLineColor(ROOT.kRed)
    g_iMon_TOP.SetMarkerColor(ROOT.kRed)
    g_iMon_TOP.Draw("SAME LP")
    
    params.DrawLatex(0.16, 0.9, "#color[4]{Current top gap at WP: %.2f #muA}" % g_iMon_TOP.Eval(WP))
    params.DrawLatex(0.16, 0.85, "#color[2]{Current bottom gap at WP: %.2f #muA}" % g_iMon_BOT.Eval(WP))
    params.DrawLatex(0.16, 0.8, "Total current at WP: %.2f #muA" % (g_iMon_BOT.Eval(WP)+g_iMon_TOP.Eval(WP)))
    out["iMonTop"] = g_iMon_TOP.Eval(WP)
    out["iMonBot"] = g_iMon_BOT.Eval(WP)
    out["iMonTot"] = (g_iMon_BOT.Eval(WP)+g_iMon_TOP.Eval(WP))
    drawAux(c)
    c.SaveAs("%s/gapCurrents.png" % outputdir)
    c.SaveAs("%s/gapCurrents.pdf" % outputdir)
    
    
    
    ############################
    # muon cluster size
    ############################
    c.Clear()
    g_muon_cls = setGraphStyle(g_muon_cls, "HV_{eff} (V)", "Muon cluster size")
    g_muon_cls.Draw("ALP")
    params.DrawLatex(0.16, 0.9, "Muon cluster size at WP: %.1f" % g_muon_cls.Eval(WP))
    out["muonCLS"] = g_muon_cls.Eval(WP)
    out["muonCLS_err"] = g_muon_cls_err.Eval(WP)
    drawAux(c)
    c.SaveAs("%s/muonCLS.png" % outputdir)
    c.SaveAs("%s/muonCLS.pdf" % outputdir)
    

    
    with open("%s/results.json" % outputdir, 'w') as fp: json.dump(out, fp, indent=4)
    
    
        
