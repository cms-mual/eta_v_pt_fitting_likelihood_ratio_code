#made by Ryan Mueller
import ROOT as r
from array import array
from scipy import stats
import time
import math 

test = r.TH1F("test","test", 10,10,10)

c1 = r.TCanvas("Canvas1", "visulaization")

f = r.TFile("histodemo.root")

TDir = f.Get("demo")

plot = TDir.Get("TH2F_pt_vs_eta_distribution")

# define custom fits:
myFit = r.TF1("myFit", "[0]+[1]*x^2+[2]*x^4",0,2)
myFit.SetParName(0, "Yintercept")
myFit.SetParName(1, "XsquaredConstant")
myFit.SetParName(2, "XcubedConstant")
myFit.SetParameter(0,1)
myFit.SetParameter(1,-1)
myFit.SetParameter(2,-1)

myFit_Null = r.TF1("myFit_Null", "[0]+[1]*x^2",0,1)
myFit_Null.SetParName(0, "Yintercept")
myFit_Null.SetParName(1, "XsquaredConstant")
myFit_Null.SetParameter(0,1)
myFit_Null.SetParameter(1,-1)

myFit_high_order= r.TF1("myFit_high_order", "[0]+[1]*x^2+[2]*x^4+[3]*x^6+[4]*x^8+[5]*x^10+[6]*x^12+[7]*x^14",0,7)
myFit_high_order.SetParName(0, "Yintercept")
myFit_high_order.SetParName(1, "XsquaredConstant")
myFit_high_order.SetParName(2, "XcubedConstant")
myFit_high_order.SetParName(3, "X6Constant")
myFit_high_order.SetParName(4, "X8Constant")
myFit_high_order.SetParName(5, "X10Constant")
myFit_high_order.SetParName(6, "X12Constant")
myFit_high_order.SetParName(7, "X14Constant")
myFit_high_order.SetParameter(0,1)
myFit_high_order.SetParameter(1,0)
myFit_high_order.SetParameter(2,0)
myFit_high_order.SetParameter(3,0)
myFit_high_order.SetParameter(4,0)
myFit_high_order.SetParameter(5,0)
myFit_high_order.SetParameter(6,0)
myFit_high_order.SetParameter(7,0)

#define a function that will be used a lot later:
def binFromPt(l_pT, l_numPtBins, l_maxPt):
    return int(float(l_numPtBins)/float(l_maxPt)*float(l_pT))


# initialize variables
ptSlices = []
ptSlicesN = 100
maxPt = 250
etaCount = plot.GetYaxis().GetNbins()
ptCount = plot.GetXaxis().GetNbins()
likelyHoods = []

#returns some fit params on graphs:
r.gStyle.SetOptFit(1100)

# analyze 0to 50 and 50 to 250 pins:
ptSlice_0to50 = r.TH1F("ptSlices_0to50", "ptSlices_0to50", etaCount, -3, 3)
ptSlice_50to250 = r.TH1F("ptSlices_50to250", "ptSlices_50to250", etaCount, -3, 3)

# get the right projections from 2d histogram
ptSlice_0to50 = plot.ProjectionY("ptSlice_0to50", binFromPt(0, ptCount, maxPt), binFromPt(50, ptCount, maxPt))
ptSlice_50to250 = plot.ProjectionY("ptSlice_50to250", binFromPt(50, ptCount, maxPt), binFromPt(250, ptCount, maxPt))

#get fits: L is for log likelihood, S is to return fitting parameters, and Q is to reduce the printouts:
fit_myFit_Null_0to50 = ptSlice_0to50.Fit("myFit_Null", "L S Q", "", -2.2,2.2)
fit_myFit_0to50 = ptSlice_0to50.Fit("myFit", "L S Q","", -2.2,2.2)
fit_myFit_Null_50to250 = ptSlice_50to250.Fit("myFit_Null", "L S Q", "", -2.2,2.2)
fit_myFit_50to250 = ptSlice_50to250.Fit("myFit", "L S Q","", -2.2,2.2)

# compute D values for log likelihood ratio
DNulltoAlt_0to50 = fit_myFit_Null_0to50.MinFcnValue() - fit_myFit_0to50.MinFcnValue()
DNulltoAlt_50to250 = fit_myFit_Null_50to250.MinFcnValue() - fit_myFit_50to250.MinFcnValue()

#print out resuts:

print "1-p value for 0 to 50: ", '\t', stats.chi2.cdf(DNulltoAlt_0to50, 1)
print "1-p value for 50 to 250: ", '\t', stats.chi2.cdf(DNulltoAlt_50to250, 1)

#draw histos
ptSlice_0to50.Draw()
c1.SaveAs("ptSlice_0to50.png")
ptSlice_50to250.Draw()
c1.SaveAs("ptSlice_50to250.png")


# analyze finer ptSlices 
for i in range(ptSlicesN):
    #create list of histograms for slices
    ptSlices.append( r.TH1F("ptSlices_" + str(250.0/ptSlicesN*i), "ptSlices_" + str(250.0/ptSlicesN*i), etaCount, -3, 3))
    
    #fill slices
    ptSlices[i] = plot.ProjectionY("projection Y " + str(i), binFromPt(250.0/ptSlicesN*i, ptCount, maxPt), binFromPt(250.0/ptSlicesN*(i+1), ptCount, maxPt))

    #fit the null and test models
    fit_myFit_high_order = ptSlices[i].Fit("myFit_high_order", "L S Q","", -2.2,2.2)
    fit_myFit = ptSlices[i].Fit("myFit", "L S Q","", -2.2,2.2)
    fit_myFit_Null = ptSlices[i].Fit("myFit_Null", "L S Q", "", -2.2,2.2)

    #get likelyhoods (-2*ln(likelyhood))
    likelyHoods.append( [fit_myFit_Null.MinFcnValue(), fit_myFit.MinFcnValue(), fit_myFit_high_order.MinFcnValue()])

    #drawHistograms
    ptSlices[i].Draw("")
    c1.SaveAs("ptSlices_" + str(250.0/ptSlicesN*i) + ".png")

#print -2ln(likelyhoods)
for i in range(len(likelyHoods)):
    print "Null -2*ln(likelyhood): for  ",'\t', 250.0/ptSlicesN*i ,'\t', " to ",'\t', 250.0/ptSlicesN*(i+1),'\t', ": ",'\t', likelyHoods[i][0],'\t',
    print "terms: 2"
    print "alt 1 -2*ln(likelyhood): for  ",'\t', 250.0/ptSlicesN*i ,'\t', " to ",'\t', 250.0/ptSlicesN*(i+1),'\t', ": ",'\t', likelyHoods[i][1],'\t',
    print "terms: 3"
    print "alt 2 -2*ln(likelyhood): for  ",'\t', 250.0/ptSlicesN*i ,'\t', " to ",'\t', 250.0/ptSlicesN*(i+1),'\t', ": ",'\t', likelyHoods[i][2],'\t',
    print "terms: 8, warning, 8 term fit does not work at high pt as is! "
    

for i in range(len(likelyHoods)):
    #compute p values and report them
    DNulltoAlt = likelyHoods[i][0]-likelyHoods[i][1] 
    DNulltohighOrder = likelyHoods[i][0]-likelyHoods[i][2] 
    #cdf not pdf! cdf is the cumilative distribution, pdf is the probabilty, last factor is the DoF
    print "1-p value of better fit of 3 order vs 2 order: ",'\t',250.0/ptSlicesN*i ,'\t', " to ",'\t', 250.0/ptSlicesN*(i+1),'\t', ": ", '\t', stats.chi2.cdf(DNulltoAlt, 1), '\t' 
