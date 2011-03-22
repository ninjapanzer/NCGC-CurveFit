/*
<h2 align=center>PUBLIC DOMAIN NOTICE<h2/>

<h3 align=center>NIH Chemical Genomics Center<h3/>
<h3 align=center>National Human Genome Research Institute<h3/>

<p>This software/database is a "United States Government Work" under the 
terms of the United States Copyright Act.  It was written as part of 
the author's official duties as United States Government employee and 
thus cannot be copyrighted.  This software/database is freely 
available to the public for use. The NIH Chemical Genomics Center 
(NCGC) and the U.S. Government have not placed any restriction on its 
use or reproduction. 

<p>Although all reasonable efforts have been taken to ensure the accuracy 
and reliability of the software and data, the NCGC and the U.S. 
Government do not and cannot warrant the performance or results that 
may be obtained by using this software or data. The NCGC and the U.S. 
Government disclaim all warranties, express or implied, including 
warranties of performance, merchantability or fitness for any 
particular purpose. 

<p>Please cite the authors in any work or product based on this material. 
*/

/*
 * CurveClass.java
 *
 * Created on May 9, 2006, 3:23 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package gov.nih.ncgc.batch;

/**
 *
 * @author southalln
 */
public class CurveClass {
    
    private boolean useMask_ = true;
    private boolean allowBell_ = true;
    private double SD_ = 10.;
    private double SDfactor_ = 3.0;
    private double robust_ = 80.;
    private double R2cutoff_ = HillConstants.R2;
    private double asymThresh_ = 0.75;
    
    /** Creates a new instance of CurveClass */
    public CurveClass() {
    }
    
    public void setUsemask(boolean useMask) {useMask_ = useMask;}
    public void setAllowBell(boolean allowBell) {allowBell_ = allowBell;}
    public void setSD(double SD) {SD_ = SD;}
    public void setSDfactor(double SDfactor) {SDfactor_ = SDfactor;}
    public void setRobust(double robust) {robust_ = robust;}
    public void setR2Cutoff(double R2cutoff) {R2cutoff_ = R2cutoff;}
    public void setAsym(double threshold) {asymThresh_ = threshold;}
    
    public double enzCurveClass(HillCurve hc) {
        double maxact = -Double.MAX_VALUE;
        double minact = Double.MAX_VALUE;
        int uSDf = 0, uSD=0, dSDf=0, dSD=0; // for sign code below ...
        for (int i=0; i<hc.getAct().length; i++) {
            if (!useMask_ || hc.getActMask(i)) {
                double act = hc.getAct(i);
                if (maxact < act)
                    maxact = act;
                if (minact > act)
                    minact = act;
                if (act > SD_) {
                    uSD++;
                    if (act > SD_*SDfactor_)
                        uSDf++;
                }
                else if (act < -1*SD_) {
                    dSD++;
                    if (act < -1*SD_*SDfactor_)
                        dSDf++;
                }
            }
        }
        
        int sign = -1; // does curve go up or down?
        if (hc.hasCurve()) {
            // watch out for super actives w/ no change in slope
            if (minact < -1*SD_*SDfactor_ && maxact < -1*SD_*SDfactor_) 
                sign = -1;
            else if (minact > SD_*SDfactor_ && maxact > SD_*SDfactor_) 
                sign = 1;
            else if (hc.getZeroAct() < hc.getInfAct())
                sign = 1;
        }
        else { // inactive, superactive, undefined, poor fits
            if (maxact < SD_*SDfactor_ && minact < -1*SD_*SDfactor_)
                sign = -1;
            else if (maxact > SD_*SDfactor_ && minact > -1*SD_*SDfactor_)
                sign = 1;
            else if (maxact < SD_*SDfactor_ && minact > -1*SD_*SDfactor_)
                return 4.0; // no significant points
            else {
                if (uSDf > dSDf) sign = 1;
                else if (dSDf > uSDf) sign = -1;
                else if (uSD > dSD) sign = 1;
                else sign = -1; // greater probability of freak up than down
            }
        }
        
        double sigval = maxact; // largest significant value and # of pts > SD
        int ssigpts = uSDf;
        int sigpts = uSD;
        if (sign == -1) {
            sigval = minact;
            ssigpts = dSDf;
            sigpts = dSD;
        }
        
        if (hc.hasCurve()) {
//        # Don't use slope of curve -- cell-based assays often have steep hill slopes
//        # Use 2 unmasked pts within 80% of max efficacy
//        #HConcSlope = 0
//        #if ((2.*ac50/HConc) < 1 and hill > -10) or ((ac50/HConc/2.) > 1 and hill < 10):
//        #    HConcAct1 = 100./(1.+pow(ac50/(HConc*2.),hill))
//        #    HConcAct2 = 100./(1.+pow(ac50/(HConc/2.),hill))
//        #    # normalize by (infact-zeroact)
//        #    HConcSlope = (HConcAct1 - HConcAct2)/(math.log(4))
//        # Don't search for lower asymptote ... unnecessary
//        # only super actives have been found so far
//        ## has lower asymptote? first points within SD*SDfactor of zero
//        #pts = 0
//        #for i in range(len(act)):
//        #    if conc[i] < ac50 and abs(act[i]) < SD*SDfactor:
//        #        pts = pts + 1
//        #if pts > 1:
//        #    asymptotes = asymptotes + 1

            // bell shaped curve fits to wrong sig pt
            // not using masks (uSDf, dSDf) makes this too messy
            if (allowBell_ && dSDf > 0 && uSDf > 0) {
                boolean hasBell = false;
                for (int i=1; i<hc.getAct().length; i++)
                    if (hc.getAct(i) == sigval)
                        continue;
                    else if (sign==1 && hc.getAct(i) < -1*SD_*SDfactor_)
                        hasBell = true;
                    else if (sign==-1 && hc.getAct(i) > SD_*SDfactor_)
                        hasBell = true;
                if (hasBell)
                    return 5.0; // requires hand-annotation
            }
            
            int asymptote = 1;
            // has upper asymptote? More than one point within 80% of max efficacy
            int pts=0;
            for (int i=0; i<hc.getAct().length; i++) {
                if (!useMask_ || hc.getActMask(i)) {
                    if (hc.getActConc()[i] > hc.getAC50() && 
                            Math.abs(hc.getAct(i)) > asymThresh_*Math.abs(hc.getInfAct()))
                        pts++;
                }
            }
            if (pts > 1)
                asymptote++;
                    
            if (ssigpts == 1) { // single point outliers
                pts = 0;
                boolean pastSig = false;
                for (int i=0; i<hc.getAct().length; i++) {
                    if (!useMask_ || hc.getActMask(i)) {
                        if (hc.getAct(i) == sigval)
                            pastSig = true;
                        else if (pastSig && (Math.abs(hc.getAct(i)) < SD_))
                            pts++;
                    }
                }
                if (pts > 1)
                    return 4.0;
                // Class 3 single point actives
                pts = 0;
                boolean atsig = false;
                for (int i=0; i<hc.getAct().length; i++) {
                    if (!useMask_ || hc.getActMask(-i-1)) {
                        pts++;
                        if (hc.getAct(-i-1) == sigval)
                            atsig = true;
                        else if (atsig) {
                            atsig = false;
                            if (hc.getAct(-i-1)*sign < SD_ && pts < 4 && asymptote < 2)
                                return 3.0*sign;
                        }
                    }
                }
            }
            // Class 4 inactive
            if (sigpts < 2 || ssigpts < 1)
                return 4.0;
            // Class 1, two asymptotes
            else if (asymptote == 2) {
                if (hc.getR2() > R2cutoff_) {
                    if (Math.abs(sigval) > robust_)
                        return 1.1*sign;
                    else if (ssigpts < 2)
                        return 1.4*sign; // might consider requiring hand-curation here
                    else return 1.2*sign;
                }
                else if (Math.abs(sigval) > robust_)
                    return 1.3*sign;
                return 1.4*sign;
            }
            // Class 2
            else {
                if (hc.getR2() > R2cutoff_) {
                    if (Math.abs(sigval) > robust_)
                        return 2.1*sign;
                    else return 2.2*sign;
                }
                else if (Math.abs(sigval) > robust_)
                    return 2.3*sign;
                return 2.4*sign;
            }
        }
        // Infer curve class from points without curve fit!
        // More work to do here !!!
        else {
            if ((sigpts < 2 || ssigpts < 1) || (ssigpts == 1 && 
                    (sigval != hc.getAct(-1) || hc.getAct(-2)*sign < SD_))) {
                // Class 3, single point curves
                if (ssigpts == 1 && sigval == hc.getAct(-1))
                    return 3.0*sign;
                // Class 4
                return 4.0;
            }
            else if (allowBell_ && dSDf > 0 && uSDf > 0)
                return 5.0;
            else if (hc.getAct(-1) > SD_*SDfactor_ || hc.getAct(-2) > SD_*SDfactor_) 
        	// test if either of last two points is significant 
        	return 5.0;
            else
        	// not able to make call, usually reversed curves (active picoM)
                return 4.0; // might require hand-annotation here
//##         # for superactive and undefined, assume R2 is poor -> class .3 + .4
//##         asymptote = 0
//##         # has lower asymptote? first points within SD*SDfactor of zero
//##         if abs(act[0]) < SD*SDfactor or abs(act[1]) < SD*SDfactor:
//##             asymptote = asymptote + 1
//##         # has upper asymptote? More than one point within SD of sigval
//##         if abs(act[-1]) > Robust or abs(act[-2]) > Robust:
//##             asymptote = asymptote + 1
//##         elif act.index(sigval) != -1 and abs(act[-1] - sigval) < SD*SDfactor and abs(act[-1]) > SD:
//##             asymptote = asymptote + 1
//##         elif act.index(sigval) != -2 and abs(act[-2] - sigval) < SD*SDfactor and abs(act[-2]) > SD:
//##             asymptote = asymptote + 1
//##         elif act.index(sigval) != -3 and abs(act[-3] - sigval) < SD*SDfactor and abs(act[-3]) > SD:
//##             asymptote = asymptote + 1
//##         if asymptote == 2:
//##             if abs(sigval) > Robust:
//##                 return 1.3*sign
//##             return 1.4*sign
//##         elif asymptote == 1:
//##             if abs(sigval) > Robust:
//##                 return 2.3*sign
//##             return 2.4*sign
        }
    }
}
