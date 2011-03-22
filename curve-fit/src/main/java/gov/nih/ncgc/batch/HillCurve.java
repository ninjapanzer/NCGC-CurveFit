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
 * HillCurve.java
 *
 * Created on May 9, 2006, 3:30 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package gov.nih.ncgc.batch;

import java.util.Vector;

/**
 *
 * @author southalln
 */
public class HillCurve {
    
    private String name_;  // Compound name
    private double qAC50_; // qAC50, AC50 + conc must be in same units!
    private double hill_; // hillCoefficient
    private double zeroAct_; // zeroActivity
    private double infAct_; // infiniteActivity
    private double r2_; // R2
    private int partialFit_; // partial fit?
    
    private double[] activity_; // titrated activity values
    private double[] conc_; // titration points
    private boolean[] mask_; // whether titration points are 'valid'
    private double[] objs_; // see HillFit class
    private boolean hasCurve_ = false;
    
    /** Creates a new instance of HillCurve */
    public HillCurve (double[] act, double[] conc, boolean[] mask, String name) {
        name_ = name;
        activity_ = act;
        conc_ = conc;
        mask_ = mask;
    }
    
    public HillCurve (double[] act, double[] conc) {
        activity_ = act;
        conc_ = conc;
    }
    
    public HillCurve (Vector act, Vector conc) {
        activity_ = new double[act.size()];
        conc_ = new double[act.size()];
        for (int i=0; i<act.size(); i++) {
            activity_[i] = ((Number)act.get(i)).doubleValue();
            conc_[i] = ((Number)conc.get(i)).doubleValue();
        }
    }
    
    public void setCurve(double qAC50, 
            double hill, 
            double zeroAct, 
            double infAct,
            double r2,
            int partialFit) {
        qAC50_ = qAC50;
        hill_ = hill;
        zeroAct_ = zeroAct;
        infAct_ = infAct;
        r2_ = r2;
        partialFit_ = partialFit;
        hasCurve_ = true;
    }
    
    public String getName() {return name_;}
    public double getAC50() {return qAC50_;}
    public double getHill() {return hill_;}
    public double getZeroAct() {return zeroAct_;}
    public double getInfAct() {return infAct_;}
    public double getR2() {return r2_;}
    public double[] getAct() {return activity_;}
    public double getAct(int index) {
        if (index < 0) {
            if (-1*index > activity_.length-1)
                return activity_[0];
            return activity_[activity_.length+index];
        }
        else if (index > activity_.length-1)
            return activity_[activity_.length-1];
        else return activity_[index];
    }
    public double[] getActConc() {return conc_;}
    public boolean[] getActMask() {return mask_;}
    public boolean getActMask(int index) {
        if (mask_ == null)
            return true;
        if (index > mask_.length - 1)
            return false;
        if (index < 0) {
            if (-1*index > mask_.length-1)
                return false;
            return mask_[mask_.length+index];
        }
        return mask_[index];
    }
    public double[] getObjs() {return objs_;}
    public void setObjs(double[] objs) {objs_ = objs;}
    public boolean hasCurve() {return hasCurve_;}
}
