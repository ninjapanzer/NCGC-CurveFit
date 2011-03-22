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
 * BatchCorrector.java
 *
 * Created on May 9, 2006, 5:32 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package gov.nih.ncgc.batch;




import java.sql.PreparedStatement;
import java.util.Comparator;

/**
 *
 *
 *
 * @author wangyuh
 *
 * plate_data: type_index 0-raw data 1-normalized 2-corrected
 *
 * doFormula
 *
 * doNormalization
 *
 * doCorrection
 *
 * doSampleData
 *
 * doHill
 *
 *28463
 */

public class BatchHill
{  
    public static int colNo = 48;   
    public static int rowNo = 32;
    
    public static boolean DEBUG = false;
    
    public static int DATA_POINTS = HillConstants.DATA_POINTS;
    
    private PreparedStatement sampleDataUS;
     
    private boolean p4Fit = true;
   
    
    private int fitIterNo = HillConstants.PARTIAL_FIT_MASK_NO;
    
    private String  fixPara;
      
    private Comparator concComparator;
   
    public static double calcRank(double curveClass, double maxResp, double range)
    {              
        if (curveClass < 0.0)
            return (5+curveClass)*1000000-maxResp;
           
        if (curveClass > 0.0 && curveClass < 4.0)
            return (5-curveClass)*10000+maxResp;
        
        return range;
    }
           
    //ic  or absolute
    public static double iccalc(double y0, double yinf, double ec50, double slope, double pct)
    {                
        if (false) //HillConstants.relativeIC50)
        {
            if (yinf-y0 < 0 || (yinf-y0)/pct < 1.0)
                return 0.0;
            
            return ec50-(1/slope)*Math.log10((yinf-y0)/pct-1);
        }
        else
        {
            if (yinf-y0 <0 || pct-y0 < 0 || (yinf-y0)/(pct-y0) < 1.0)
                return 0.0;
            
            return ec50-(1/slope)*Math.log10((yinf-y0)/(pct-y0)-1);
        }
    }
        
    
     public String composeMaskFlag(boolean[] flags)
     {
            if (flags == null)
                return null;
            
            StringBuffer buffer = new StringBuffer();
            for (int i=0; i<flags.length; i++)
            {
                if (flags[i])
                    buffer.append("0");
                else
                    buffer.append("1");
            }
            
            return buffer.toString().trim();
     }
     
     
     
    
}
 
