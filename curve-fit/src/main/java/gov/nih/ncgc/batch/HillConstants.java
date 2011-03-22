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
 * NCGCConstants.java
 *
 * Created on February 22, 2006, 1:12 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package gov.nih.ncgc.batch;

/**
 *
 * @author yuhong
 */
public class HillConstants 
{
        public static int DATA_POINTS = 25;
        
	//charting
        //public static boolean Y0_AS_100 = false;
        
	public static boolean AUTO_RANGE = true;
        public static boolean AUTO_RANGE2 = false; //for hill cell renderer
        public static double MIN_RANGE = -120.0;
        public static double MAX_RANGE = 120.0;

        //fitting
        public static double MIN_Y_RANGE = 20.0;
        
        public static double Y0 = 10.0;
                
	public static double PS_MIN = 20.0;  //to be removed
	public static double PS_MAX = 30.0; //not used
	public static double PI_MIN = -30.0;  //not used
	public static double PI_MAX = -20.0; //to be removed
        
        public static double SUPER_Y = 30.0;  //do not check range
        
        public static int FIT_ITER_NO  = 2;
        
        public static int PARTIAL_FIT_NO = 5;
        public static int PARTIAL_FIT_MASK_NO = 1;
        
	public static double R2 = 0.30;
	public static double MIN_SLOPE = 0.1;
	public static double MAX_SLOPE = 5.0;
        
        public static double Y0_INF_COEF = 1.2; //limit of y0 and yinf 

        public static boolean P4_FIT = true;
        public static String FIT_TYPE = "P4";
        
        public static double EC50_DELTA = 0.05;
        
        public static boolean REFITTING = false;
        
        //fitting parameters
        public static String FIXED_Y0;
        public static String FIXED_YINF;
        public static String FIXED_SLOPE;
        
        public static boolean UPDATE_FIT_DATA = false;
        
        
        //masking
	public static double P0 = 20.0; //% activity allowed for first data point
        
	public static double PMAX = -250.0; // max % activity; not used

	public static double THETA = 120.0;

        public static double V_DEPTH = 20.0;
        
        public static double TP0 = 20.0; //mask delta cutoff for first point 
	public static double TPV = 30.0; //mask delta cutoff for v shape
	public static double TP = 70.0; //mask delta cutoff
	public static double TPN = 80.0; //mask delta cutoff for last point
	public static double TR = 4.0; //(y-mean)/std
	
	public static int MASK = 8; //frequency
        
        public static boolean bellMask = false; //data point lower than peak
        
        //confidence interval
        public static boolean ciFlag = false;
        
        public static boolean rawData = false; //if true, normalize to 100
       
        
        //curve classification
        public static double CLASSIFICATION_SD = 9;
        
        public static double CLASSIFICATION_SD_FACTOR = 4;
}

	
