package edu.scripps.test;

import java.net.URL;

import edu.scripps.fl.curves.Curve;
import edu.scripps.fl.curves.CurveFit;
import edu.scripps.fl.hibernate.HibernateUtil;

public class HibTest {

	public static void main(String[] args) throws Exception {
		URL url = HibTest.class.getResource("/hibernate.cfg.xml");
        HibernateUtil.build(url);
		Curve curve = CurveTest.fullCurveDecrease();
		CurveFit.fit(curve);
		HibernateUtil.getSession().save(curve);
		HibernateUtil.getSession().getTransaction().commit();
	}
}
