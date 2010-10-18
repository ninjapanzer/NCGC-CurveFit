/**
 * 
 */
package edu.scripps.fl.hibernate;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections.list.GrowthList;
import org.apache.commons.lang.StringUtils;

public class DoubleListStringType extends ListStringType<Double> {

	public List<Double> getListFromString(String list) {
		String strs[] = list.split("\r?\n");
		List<Double> ids = newList(strs.length);
		for (int ii = 0; ii < strs.length; ii++) {
			if (null != strs[ii] && !strs[ii].equals("")) {
				Double dbl = Double.parseDouble(strs[ii]);
				ids.set(ii, dbl);
			}
		}
		return ids;
	}

	public String getStringFromList(List<Double> ids) {
		if (null == ids)
			return null;
		return StringUtils.join(ids, "\n");
	}
	
	public List<Double> newList(int initialCapacity) {
		return (List<Double>) GrowthList.decorate(new ArrayList<Double>(initialCapacity));
	}	
}