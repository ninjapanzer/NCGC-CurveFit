package edu.scripps.fl.hibernate;

import java.io.Serializable;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.Collection;
import java.util.List;

import org.hibernate.Hibernate;
import org.hibernate.HibernateException;
import org.hibernate.usertype.UserType;

public abstract class ListStringType<E> implements UserType {
	
	protected abstract String getStringFromList(List<E> list);
	   
	protected abstract List<E> getListFromString(String str);
	
	protected abstract List<E> newList(int capacity);

	public int[] sqlTypes() {
		return new int[] { Types.VARCHAR };
	}

	public Class<? extends Object> returnedClass() {
		return List.class;
	}

	public boolean equals(Object arg0, Object arg1) throws HibernateException {
		return collectionsAreEqual((List<E>) arg0, (List<E>) arg1);
	}

	public boolean collectionsAreEqual(Collection<E> a, Collection<E> a2) {
		return a == null ? a2==null : a.equals(a2);
	}

	public int hashCode(Object arg0) throws HibernateException {
		return ((List) arg0).hashCode();
	}

	public Object nullSafeGet(ResultSet arg0, String[] arg1, Object arg2) throws HibernateException, SQLException {
		String list = (String) Hibernate.STRING.nullSafeGet(arg0, arg1[0]);
		if (list == null)
			return newList(16);
		return getListFromString(list);
	}

	public void nullSafeSet(PreparedStatement arg0, Object arg1, int arg2) throws HibernateException, SQLException {
		String list = getStringFromList((List<E>) arg1); // Arrays.toString((int[])arg1);
		Hibernate.STRING.nullSafeSet(arg0, list, arg2);
	}

	public Object deepCopy(Object arg0) throws HibernateException {
		if (arg0 == null)
			return null;
		List list = (List) arg0;
		List newList = newList(list.size());
		newList.addAll(list);
		return newList;
	}

	public boolean isMutable() {
		return true;
	}

	public Serializable disassemble(Object arg0) throws HibernateException {
		return (Serializable) deepCopy(arg0);
	}

	public Object assemble(Serializable arg0, Object arg1) throws HibernateException {
		return (Serializable) deepCopy(arg0);
	}

	public Object replace(Object arg0, Object arg1, Object arg2) throws HibernateException {
		return (Serializable) deepCopy(arg0);
	}	
}