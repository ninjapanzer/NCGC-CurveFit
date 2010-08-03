package edu.scripps.fl.hibernate;

import java.util.HashMap;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.AnnotationConfiguration;

public class HibernateStaticService {
    private static HashMap<Object, SessionFactory> hibernateSessionFactories = new HashMap<Object, SessionFactory>();

    public static void buildHibernateSessionFactory(AnnotationConfiguration config) {
        setHibernateSessionFactory(null, config.buildSessionFactory());
    }

    public static void buildHibernateSessionFactory(Object key, AnnotationConfiguration config) {
        setHibernateSessionFactory(key, config.buildSessionFactory());
    }

    protected static void setHibernateSessionFactory(Object key, SessionFactory sf) {
        hibernateSessionFactories.put(key, sf);
    }

    public static SessionFactory getHibernateSessionFactory() {
        return getHibernateSessionFactory(null);
    }

    public static SessionFactory getHibernateSessionFactory(Object key) {
        return hibernateSessionFactories.get(key);
    }

    public static Session getHibernateSession() {
        return getHibernateSession(null);
    }

    public static Session getHibernateSession(Object key) {
        Session session = getHibernateSessionFactory(key).getCurrentSession();
        if( ! session.getTransaction().isActive() )
            session.beginTransaction();
        return session;
    }
}