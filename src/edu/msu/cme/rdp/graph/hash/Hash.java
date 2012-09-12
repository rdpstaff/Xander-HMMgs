/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.graph.hash;

/**
 *
 * @author wangqion
 */
public interface Hash {

    public long eatRight(long hashvalue, int c);
    public long eatLeft(long hashvalue, int c);

    //public long getInitialHashvalue(String s);
    public long updateLeft(long hashvalue , int outint, int inint );
    public long updateRight(long hashvalue , int outint, int inint );

    public long replaceRight(long hashvalue, int oldC, int newC) ;
    public long replaceLeft( long hashvalue, int oldC, int newC);
}
