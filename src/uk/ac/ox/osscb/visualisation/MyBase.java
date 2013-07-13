package uk.ac.ox.osscb.visualisation;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


/**
 *
 * @author Michael Golden <michaelgolden0@gmail.com>
 */
public class MyBase extends ModeleBase {

    @Override
    public int getIndex() {
        return this._realIndex;
    }

    @Override
    public String getContent() {
        return "";
    }

    @Override
    public void setContent(String s) {
        
    }
    
    public String toString()
    {
        return "" +this._coords.toPoint2D();
    }
    
}
