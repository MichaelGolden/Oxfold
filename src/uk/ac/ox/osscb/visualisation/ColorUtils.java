package uk.ac.ox.osscb.visualisation;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.awt.Color;

/**
 *
 * @author Michael Golden
 */
public class ColorUtils {

    public static Color selectBestForegroundColor(Color background, Color c1, Color c2)
    {
        float[] hsbBackground = new float[3];
        float[] hsbC1 = new float[3];
        float[] hsbC2 = new float[3];
        Color.RGBtoHSB(background.getRed(), background.getGreen(), background.getBlue(), hsbBackground);
        Color.RGBtoHSB(c1.getRed(), c1.getGreen(), c1.getBlue(), hsbC1);
        Color.RGBtoHSB(c2.getRed(), c2.getGreen(), c2.getBlue(), hsbC2);

        float c1Diff = Math.abs(hsbBackground[2]-hsbC1[2]);
        float c2Diff = Math.abs(hsbBackground[2]-hsbC2[2]);
        if(c1Diff > c2Diff)
        {
            return c1;
        }
        else
        {
            return c2;
        }
    }

    public static String colorToString(Color c) {
        return "rgba(" + c.getRed() + "," + c.getGreen() + "," + c.getBlue() + "," + c.getAlpha() + ")";
    }

    public static Color getColorFromString(String s) {
        String[] rgba = s.replaceAll("([^0-9,])+", "").split(",");
        return new Color(Integer.parseInt(rgba[0]), Integer.parseInt(rgba[1]), Integer.parseInt(rgba[2]), Integer.parseInt(rgba[3]));
    }
}
