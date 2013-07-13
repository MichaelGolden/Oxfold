package uk.ac.ox.osscb.visualisation;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import fr.orsay.lri.varna.models.CubicBezierCurve;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.templates.*;
import fr.orsay.lri.varna.models.treealign.Tree;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Vector;

/**
 *
 * @author Michael Golden <michaelgolden0@gmail.com>
 */
public class RadiateView {

    /**
     * Selects the "Feynman diagram" drawing algorithm that places the bases on
     * a circle and draws the base-pairings as chords of the circle graph.
     */
    public static final int DRAW_MODE_CIRCULAR = 1;
    /**
     * Selects the "tree drawing" algorithm. Draws each loop on a circle whose
     * radius depends on the number of bases involved in the loop. As some
     * helices can be overlapping in the result, basic interaction is provided
     * so that the user can "disentangle" the drawing by spinning the helices
     * around the axis defined by their multiloop (bulge or internal loop)
     * origin. This is roughly the initial placement strategy of RNAViz.
     *
     * @see <a href="http://rnaviz.sourceforge.net/">RNAViz</a>
     */
    public static final int DRAW_MODE_RADIATE = 2;
    /**
     * Selects the NAView algorithm.
     */
    public static final int DRAW_MODE_NAVIEW = 3;
    /**
     * Selects the linear algorithm.
     */
    public static final int DRAW_MODE_LINEAR = 4;
    public static final int DRAW_MODE_VARNA_VIEW = 5;
    /**
     * Selects the RNAView algorithm.
     */
    public static final int DRAW_MODE_MOTIFVIEW = 6;
    public static final int DRAW_MODE_TEMPLATE = 7;
    public static final int DEFAULT_DRAW_MODE = DRAW_MODE_RADIATE;
    private Double _spaceBetweenBases = 1.0;
    /**
     * The draw algorithm mode
     */
    private int _drawMode = DRAW_MODE_RADIATE;
    public int BASE_RADIUS = 10;
    public static final double LOOP_DISTANCE = 40.0; // distance between base pairs in an helix
   // public static final double BASE_PAIR_DISTANCE = 65.0; // distance between the two bases of a pair
     public static final double BASE_PAIR_DISTANCE = 45.0; // distance between the two bases of a pair
    public static final double MULTILOOP_DISTANCE = 40.0;
    public static final double VIRTUAL_LOOP_RADIUS = 40.0;
    public double CHEM_PROB_DIST = 14;
    public double CHEM_PROB_BASE_LENGTH = 30;
    public double CHEM_PROB_ARROW_HEIGHT = 10;
    public double CHEM_PROB_ARROW_WIDTH = 5;
    public double CHEM_PROB_TRIANGLE_WIDTH = 2.5;
    public double CHEM_PROB_PIN_SEMIDIAG = 6;
    public double CHEM_PROB_DOT_RADIUS = 6.;
    public GeneralPath _debugShape = null;
    private boolean _drawn = false;
    public double _bpHeightIncrement = VARNAConfig.DEFAULT_BP_INCREMENT;
    /**
     * the base list
     */
    protected ArrayList<ModeleBase> _listeBases;

    public static ArrayList<Point2D.Double> radiateview_xy_coordinates(int[] pairedSites, boolean flatExteriorLoop) throws Exception {

        ArrayList<ModeleBase> pairs = new ArrayList<ModeleBase>();
        for (int i = 0; i < pairedSites.length; i++) {
            MyBase base = new MyBase();
            base.setBaseNumber(i);
            base.setElementStructure(pairedSites[i] - 1, new ModeleBP(""));
            pairs.add(base);
        }


        RadiateView rna = new RadiateView();
        rna._listeBases = pairs;
        rna.drawRNARadiate(-1, flatExteriorLoop);
        ArrayList<Point2D.Double> coordinates = new ArrayList<Point2D.Double>();
        for (ModeleBase base : rna._listeBases) {
            Point2D.Double coord = base._coords.toPoint2D();
            coord.x = coord.x/2;
            coord.y = coord.y/2;
            coordinates.add(coord);
        }
        return coordinates;
    }

    public void drawRNARadiate(double dirAngle, boolean flatExteriorLoop) {
        _drawn = true;
        _drawMode = DRAW_MODE_RADIATE;
        Point2D.Double[] coords = new Point2D.Double[_listeBases.size()];
        Point2D.Double[] centers = new Point2D.Double[_listeBases.size()];
        for (int i = 0; i < _listeBases.size(); i++) {
            coords[i] = new Point2D.Double(0, 0);
            centers[i] = new Point2D.Double(0, 0);
        }
        if (flatExteriorLoop) {
            dirAngle += 1.0 - Math.PI / 2.0;
            int i = 0;
            double x = 0.0;
            double y = 0.0;
            double vx = -Math.sin(dirAngle);
            double vy = Math.cos(dirAngle);
            while (i < _listeBases.size()) {
                coords[i].x = x;
                coords[i].y = y;
                centers[i].x = x + BASE_PAIR_DISTANCE * vy;
                centers[i].y = y - BASE_PAIR_DISTANCE * vx;
                int j = _listeBases.get(i).getElementStructure();
                if (j > i) {
                    drawLoop(i, j, x + (BASE_PAIR_DISTANCE * vx / 2.0), y + (BASE_PAIR_DISTANCE * vy / 2.0), dirAngle, coords, centers);
                    centers[i].x = coords[i].x + BASE_PAIR_DISTANCE * vy;
                    centers[i].y = y - BASE_PAIR_DISTANCE * vx;
                    i = j;
                    x += BASE_PAIR_DISTANCE * vx;
                    y += BASE_PAIR_DISTANCE * vy;
                    centers[i].x = coords[i].x + BASE_PAIR_DISTANCE * vy;
                    centers[i].y = y - BASE_PAIR_DISTANCE * vx;
                }
                x += MULTILOOP_DISTANCE * vx;
                y += MULTILOOP_DISTANCE * vy;
                i += 1;
            }
        } else {
            drawLoop(0, _listeBases.size() - 1, 0, 0, dirAngle, coords, centers);
        }
        for (int i = 0; i < _listeBases.size(); i++) {
            _listeBases.get(i).setCoords(
                    new Point2D.Double(coords[i].x * _spaceBetweenBases,
                    coords[i].y * _spaceBetweenBases));
            _listeBases.get(i).setCenter(
                    new Point2D.Double(centers[i].x * _spaceBetweenBases,
                    centers[i].y * _spaceBetweenBases));
        }

        // TODO
        // change les centres des bases de la premiere helice vers la boucle la
        // plus proche
    }

    private double drawHelixLikeTemplateHelix(
            int[] basesInHelixArray, // IN
            RNATemplate.RNATemplateHelix helix, // IN  (optional, ie. may be null)
            Point2D.Double[] coords, // OUT (optional, ie. may be null)
            Point2D.Double[] centers, // OUT (optional, ie. may be null)
            double scaleHelixOrigin, // IN
            Map<RNATemplate.RNATemplateHelix, Point2D.Double> translateVectors // IN (optional, ie. may be null)
            ) {
        int n = basesInHelixArray.length / 2;
        if (n == 0) {
            return 0;
        }
        // Default values when not template helix is provided:
        Point2D.Double o = new Point2D.Double(0, 0);
        Point2D.Double i = new Point2D.Double(1, 0);
        Point2D.Double j = new Point2D.Double(0, 1);
        boolean flipped = false;
        if (helix != null) {
            computeTemplateHelixVectors(helix, o, i, j);
            flipped = helix.isFlipped();
        }
        Point2D.Double li = new Point2D.Double(i.x * LOOP_DISTANCE, i.y * LOOP_DISTANCE);
        // We want o to be the point where the first base (5' end) is
        o.x = (o.x - j.x * BASE_PAIR_DISTANCE / 2) * scaleHelixOrigin;
        o.y = (o.y - j.y * BASE_PAIR_DISTANCE / 2) * scaleHelixOrigin;
        if (translateVectors != null && translateVectors.containsKey(helix)) {
            Point2D.Double v = translateVectors.get(helix);
            o.x = o.x + v.x;
            o.y = o.y + v.y;
        }

        // We need this array so that we can store positions even if coords == null
        Point2D.Double[] helixBasesPositions = new Point2D.Double[basesInHelixArray.length];
        for (int k = 0; k < helixBasesPositions.length; k++) {
            helixBasesPositions[k] = new Point2D.Double();
        }
        Point2D.Double accDelta = new Point2D.Double(0, 0);
        for (int k = 0; k < n; k++) {
            int kp = 2 * n - k - 1;
            Point2D.Double p1 = helixBasesPositions[k]; // we assign the point *reference*
            Point2D.Double p2 = helixBasesPositions[kp];
            // Do we have a bulge between previous base pair and this one?
            boolean bulge = k >= 1 && (basesInHelixArray[k] != basesInHelixArray[k - 1] + 1
                    || basesInHelixArray[kp + 1] != basesInHelixArray[kp] + 1);
            if (k >= 1) {
                if (basesInHelixArray[k] < basesInHelixArray[k - 1]
                        || basesInHelixArray[kp + 1] < basesInHelixArray[kp]) {
                    throw new Error("Internal bug: basesInHelixArray must be sorted");
                }
                if (bulge) {
                    // Estimate a good distance (delta) between the previous base pair and this one
                    double delta1 = estimateBulgeWidth(basesInHelixArray[k - 1], basesInHelixArray[k]);
                    double delta2 = estimateBulgeWidth(basesInHelixArray[kp], basesInHelixArray[kp + 1]);
                    // The biggest bulge defines the width
                    double delta = Math.max(delta1, delta2);

                    if (coords != null) {
                        // Now, where do we put the bases that are part of the bulge?
                        for (int side = 0; side < 2; side++) {
                            Point2D.Double pstart = new Point2D.Double();
                            Point2D.Double pend = new Point2D.Double();
                            Point2D.Double bisectVect = new Point2D.Double();
                            Point2D.Double is = new Point2D.Double();
                            int firstBase, lastBase;
                            double alphasign = flipped ? -1 : 1;
                            if (side == 0) {
                                firstBase = basesInHelixArray[k - 1];
                                lastBase = basesInHelixArray[k];
                                pstart.setLocation(o.x + accDelta.x,
                                        o.y + accDelta.y);
                                pend.setLocation(o.x + accDelta.x + i.x * delta,
                                        o.y + accDelta.y + i.y * delta);
                                bisectVect.setLocation(-j.x, -j.y);
                                is.setLocation(i);
                            } else {
                                firstBase = basesInHelixArray[kp];
                                lastBase = basesInHelixArray[kp + 1];
                                pstart.setLocation(o.x + accDelta.x + i.x * delta + j.x * BASE_PAIR_DISTANCE,
                                        o.y + accDelta.y + i.y * delta + j.y * BASE_PAIR_DISTANCE);
                                pend.setLocation(o.x + accDelta.x + j.x * BASE_PAIR_DISTANCE,
                                        o.y + accDelta.y + j.y * BASE_PAIR_DISTANCE);

                                bisectVect.setLocation(j);
                                is.setLocation(-i.x, -i.y);
                            }
                            double arclen = estimateBulgeArcLength(firstBase, lastBase);
                            double centerOnBisect = RadiateView.ComputeArcCenter.computeArcCenter(delta, arclen);

                            // Should we draw the base on an arc or simply use a line?
                            if (centerOnBisect > -1000) {
                                Point2D.Double center = new Point2D.Double(pstart.x + is.x * delta / 2 + bisectVect.x * centerOnBisect,
                                        pstart.y + is.y * delta / 2 + bisectVect.y * centerOnBisect);
                                int b = firstBase;
                                double len = 0;
                                double r = Math.hypot(pstart.x - center.x, pstart.y - center.y);
                                double alpha0 = angleFromVector(pstart.x - center.x, pstart.y - center.y);
                                while (b < lastBase) {
                                    int l = _listeBases.get(b).getElementStructure();
                                    if (b < l && l < lastBase) {
                                        len += BASE_PAIR_DISTANCE;
                                        b = l;
                                    } else {
                                        len += LOOP_DISTANCE;
                                        b++;
                                    }
                                    if (b < lastBase) {
                                        coords[b].x = center.x + r * Math.cos(alpha0 + alphasign * len / r);
                                        coords[b].y = center.y + r * Math.sin(alpha0 + alphasign * len / r);
                                    }
                                }
                            } else {
                                // Draw on a line
                                double nBP = 0;
                                double nLD = 0;
                                {
                                    int b = firstBase;
                                    while (b < lastBase) {
                                        int l = _listeBases.get(b).getElementStructure();
                                        if (b < l && l < lastBase) {
                                            nBP++;
                                            b = l;
                                        } else {
                                            nLD++;
                                            b++;
                                        }
                                    }
                                }
                                // Distance between paired bases cannot be changed
                                // (imposed by helix width) but distance between other
                                // bases can be adjusted.
                                double LD = Math.max((delta - nBP * BASE_PAIR_DISTANCE) / nLD, 0);
                                //System.out.println("nBP=" + nBP + " nLD=" + nLD);
                                double len = 0;
                                {
                                    int b = firstBase;
                                    while (b < lastBase) {
                                        int l = _listeBases.get(b).getElementStructure();
                                        if (b < l && l < lastBase) {
                                            len += BASE_PAIR_DISTANCE;
                                            b = l;
                                        } else {
                                            len += LD;
                                            b++;
                                        }
                                        if (b < lastBase) {
                                            coords[b].x = pstart.x + is.x * len;
                                            coords[b].y = pstart.y + is.y * len;
                                        }
                                    }
                                }
                                //System.out.println("len=" + len + " delta=" + delta + " d(pstart,pend)=" + Math.hypot(pend.x-pstart.x, pend.y-pstart.y));
                            }

                            // Does the bulge contain an helix?
                            // If so, use drawLoop() to draw it.
                            {
                                int b = firstBase;
                                while (b < lastBase) {
                                    int l = _listeBases.get(b).getElementStructure();
                                    if (b < l && l < lastBase) {
                                        // Helix present in bulge
                                        Point2D.Double b1pos = coords[b];
                                        Point2D.Double b2pos = coords[l];
                                        double beta = angleFromVector(b2pos.x - b1pos.x, b2pos.y - b1pos.y) - Math.PI / 2 + (flipped ? Math.PI : 0);
                                        Point2D.Double loopCenter = new Point2D.Double((b1pos.x + b2pos.x) / 2, (b1pos.y + b2pos.y) / 2);
                                        drawLoop(b,
                                                l,
                                                loopCenter.x,
                                                loopCenter.y,
                                                beta,
                                                coords,
                                                centers);
                                        // If the helix is flipped, we need to compute the symmetric
                                        // of the whole loop.
                                        if (helix.isFlipped()) {
                                            Point2D.Double v = new Point2D.Double(Math.cos(beta), Math.sin(beta));
                                            Point2D.Double[] points1 = new Point2D.Double[l - b + 1];
                                            Point2D.Double[] points2 = new Point2D.Double[l - b + 1];
                                            for (int c = b; c <= l; c++) {
                                                // This is an assignment by reference.
                                                points1[c - b] = coords[c];
                                                points2[c - b] = centers[c];
                                            }
                                            symmetric(loopCenter, v, points1);
                                            symmetric(loopCenter, v, points2);
                                        }
                                        // Continue
                                        b = l;
                                    } else {
                                        b++;
                                    }
                                }
                            }
                        }
                    }

                    accDelta.x += i.x * delta;
                    accDelta.y += i.y * delta;
                    p1.x = o.x + accDelta.x;
                    p1.y = o.y + accDelta.y;
                    p2.x = p1.x + j.x * BASE_PAIR_DISTANCE;
                    p2.y = p1.y + j.y * BASE_PAIR_DISTANCE;

                } else {
                    accDelta.x += li.x;
                    accDelta.y += li.y;
                    p1.x = o.x + accDelta.x;
                    p1.y = o.y + accDelta.y;
                    p2.x = p1.x + j.x * BASE_PAIR_DISTANCE;
                    p2.y = p1.y + j.y * BASE_PAIR_DISTANCE;
                }
            } else {
                // First base pair
                p1.x = o.x;
                p1.y = o.y;
                p2.x = p1.x + j.x * BASE_PAIR_DISTANCE;
                p2.y = p1.y + j.y * BASE_PAIR_DISTANCE;
            }
        }

        Point2D.Double p1 = helixBasesPositions[0];
        Point2D.Double p2 = helixBasesPositions[n - 1];

        if (coords != null) {
            for (int k = 0; k < helixBasesPositions.length; k++) {
                coords[basesInHelixArray[k]] = helixBasesPositions[k];
            }
        }

        return Math.hypot(p2.x - p1.x, p2.y - p1.y);
    }

    private void drawHelixLikeTemplateHelixOLD(
            int[] basesInHelixArray, // IN
            RNATemplate.RNATemplateHelix helix, // IN
            Point2D.Double[] coords, // OUT
            Point2D.Double[] centers, // OUT
            double scaleHelixOrigin, // IN
            Map<RNATemplate.RNATemplateHelix, Point2D.Double> translateVectors // IN
            ) {
        int n = basesInHelixArray.length / 2;
        Point2D.Double o = new Point2D.Double();
        Point2D.Double i = new Point2D.Double();
        Point2D.Double j = new Point2D.Double();
        computeTemplateHelixVectors(helix, o, i, j);
        i.x = i.x * LOOP_DISTANCE;
        i.y = i.y * LOOP_DISTANCE;
        o.x = (o.x - j.x * BASE_PAIR_DISTANCE / 2) * scaleHelixOrigin;
        o.y = (o.y - j.y * BASE_PAIR_DISTANCE / 2) * scaleHelixOrigin;
        if (translateVectors != null && translateVectors.containsKey(helix)) {
            Point2D.Double v = translateVectors.get(helix);
            o.x = o.x + v.x;
            o.y = o.y + v.y;
        }

        for (int k = 0; k < n; k++) {
            int b1 = basesInHelixArray[k];
            int b2 = basesInHelixArray[2 * n - k - 1];
            Point2D.Double p1 = coords[b1]; // we assign the point *reference*
            p1.x = o.x + k * i.x;
            p1.y = o.y + k * i.y;
            Point2D.Double p2 = coords[b2];
            p2.x = p1.x + j.x * BASE_PAIR_DISTANCE;
            p2.y = p1.y + j.y * BASE_PAIR_DISTANCE;
        }

        for (int k = 0; k < 2 * n - 1; k++) {
            if (k == n - 1) {
                continue;
            }
            int b1 = basesInHelixArray[k];
            int b2 = basesInHelixArray[k + 1];
            if (b1 + 1 != b2) {
                // There is a loop between these 2 bases
                Point2D.Double b1pos = coords[b1];
                Point2D.Double b2pos = coords[b2];
                // ||v|| = 1
                Point2D.Double v = new Point2D.Double();
                if (k >= n) {
                    v.x = j.x;
                    v.y = j.y;
                } else {
                    v.x = -j.x;
                    v.y = -j.y;
                }
                Point2D.Double loopCenter = new Point2D.Double();
                loopCenter.x = (b1pos.x + b2pos.x) / 2 + v.x * LOOP_DISTANCE;
                loopCenter.y = (b1pos.y + b2pos.y) / 2 + v.y * LOOP_DISTANCE;
                drawLoop(b1 + 1,
                        b2 - 1,
                        loopCenter.x,
                        loopCenter.y,
                        angleFromVector(v),
                        coords,
                        centers);
                // If the helix is flipped, we need to compute the symmetric
                // of the whole loop.
                if (helix.isFlipped()) {
                    int from = b1 + 1;
                    int to = b2 - 1;
                    Point2D.Double[] points1 = new Point2D.Double[to - from + 1];
                    Point2D.Double[] points2 = new Point2D.Double[to - from + 1];
                    for (int b = from; b <= to; b++) {
                        // This is an assignment by reference.
                        points1[b - from] = coords[b];
                        points2[b - from] = centers[b];
                    }
                    symmetric(loopCenter, v, points1);
                    symmetric(loopCenter, v, points2);
                }
            }
        }
    }

    /**
     * A Bezier curve can be defined by four points, see
     * http://en.wikipedia.org/wiki/Bezier_curve#Cubic_B.C3.A9zier_curves Here
     * we give this four points and an array of bases indexes (which must be
     * indexes in this RNA sequence) which will be moved to be on the Bezier
     * curve. The bases positions are not changed in fact, instead the coords
     * and centers arrays are modified.
     */
    private void drawOnBezierCurve(int[] basesInSequence,
            Point2D.Double P0,
            Point2D.Double P1,
            Point2D.Double P2,
            Point2D.Double P3,
            Point2D.Double[] coords,
            Point2D.Double[] centers) {
        // Draw the bases of the sequence along a Bezier curve
        int n = basesInSequence.length;
        // We choose to approximate the Bezier curve by 10*n straight lines.
        CubicBezierCurve bezier = new CubicBezierCurve(P0, P1, P2, P3, 10 * n);
        double curveLength = bezier.getApproxCurveLength();
        double delta_t = curveLength / (n + 1);
        double[] t = new double[n];
        for (int k = 0; k < n; k++) {
            t[k] = (k + 1) * delta_t;
        }
        Point2D.Double[] sequenceBasesCoords = bezier.uniformParam(t);
        for (int k = 0; k < n; k++) {
            coords[basesInSequence[k]] = sequenceBasesCoords[k];
        }
    }

    /**
     * Like drawOnBezierCurve(), but on a straight line.
     */
    private void drawOnStraightLine(int[] basesInSequence,
            Point2D.Double P0,
            Point2D.Double P3,
            Point2D.Double[] coords,
            Point2D.Double[] centers) {
        // Draw the bases of the sequence along a Bezier curve
        int n = basesInSequence.length;
        Point2D.Double v = new Point2D.Double(P3.x - P0.x, P3.y - P0.y);
        for (int k = 0; k < n; k++) {
            coords[basesInSequence[k]].x = P0.x + (k + 1) / (double) (n + 1) * v.x;
            coords[basesInSequence[k]].y = P0.y + (k + 1) / (double) (n + 1) * v.y;
        }
    }

    /**
     * This functions draws the RNA sequence between (including) firstBase and
     * lastBase along a curve. The sequence may contain helices.
     *
     * Bezier curve: A Bezier curve can be defined by four points, see
     * http://en.wikipedia.org/wiki/Bezier_curve#Cubic_B.C3.A9zier_curves
     *
     * Straight line: If P1 and P2 are null, the bases are drawn on a straight
     * line.
     *
     * OUT: The bases positions are not changed in fact, instead the coords and
     * centers arrays are modified.
     *
     */
    private void drawAlongCurve(
            int firstBase,
            int lastBase,
            Point2D.Double P0,
            Point2D.Double P1,
            Point2D.Double P2,
            Point2D.Double P3,
            Point2D.Double[] coords,
            Point2D.Double[] centers) {

        // First we find the bases which are directly on the Bezier curve
        ArrayList<Integer> alongBezierCurve = new ArrayList<Integer>();
        for (int depth = 0, i = firstBase; i <= lastBase; i++) {
            int k = _listeBases.get(i).getElementStructure();
            if (k < 0 || k > lastBase || k < firstBase) {
                if (depth == 0) {
                    alongBezierCurve.add(i);
                }
            } else {
                if (i < k) {
                    if (depth == 0) {
                        alongBezierCurve.add(i);
                        alongBezierCurve.add(k);
                    }
                    depth++;
                } else {
                    depth--;
                }
            }
        }
        // number of bases along the Bezier curve
        int n = alongBezierCurve.size();
        int[] alongBezierCurveArray = RNATemplateAlign.intArrayFromList(alongBezierCurve);
        if (n > 0) {
            if (P1 != null && P2 != null) {
                drawOnBezierCurve(alongBezierCurveArray, P0, P1, P2, P3, coords, centers);
            } else {
                drawOnStraightLine(alongBezierCurveArray, P0, P3, coords, centers);
            }
        }
        // Now use the radiate algorithm to draw the helixes
        for (int k = 0; k < n - 1; k++) {
            int b1 = alongBezierCurveArray[k];
            int b2 = alongBezierCurveArray[k + 1];
            if (_listeBases.get(b1).getElementStructure() == b2) {
                Point2D.Double b1pos = coords[b1];
                Point2D.Double b2pos = coords[b2];
                double alpha = angleFromVector(b2pos.x - b1pos.x, b2pos.y - b1pos.y);
                drawLoop(b1,
                        b2,
                        (b1pos.x + b2pos.x) / 2,
                        (b1pos.y + b2pos.y) / 2,
                        alpha - Math.PI / 2,
                        coords,
                        centers);
            }
        }
    }

    /**
     * Compute the symmetric of all the points in the points array relative to
     * the line that goes through p and has director vector v. The array is
     * modified in-place.
     */
    private static void symmetric(
            Point2D.Double p,
            Point2D.Double v,
            Point2D.Double[] points) {
        // ||v||^2
        double lv = v.x * v.x + v.y * v.y;
        for (int i = 0; i < points.length; i++) {
            // A is the coordinates of points[i] after moving the origin at p
            Point2D.Double A = new Point2D.Double(points[i].x - p.x, points[i].y - p.y);
            // Symmetric of A
            Point2D.Double B = new Point2D.Double(
                    -(A.x * v.y * v.y - 2 * A.y * v.x * v.y - A.x * v.x * v.x) / lv,
                    (A.y * v.y * v.y + 2 * A.x * v.x * v.y - A.y * v.x * v.x) / lv);
            // Change the origin back
            points[i].x = B.x + p.x;
            points[i].y = B.y + p.y;
        }
    }

    private void computeHelixTranslations(
            Tree<RNANodeValueTemplate> tree, // IN
            Map<RNATemplate.RNATemplateHelix, Point2D.Double> translateVectors, // OUT (must be initialized)
            RNATemplateMapping mapping, // IN
            Point2D.Double parentDeltaVector // IN
            ) {
        RNANodeValueTemplate nvt = tree.getValue();
        Point2D.Double newDeltaVector = parentDeltaVector;
        if (nvt instanceof RNANodeValueTemplateBasePair) {
            RNATemplate.RNATemplateHelix helix = ((RNANodeValueTemplateBasePair) nvt).getHelix();
            if (!translateVectors.containsKey(helix)) {
                translateVectors.put(helix, parentDeltaVector);
                int[] basesInHelixArray;
                if (mapping.getAncestor(helix) != null) {
                    basesInHelixArray = RNATemplateAlign.intArrayFromList(mapping.getAncestor(helix));
                    Arrays.sort(basesInHelixArray);
                } else {
                    basesInHelixArray = new int[0];
                }
                Point2D.Double helixDeltaVector = computeLengthIncreaseDelta(basesInHelixArray, helix);
                newDeltaVector = new Point2D.Double(parentDeltaVector.x + helixDeltaVector.x, parentDeltaVector.y + helixDeltaVector.y);
            }
        }
        for (Tree<RNANodeValueTemplate> subtree : tree.getChildren()) {
            computeHelixTranslations(subtree, translateVectors, mapping, newDeltaVector);
        }
    }

    private Point2D.Double computeLengthIncreaseDelta(
            int[] basesInHelixArray, // IN
            RNATemplate.RNATemplateHelix helix // IN
            ) {
        double templateLength = computeHelixTemplateLength(helix);
        double realLength = computeHelixRealLength(basesInHelixArray);
        Point2D.Double i = new Point2D.Double();
        computeTemplateHelixVectors(helix, null, i, null);
        return new Point2D.Double(i.x * (realLength - templateLength), i.y * (realLength - templateLength));
    }

    private void drawLoop(int i, int j, double x, double y, double dirAngle,
            Point2D.Double[] coords, Point2D.Double[] centers) {
        if (i > j) {
            return;
        }

        // BasePaired
        if (_listeBases.get(i).getElementStructure() == j) {
            double normalAngle = Math.PI / 2.0;
            centers[i] = new Point2D.Double(x, y);
            centers[j] = new Point2D.Double(x, y);
            coords[i].x = (x + BASE_PAIR_DISTANCE
                    * Math.cos(dirAngle - normalAngle) / 2.0);
            coords[i].y = (y + BASE_PAIR_DISTANCE
                    * Math.sin(dirAngle - normalAngle) / 2.0);
            coords[j].x = (x + BASE_PAIR_DISTANCE
                    * Math.cos(dirAngle + normalAngle) / 2.0);
            coords[j].y = (y + BASE_PAIR_DISTANCE
                    * Math.sin(dirAngle + normalAngle) / 2.0);
            drawLoop(i + 1, j - 1, x + LOOP_DISTANCE * Math.cos(dirAngle), y
                    + LOOP_DISTANCE * Math.sin(dirAngle), dirAngle, coords,
                    centers);
        } else {
            int k = i;
            Vector<Integer> basesMultiLoop = new Vector<Integer>();
            Vector<Integer> helices = new Vector<Integer>();
            int l;
            while (k <= j) {
                l = _listeBases.get(k).getElementStructure();
                if (l > k) {
                    basesMultiLoop.add(new Integer(k));
                    basesMultiLoop.add(new Integer(l));
                    helices.add(new Integer(k));
                    k = l + 1;
                } else {
                    basesMultiLoop.add(new Integer(k));
                    k++;
                }
            }
            int mlSize = basesMultiLoop.size() + 2;
            int numHelices = helices.size() + 1;
            double totalLength = MULTILOOP_DISTANCE * (mlSize - numHelices)
                    + BASE_PAIR_DISTANCE * numHelices;
            double multiLoopRadius;
            double angleIncrementML;
            double angleIncrementBP;
            if (mlSize > 3) {
                multiLoopRadius = determineRadius(numHelices, mlSize
                        - numHelices, (totalLength) / (2.0 * Math.PI),
                        BASE_PAIR_DISTANCE,
                        MULTILOOP_DISTANCE);
                angleIncrementML = -2.0
                        * Math.asin(((float) MULTILOOP_DISTANCE)
                        / (2.0 * multiLoopRadius));
                angleIncrementBP = -2.0
                        * Math.asin(((float) BASE_PAIR_DISTANCE)
                        / (2.0 * multiLoopRadius));
            } else {
                multiLoopRadius = 35.0;
                angleIncrementBP = -2.0
                        * Math.asin(((float) BASE_PAIR_DISTANCE)
                        / (2.0 * multiLoopRadius));
                angleIncrementML = (-2.0 * Math.PI - angleIncrementBP) / 2.0;
            }
            // System.out.println("MLr:"+multiLoopRadius+" iBP:"+angleIncrementBP+" iML:"+angleIncrementML);

            double centerDist = Math.sqrt(Math.max(Math.pow(multiLoopRadius, 2)
                    - Math.pow(BASE_PAIR_DISTANCE / 2.0, 2), 0.0))
                    - LOOP_DISTANCE;
            Point2D.Double mlCenter = new Point2D.Double(
                    (x + (centerDist * Math.cos(dirAngle))),
                    (y + (centerDist * Math.sin(dirAngle))));

            // Base directing angle for (multi|hairpin) loop, from the center's
            // perspective
            double baseAngle = dirAngle
                    // U-turn
                    + Math.PI
                    // Account for already drawn supporting base-pair
                    + 0.5 * angleIncrementBP
                    // Base cannot be paired twice, so next base is at
                    // "unpaired base distance"
                    + 1.0 * angleIncrementML;
            double[] angles = new double[_listeBases.size()];
            int n1 = 1;
            int n2 = 1;
            for (k = basesMultiLoop.size() - 1; k >= 0; k--) {
                l = basesMultiLoop.get(k).intValue();
                centers[l] = mlCenter;
                angles[l] = baseAngle;
                coords[l].x = mlCenter.x + multiLoopRadius
                        * Math.cos(baseAngle);
                coords[l].y = mlCenter.y + multiLoopRadius
                        * Math.sin(baseAngle);
                if ((_listeBases.get(l).getElementStructure() < l)
                        && (_listeBases.get(l).getElementStructure() != -1)) {
                    baseAngle += angleIncrementBP;
                    n1++;
                } else {
                    baseAngle += angleIncrementML;
                    n2++;
                }
            }
            // System.out.println("n1:"+n1+" n2:"+n2);
            double newAngle;
            int m, n;
            for (k = 0; k < helices.size(); k++) {
                m = helices.get(k).intValue();
                n = _listeBases.get(m).getElementStructure();
                newAngle = (angles[m] + angles[n]) / 2.0;
                drawLoop(m + 1, n - 1, (LOOP_DISTANCE * Math.cos(newAngle))
                        + (coords[m].x + coords[n].x) / 2.0,
                        (LOOP_DISTANCE * Math.sin(newAngle))
                        + (coords[m].y + coords[n].y) / 2.0, newAngle,
                        coords, centers);
            }
        }
    }

    /**
     * Get helix length in template.
     */
    private double computeHelixTemplateLength(RNATemplate.RNATemplateHelix helix) {
        return Math.hypot(helix.getStartPosition().x - helix.getEndPosition().x,
                helix.getStartPosition().y - helix.getEndPosition().y);
    }

    /**
     * Compute helix actual length (as drawHelixLikeTemplateHelix() would draw
     * it).
     */
    private double computeHelixRealLength(int[] basesInHelixArray) {
        return drawHelixLikeTemplateHelix(basesInHelixArray, null, null, null, 0, null);
    }

    private void computeBezierTangentVectorTarget(
            RNATemplate.RNATemplateElement.EdgeEndPoint endPoint,
            Point2D.Double curveEndPoint,
            Point2D.Double curveVectorOtherPoint)
            throws RNATemplateDrawingAlgorithmException {

        boolean sequenceEndPointIsIn;
        RNATemplate.RNATemplateUnpairedSequence sequence;

        if (endPoint.getElement() instanceof RNATemplate.RNATemplateHelix) {
            sequence = (RNATemplate.RNATemplateUnpairedSequence) endPoint.getOtherElement();
            RNATemplate.EdgeEndPointPosition endPointPositionOnHelix = endPoint.getPosition();
            switch (endPointPositionOnHelix) {
                case IN1:
                case IN2:
                    sequenceEndPointIsIn = false;
                    break;
                default:
                    sequenceEndPointIsIn = true;
            }

            RNATemplate.RNATemplateElement.EdgeEndPoint endPointOnHelix =
                    sequenceEndPointIsIn
                    ? sequence.getIn().getOtherEndPoint()
                    : sequence.getOut().getOtherEndPoint();
            if (endPointOnHelix == null) {
                throw (new RNATemplateDrawingAlgorithmException("Sequence is not connected to an helix."));
            }
        } else {
            // The endpoint is on an unpaired sequence.
            sequence = (RNATemplate.RNATemplateUnpairedSequence) endPoint.getElement();
            if (endPoint == sequence.getIn()) {
                // endpoint is 5'
                sequenceEndPointIsIn = true;
            } else {
                sequenceEndPointIsIn = false;
            }
        }

        double l =
                sequenceEndPointIsIn
                ? sequence.getInTangentVectorLength()
                : sequence.getOutTangentVectorLength();

        // Compute the absolute angle our line makes to the helix
        double theta =
                sequenceEndPointIsIn
                ? sequence.getInTangentVectorAngle()
                : sequence.getOutTangentVectorAngle();

        // Compute v, the tangent vector of the Bezier curve
        Point2D.Double v = new Point2D.Double();
        v.x = l * Math.cos(theta);
        v.y = l * Math.sin(theta);
        curveVectorOtherPoint.x = curveEndPoint.x + v.x;
        curveVectorOtherPoint.y = curveEndPoint.y + v.y;
    }

    /**
     * Compute the angle made by a vector.
     */
    private static double angleFromVector(Point2D.Double v) {
        return angleFromVector(v.x, v.y);
    }

    private static double angleFromVector(double x, double y) {
        double l = Math.hypot(x, y);
        if (y > 0) {
            return Math.acos(x / l);
        } else if (y < 0) {
            return -Math.acos(x / l);
        } else {
            return x > 0 ? 0 : Math.PI;
        }
    }

    private static class ComputeArcCenter {

        /**
         * Given an arc length (l) and segment length (delta) of the arc, find
         * where to put the center, returned as a position of the perpendicular
         * bisector of the segment. The positive side is the one where the arc
         * is drawn. It works using Newton's method.
         */
        public static double computeArcCenter(double delta, double l) {
            double x_n = 0;
            double x_n_plus_1, f_x_n, f_x_n_plus_1;
            int steps = 0;
            while (true) {
                f_x_n = f(x_n, delta);
                x_n_plus_1 = x_n - (f_x_n - l) / fprime(x_n, delta);
                f_x_n_plus_1 = f(x_n_plus_1, delta);
                steps++;
                // We want a precision of 0.1 on arc length
                if (x_n_plus_1 == Double.NEGATIVE_INFINITY || Math.abs(f_x_n_plus_1 - f_x_n) < 0.1) {
                    //System.out.println("computeArcCenter: steps = " + steps + "    result = " + x_n_plus_1);
                    return x_n_plus_1;
                }
                x_n = x_n_plus_1;
                f_x_n = f_x_n_plus_1;
            }
        }

        public static double f(double c, double delta) {
            if (c < 0) {
                return 2 * Math.atan(delta / (-2 * c)) * Math.sqrt(delta * delta / 4 + c * c);
            } else if (c != 0) { // c > 0
                return (2 * Math.PI - 2 * Math.atan(delta / (2 * c))) * Math.sqrt(delta * delta / 4 + c * c);
            } else { // c == 0
                return Math.PI * Math.sqrt(delta * delta / 4 + c * c);
            }
        }

        /**
         * d/dc f(c,delta)
         */
        public static double fprime(double c, double delta) {
            if (c < 0) {
                return delta / (c * c + delta / 4) * Math.sqrt(delta * delta / 4 + c * c) + 2 * Math.atan(delta / (-2 * c)) * c / Math.sqrt(delta * delta / 4 + c * c);
            } else if (c != 0) { // c > 0
                return delta / (c * c + delta / 4) * Math.sqrt(delta * delta / 4 + c * c) + (2 * Math.PI - 2 * Math.atan(delta / (-2 * c))) * c / Math.sqrt(delta * delta / 4 + c * c);
            } else { // c == 0
                return 2;
            }
        }
    }

    /**
     * Estimate bulge arc length.
     */
    private double estimateBulgeArcLength(int firstBase, int lastBase) {
        if (firstBase + 1 == lastBase) {
            return LOOP_DISTANCE; // there is actually no bulge
        }
        double len = 0.0;
        int k = firstBase;
        while (k < lastBase) {
            int l = _listeBases.get(k).getElementStructure();
            if (k < l && l < lastBase) {
                len += BASE_PAIR_DISTANCE;
                k = l;
            } else {
                len += LOOP_DISTANCE;
                k++;
            }
        }
        return len;
    }

    /**
     * Estimate bulge width, the given first and last bases must be those in the
     * helix.
     */
    private double estimateBulgeWidth(int firstBase, int lastBase) {
        double len = estimateBulgeArcLength(firstBase, lastBase);
        return 2 * (len / Math.PI);
    }

    private void computeTemplateHelixVectors(
            RNATemplate.RNATemplateHelix helix, // IN
            Point2D.Double o, // OUT
            Point2D.Double i, // OUT
            Point2D.Double j // OUT
            ) {
        Point2D.Double startpos, endpos;
        if (helix.getIn1Is() == RNATemplate.In1Is.IN1_IS_5PRIME) {
            startpos = helix.getStartPosition();
            endpos = helix.getEndPosition();
        } else {
            endpos = helix.getStartPosition();
            startpos = helix.getEndPosition();
        }
        if (o != null) {
            o.x = startpos.x;
            o.y = startpos.y;
        }
        if (i != null || j != null) {
            // (i_x,i_y) is the vector between two consecutive bases of the same side of an helix
            if (i == null) {
                i = new Point2D.Double();
            }
            i.x = (endpos.x - startpos.x);
            i.y = (endpos.y - startpos.y);
            double i_original_norm = Math.hypot(i.x, i.y);
            // change its norm to 1
            i.x = i.x / i_original_norm;
            i.y = i.y / i_original_norm;
            if (j != null) {
                j.x = -i.y;
                j.y = i.x;
                if (helix.isFlipped()) {
                    j.x = -j.x;
                    j.y = -j.y;
                }
                double j_original_norm = Math.hypot(j.x, j.y);
                // change (j_x,j_y) so that its norm is 1
                j.x = j.x / j_original_norm;
                j.y = j.y / j_original_norm;
            }
        }
    }

    private static double objFun(int n1, int n2, double r, double bpdist, double multidist) {
        return (((double) n1) * 2.0
                * Math.asin(((double) bpdist) / (2.0 * r))
                + ((double) n2) * 2.0
                * Math.asin(((double) multidist) / (2.0 * r)) - (2.0 * Math.PI));
    }

    public double determineRadius(int nbHel, int nbUnpaired, double startRadius) {
        return determineRadius(nbHel, nbUnpaired, startRadius, BASE_PAIR_DISTANCE, MULTILOOP_DISTANCE);
    }

    public static double determineRadius(int nbHel, int nbUnpaired, double startRadius, double bpdist, double multidist) {
        double xmin = bpdist / 2.0;
        double xmax = 3.0 * multidist + 1;
        double x = (xmin + xmax) / 2.0;
        double y = 10000.0;
        double ymin = -1000.0;
        double ymax = 1000.0;
        int numIt = 0;
        double precision = 0.00001;
        while ((Math.abs(y) > precision) && (numIt < 10000)) {
            x = (xmin + xmax) / 2.0;
            y = objFun(nbHel, nbUnpaired, x, bpdist, multidist);
            ymin = objFun(nbHel, nbUnpaired, xmax, bpdist, multidist);
            ymax = objFun(nbHel, nbUnpaired, xmin, bpdist, multidist);
            if (ymin > 0.0) {
                xmax = xmax + (xmax - xmin);
            } else if ((y <= 0.0) && (ymax > 0.0)) {
                xmax = x;
            } else if ((y >= 0.0) && (ymin < 0.0)) {
                xmin = x;
            } else if (ymax < 0.0) {
                xmin = Math.max(xmin - (x - xmin), Math.max(
                        bpdist / 2.0, multidist / 2.0));
                xmax = x;
            }
            numIt++;
        }
        return x;
    }
}
