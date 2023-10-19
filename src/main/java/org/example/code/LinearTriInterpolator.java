package org.example.code;

public class LinearTriInterpolator {

    private double[] x;
    private double[] y;
    private double[] z;



    public LinearTriInterpolator(double[] xi, double[] yi, double[] array) {
        this.x = x.clone(); // Utilisation de clone pour éviter la modification des tableaux originaux à l'extérieur de la classe
        this.y = y.clone();
        this.z = z.clone();
    }

    public double interpolate(double xInterp, double yInterp) {
        int n = x.length;
        double result = 0;

        for (int i = 0; i < n; i++) {
            double xi = x[i];
            double yi = y[i];
            double zi = z[i];

            if (isInsideTriangle(xInterp, yInterp, xi, yi, x[(i + 1) % n], y[(i + 1) % n], x[(i + 2) % n], y[(i + 2) % n])) {
                // Interpolation linéaire
                double area = calculateTriangleArea(xi, yi, x[(i + 1) % n], y[(i + 1) % n], x[(i + 2) % n], y[(i + 2) % n]);
                double weight = calculateTriangleArea(xInterp, yInterp, x[(i + 1) % n], y[(i + 1) % n], x[(i + 2) % n], y[(i + 2) % n]) / area;
                result = zi + weight * (z[(i + 1) % n] - zi) + weight * (z[(i + 2) % n] - zi);
                break;
            }
        }

        return result;
    }

    private boolean isInsideTriangle(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3) {
        // Vérifie si le point (x, y) est à l'intérieur du triangle défini par les coordonnées (x1, y1), (x2, y2) et (x3, y3)
        // Cette logique peut être implémentée selon l'algorithme du produit vectoriel ou de barycentric
        // Retourne true si le point est à l'intérieur, sinon false
        return false;
    }

    private double calculateTriangleArea(double x1, double y1, double x2, double y2, double x3, double y3) {
        // Calculer l'aire du triangle défini par les coordonnées (x1, y1), (x2, y2) et (x3, y3)
        return Math.abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0);
    }
}
