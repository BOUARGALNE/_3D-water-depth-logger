package org.example.code;

import org.example.code.LinearTriInterpolator;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.apache.commons.math3.analysis.function.Exp;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TriangulationExample {

    public static void main(String[] args) {
        ArrayList<Double> x = new ArrayList<>();
        ArrayList<Double> y = new ArrayList<>();
        ArrayList<Double> z = new ArrayList<>();
        ArrayList<String[]> tidalData = new ArrayList<>();


        try (BufferedReader br = new BufferedReader(new FileReader("C:\\Users\\Hamid\\Documents\\Semestre5\\stream-kafka\\calcul-parallele-distribue\\TP2\\_3D-water-depth-logger\\src\\main\\java\\org\\example\\tidal_data.csv"))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] data = line.split(",");
                tidalData.add(data);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        try (BufferedReader br = new BufferedReader(new FileReader("C:\\Users\\Hamid\\Documents\\Semestre5\\stream-kafka\\calcul-parallele-distribue\\TP2\\_3D-water-depth-logger\\src\\main\\java\\org\\example\\grohn26082021_raw.csv"))) {
            String line;
       while ((line = br.readLine()) != null) {
                String[] row = line.split(",");
                for (String[] tidal : tidalData) {

                    StringBuilder sb = new StringBuilder(tidal[0].substring(1,tidal[0].length()-1));
                    StringBuilder sb2 = new StringBuilder(row[1].substring(0, row[1].length() - 3));
                    if (sb.toString().equals(sb2.toString()) && 0.3 < Double.parseDouble(row[2]) && Double.parseDouble(row[2]) < 2) {
                        x.add(Double.parseDouble(row[4]));
                        y.add(Double.parseDouble(row[3]));
                        z.add(4 - (Double.parseDouble(tidal[1].substring(1,tidal[1].length()-1)) - Double.parseDouble(row[2])));
                        break;
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        double latmin = findMinValue(x);
        double latmax = findMaxValue(x);
        double lngmin = findMinValue(y);
        double lngmax = findMaxValue(y);

        int npts = 150;
        int ngridx = 150;
        int ngridy = 150;

        double[] xi = linspace(latmin, latmax, ngridx);
        double[] yi = linspace(lngmin, lngmax, ngridy);

        double[][][] Xi = createMeshGrid(xi, yi);
/*        // Convertir les coordonnées en objets Coordinate pour la triangulation
        Coordinate[] coordinates = new Coordinate[x.size()];
        for (int i = 0; i < x.size(); i++) {
            coordinates[i] = new Coordinate(x.get(i), y.get(i));
        }

        // Créer un objet MultiPoint pour la triangulation
        GeometryFactory geometryFactory = new GeometryFactory();
        MultiPoint multiPoint = geometryFactory.createMultiPointFromCoords(coordinates);

        // Effectuer la triangulation de Delaunay
        DelaunayTriangulationBuilder triangulator = new DelaunayTriangulationBuilder();
        triangulator.setSites(multiPoint);
        Polygon[] triangles = getTriangles(triangulator.getTriangles(geometryFactory));
      for (Polygon triangle : triangles) {
            // Faites quelque chose avec le triangle, par exemple, accédez aux sommets du triangle
            Coordinate[] triangleCoordinates = triangle.getCoordinates();
            for (Coordinate vertex : triangleCoordinates) {
                double xVertex = vertex.getX();
                double yVertex = vertex.getY();

            }
        }*/

 /*       LinearTriInterpolator interpolator = new LinearTriInterpolator(xi, yi, z.stream().mapToDouble(Double::doubleValue).toArray());
        double[][] zi = new double[ngridx][ngridy];

        for (int i = 0; i < ngridx; i++) {
            for (int j = 0; j < ngridy; j++) {
                zi[i][j] = interpolator.interpolate(xi[i], yi[j]);
            }
        }*/
        double[] x1=listToArray(x);
        double[] y1=listToArray(x);
        double[] z1=listToArray(x);
        double[][] zi = interpolateValues(xi, yi, x1, y1, z1);

    }
    public static double[] listToArray(ArrayList<Double> list) {
        double[] array = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }
/*    private static Polygon[] getTriangles(org.locationtech.jts.geom.Geometry trianglesGeometry) {
        List<Polygon> triangles = new ArrayList<>();
        for (int i = 0; i < trianglesGeometry.getNumGeometries(); i++) {
            if (trianglesGeometry.getGeometryN(i) instanceof Polygon) {
                triangles.add((Polygon) trianglesGeometry.getGeometryN(i));
            }
        }
        return triangles.toArray(new Polygon[0]);
    }*/
    private static double[][] interpolateValues(double[] xi, double[] yi, double[] x, double[] y, double[] z) {
        double[][] zi = new double[xi.length][yi.length];

        LinearInterpolator interpolator = new LinearInterpolator();
        for (int i = 0; i < xi.length; i++) {
            for (int j = 0; j < yi.length; j++) {
                // Interpolation 1D pour obtenir les valeurs zi
                PolynomialSplineFunction splineFunctionX = interpolator.interpolate(x, z);
                PolynomialSplineFunction splineFunctionY = interpolator.interpolate(y, z);

                // Interpolation 2D en utilisant les fonctions 1D
                double interpolatedValueX = splineFunctionX.value(xi[i]);
                double interpolatedValueY = splineFunctionY.value(yi[j]);

                // Calculer la valeur interpolée en (xi, yi)
                zi[i][j] = interpolatedValueX + interpolatedValueY;
            }
        }

        return zi;
    }
    public static double[][][] createMeshGrid(double[] xi, double[] yi) {
        double[][][] grids = new double[2][yi.length][xi.length];

        for (int i = 0; i < yi.length; i++) {
            for (int j = 0; j < xi.length; j++) {
                grids[0][i][j] = xi[j];
                grids[1][i][j] = yi[i];
            }
        }

        return grids;
    }
    private static double findMinValue(ArrayList<Double> list) {
        double minValue = list.get(0);
        for (Double value : list) {
            if (value < minValue) {
                minValue = value;
            }
        }
        return minValue;
    }

    private static double findMaxValue(ArrayList<Double> list) {
        double maxValue = list.get(0);
        for (Double value : list) {
            if (value > maxValue) {
                maxValue = value;
            }
        }
        return maxValue;
    }

    private static double[] linspace(double start, double end, int n) {
        double[] array = new double[n];
        double step = (end - start) / (n - 1);
        for (int i = 0; i < n; i++) {
            array[i] = start + i * step;
        }
        return array;
    }


        }
