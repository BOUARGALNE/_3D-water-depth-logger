package org.example.code;

import org.apache.commons.math3.analysis.interpolation.BicubicInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolator;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Main {
    private static final GeometryFactory geometryFactory = new GeometryFactory();
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

        double[][][] Mi = createMeshGrid(xi, yi);
        double[][]Xi=new double[Mi[0].length][Mi[0][0].length];
        double[][]Yi=new double[Mi[0].length][Mi[1][0].length];
        for (int i = 0; i < Mi[0].length; i++) {
            for (int j = 0; j < Mi[0][0].length; j++) {
                Xi[i][j]=Mi[0][i][j];
                Yi[i][j]=Mi[1][i][j];
            }

        }
      // Convertir les coordonnées en objets Coordinate pour la triangulation
        Coordinate[] coordinates = new Coordinate[x.size()];
        for (int i = 0; i < x.size(); i++) {
            coordinates[i] = new Coordinate(x.get(i), y.get(i));
        }

        // Créer un objet MultiPoint pour la triangulation
        GeometryFactory geometryFactory1 = new GeometryFactory();
        MultiPoint multiPoint = geometryFactory1.createMultiPointFromCoords(coordinates);

        // Effectuer la triangulation de Delaunay
        DelaunayTriangulationBuilder triangulator = new DelaunayTriangulationBuilder();
        triangulator.setSites(multiPoint);
        Polygon[] triangles = getTriangles(triangulator.getTriangles(geometryFactory1));
      for (Polygon triangle : triangles) {
            // Faites quelque chose avec le triangle, par exemple, accédez aux sommets du triangle
            Coordinate[] triangleCoordinates = triangle.getCoordinates();
            for (Coordinate vertex : triangleCoordinates) {
                double xVertex = vertex.getX();
                double yVertex = vertex.getY();

            }
        }
        double[]Z=listToArray(z);
        double[][] zi = linearInterpolate(triangles,Z, Xi, Yi);
    }

public static double[][] linearInterpolate(Polygon[] triangles, double[] z, double[][] Xi, double[][] Yi) {
    int n = Xi.length;
    int m = Xi[0].length;
    double[][] zi = new double[n][m];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double interpolatedX = Xi[i][j];
            double interpolatedY = Yi[i][j];

            // Trouver le triangle contenant les coordonnées interpolées
            for (int k = 0; k < triangles.length; k++) {
                Polygon triangle = triangles[k];
                if (triangle.contains(geometryFactory.createPoint(new Coordinate(interpolatedX, interpolatedY)))) {
                    // Calculer la valeur interpolée zi[i][j] à partir des coordonnées interpolées et du triangle
                    zi[i][j] = interpolateInsideTriangle(interpolatedX, interpolatedY, triangle, z, k);
                    break;
                }
            }
        }
    }

    return zi;
}
    public static double interpolateInsideTriangle(double x, double y, Polygon triangle, double[] z, int triangleIndex) {
        Coordinate[] vertices = triangle.getCoordinates();
        Coordinate v0 = vertices[0];
        Coordinate v1 = vertices[1];
        Coordinate v2 = vertices[2];

        // Calculer les coordonnées barycentriques du point (x, y) par rapport au triangle
        double denominator = ((v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y));
        double alpha = ((v1.y - v2.y) * (x - v2.x) + (v2.x - v1.x) * (y - v2.y)) / denominator;
        double beta = ((v2.y - v0.y) * (x - v2.x) + (v0.x - v2.x) * (y - v2.y)) / denominator;
        double gamma = 1.0 - alpha - beta;

        // Vérifier les indices avant d'accéder au tableau z
        int index0 = triangleIndex * 3;
        int index1 = index0 + 1;
        int index2 = index0 + 2;

        // Assurer que les indices sont valides
        if (index0 >= 0 && index0 < z.length && index1 >= 0 && index1 < z.length && index2 >= 0 && index2 < z.length) {
            // Utiliser les coordonnées barycentriques pour interpoler la valeur z à l'intérieur du triangle
            double interpolatedZ = alpha * z[index0] + beta * z[index1] + gamma * z[index2];
            return interpolatedZ;
        } else {
            // Si les indices ne sont pas valides, renvoyer une valeur par défaut ou gérer l'erreur selon votre besoin
            return 0.0; // Valeur par défaut (ajustez selon votre logique)
        }
    }



    public static double[] listToArray(ArrayList<Double> list) {
        double[] array = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }
  private static Polygon[] getTriangles(org.locationtech.jts.geom.Geometry trianglesGeometry) {
        List<Polygon> triangles = new ArrayList<>();
        for (int i = 0; i < trianglesGeometry.getNumGeometries(); i++) {
            if (trianglesGeometry.getGeometryN(i) instanceof Polygon) {
                triangles.add((Polygon) trianglesGeometry.getGeometryN(i));
            }
        }
        return triangles.toArray(new Polygon[0]);
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
