package com.company;
import java.lang.*;
 
public class Main {
 
   private static final int DEFAULT_EQUATIONS_NUMBER = 2;
 
    public static void main(String[] args) {
 
	double eps = 10e-9;
	double x11 = 3;
	double x12 = 2;
	double x21 = 3;
	double x22 = -2;
	int k = 0;
 
	int kMax = 50;
	double[] a = {x11, x12, x21, x22};
 
 
 
        double[] deltaX = lsolve(matrix(a), vectSvob(a));
        System.out.println(k + ": " + deltaX[0] + "; " + deltaX[1]);
 
        if((deltaX[0] > eps) || (deltaX[1] > eps)) {
            do {
 
                a[0] += deltaX[0];
                a[2] += deltaX[0];
                a[1] += deltaX[1];
                a[3] += deltaX[1];
 
                deltaX = lsolve(matrix(a), vectSvob(a));
 
                k++;
                System.out.println(k + ": " + deltaX[0] + "; " + deltaX[1]);
                if (k == kMax) {
                    System.out.println("превышен лимит итераций!");
                    break;
                }
            } while ((deltaX[0] > eps) && (deltaX[1] > eps));
        }
 
  int k1 = k + 1;
        System.out.println("количество итераций: " + k1);
        System.out.println("Решение: ");
        for (int i = 0; i < DEFAULT_EQUATIONS_NUMBER; i++){
            System.out.print(deltaX[i] + " ");
        }
    }
 
    public static double func1(double x1, double x2){
        return (2*x1*x1 - x1*x2 - 5*x1 + 1);
    }
    public static double func2(double x1, double x2){
        return (x1 + 3*(Math.log10(x1)) - x2*x2);
    }
 
 
    public static double func1DiffX1(double x1, double x2){
        return (4*x1 - x2 - 5);
    }
    public static double func1DiffX2(double x1, double x2){
        return (-1*x1);
    }
 
 
    public static double func2DiffX1(double x1, double x2){
        return (1 + 1.3*x1);
    }
    public static double func2DiffX2(double x1, double x2){
        return (-2*x2);
    }
 
    public static double[][] matrix(double[] a){
 
        double[][] matrix = new double[DEFAULT_EQUATIONS_NUMBER][DEFAULT_EQUATIONS_NUMBER];
        matrix[0][0] = func1DiffX1(a[0], a[1]);
        matrix[0][1] = func1DiffX2(a[0], a[1]);
        matrix[1][0] = func2DiffX1(a[2], a[3]);
        matrix[1][1] = func2DiffX2(a[2], a[3]);
 
        return matrix;
    }
    public static double[] vectSvob(double[] a){
        double[] vectSvobb = new double[DEFAULT_EQUATIONS_NUMBER];
        vectSvobb[0] = -1 * func1(a[0], a[1]);
        vectSvobb[1] = -1 * func2(a[2], a[3]);
        return vectSvobb;
    }
 
 
    public static double[] lsolve(double[][] A, double[] b) {
 
        for (int p = 0; p < DEFAULT_EQUATIONS_NUMBER; p++) {
            int max = p;
            for (int i = p + 1; i < DEFAULT_EQUATIONS_NUMBER; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p];
            A[p] = A[max];
            A[max] = temp;
 
            double t = b[p];
            b[p] = b[max];
            b[max] = t;
 
            for (int i = p + 1; i < DEFAULT_EQUATIONS_NUMBER; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < DEFAULT_EQUATIONS_NUMBER; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }
 
        double[] x = new double[DEFAULT_EQUATIONS_NUMBER];
        for (int i = DEFAULT_EQUATIONS_NUMBER - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < DEFAULT_EQUATIONS_NUMBER; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }
}