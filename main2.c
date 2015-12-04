/**
 * Aufgabe 1
 *
 * Student: Jonathan Pirnay
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LRmP.h"

/**
 * Teilaufgabe i.)
 */

double Bisektionsverfahren(double (*f)(double x), double a, double b, int maxiter, double tol, FILE *file) {
	// Setze Startwerte
	double x_minus = a;
	double x_plus = b;
	int numiter = 0; // Zähler für Iterationen
	double root = NAN;

	while (numiter < 100) {
		numiter++;
		double x_middle = (x_minus + x_plus)/2.0;
		double func_value = f(x_middle);
		
		if (fabs(func_value) < tol) {
			root = x_middle;
			break;
		}
		else if (func_value < 0) {
			x_minus = x_middle;
		}
		else if (func_value > 0) {
			x_plus = x_middle;
		}

		fprintf(file, "%f \n", fabs(func_value));
	}

	return root;
}

double Newtonverfahren(double (*f)(double x), double (*df)(double x), double x0, int maxiter, double tol, FILE *file) {
	int numiter = 0; // Zähler für Iterationen
	double root = NAN;

	while (numiter < 100) {
		numiter++;
		double func_val = f(x0);
		double df_val = df(x0);
		if (fabs(func_val) < tol) {
			root = x0;
			break;
		}
		else {
			if (df_val < tol) {
				// Ableitung = 0; Abbruch
				break;
			}
			else {
				x0 = x0 - (1/df_val)*func_val;
			}
		}

		fprintf(file, "%f \n", fabs(func_val));
	}

	return root;
}

// Testfunktionen für Teilaufgabe 1
double f1(double x) {
	return (x-2.0) * x;
}

double f2(double x) {
	return pow(x-2.0, 5);
}

double df1(double x) {
	return 2*x-2; 
}

double df2(double x) {
	return 5 * pow(x-2.0, 4); 
}

/**
 *
 * Teilaufgabe ii.)
 * 
 */

// Hilfsfunktion, um euklidische Norm zu bestimmen
// Parameter: Dimension, Vektor
double euklid_norm(int n, double *vector) {
	double norm = 0.0;
	for (int i=0; i<n; i++) {
		norm += pow(vector[i], 2);
	}

	return sqrt(norm);
}

// Hilfsfunktion, die Matrix aufbaut.
double** makeMatrix(int n) {
	double** matrix = malloc(sizeof(double) * n * n);
	for (int i=0; i<n; i++) {
		matrix[i] = malloc(sizeof(double) * n);
	}

	return matrix;
}

void freeMatrix(double **matrix, int n) {
	for (int i=0; i<n; i++) {
		free(matrix[i]);
	}

	free(matrix);
}

// Hilfsfunktion von Skalarmultiplikation von Vektoren
void scalarMult(double *vector, double scalar, int n) {
	for (int i=0; i<n; i++) {
		vector[i] = vector[i] * scalar;
	}
}

// Hilfsfunktion, die den zweiten Vektor auf den ersten addiert. Dabei werden die Werte des ersten überschrieben
void addVectors(double *a, double *b, int n) {
	for (int i=0; i<n; i++) {
		a[i] += b[i];
	}
}

int Newton_Multidim(int N, double (*f)(double *, double *, int), double (*df)(double *, double **, int), double *x0, double *fx, int maxiter, double tol) {
	int numiter = 0;

	// Differentialmatrix aufbauen
	double **jacobi = makeMatrix(N);	

	while (numiter < maxiter) {
		numiter++;

		// Überprüfe zunächst, ob der aktuelle Funktionswert schon unsere Nullstellenschwelle unterschreitet
		f(x0, fx, N);
		
		if (euklid_norm(N, fx) < tol) {
			// Ist unterschritten. Gib 0 (= Erfolg) zurück
			freeMatrix(jacobi, N);
			return 0;
		}
		else {
			// f(x) wurde bereits berechnet. Berechne nun den neuen Wert von x, also x^k+1
			df(x0, jacobi, N); 			// Zunächst Differential am Punkt x^k
			scalarMult(fx, -1, N);		// f(x^k) * -1
			Solve(N, jacobi, fx, 1); 	// Löse GLS auf f(x^k). Die Lösung wird wie gehabt auf fx geschrieben
			addVectors(x0, fx, N);		// x^k+1 setzen (= x0 + Lösung von verherigem Schritt)
		}
	}

	freeMatrix(jacobi, N);
	return 1;
}

// Testfunktionen für Teilaufgabe ii.)
double f(double *x, double *fx, int n) {
	fx[0] = (1.0/10.0)*pow(x[0],2) + sin(x[1]);
	fx[1] = cos(x[0]) + (1.0/10.0)*pow(x[1],2);

	return 0;
}

double df(double *x, double **jacobi, int n) {
	jacobi[0][0] = 0.2 * x[0];
	jacobi[0][1] = cos(x[1]);
	jacobi[1][0] = -1*sin(x[0]);
	jacobi[1][1] = 0.2 * x[1];

	return 0;
}

// Hilfsfunktion für das Testen der Newtonfunktion
void testNewtonMultiDim(double *x0, double *fx, double x, double y) {
	x0[0] = x;
	x0[1] = y;
	printf("\n Startwert x0 = (%.1f, %.1f)\n", x, y);
	double success = Newton_Multidim(2, f, df, x0, fx, 100, 1.0e-7);
	if (success < 1.0e-7) {
		printf("Nullstelle gefunden am Punkt (%f, %f) \n", x0[0], x0[1]);
	}
	else {
		printf("Keine Nullstelle gefunden.\n");
	}
}

int main(void) {
	// Teilaufgabe 1: Bisektionsverfahren
	FILE *file_bisektion_f1;
	FILE *file_bisektion_f2;
	file_bisektion_f1 = fopen("bisektion_f1", "w");
	file_bisektion_f2 = fopen("bisektion_f2", "w");
	printf("Teilaufgabe i.) \n\n Bisektionsverfahren \n");
	double rf1 = Bisektionsverfahren(f1, 1.5, 100, 100, 1.0e-7, file_bisektion_f1);
	double rf2 = Bisektionsverfahren(f2, 1.5, 100, 100, 1.0e-7, file_bisektion_f2);
	printf("Nullstelle f1: %f \n", rf1);
	printf("Nullstelle f2: %f \n", rf2);
	fclose(file_bisektion_f1);
	fclose(file_bisektion_f2);
	//Teilaufgabe 1: Newtonverfahren
	FILE *file_newton_f1;
	FILE *file_newton_f2;
	file_newton_f1 = fopen("newton_f1", "w");
	file_newton_f2 = fopen("newton_f2", "w");
	printf("Teilaufgabe i.) \n\n Newtonverfahren \n");
	double rnf1 = Newtonverfahren(f1, df1, 100, 100, 1.0e-7, file_newton_f1);
	double rnf2 = Newtonverfahren(f2, df2, 100, 100, 1.0e-7, file_newton_f2);
	printf("Nullstelle f1: %f \n", rnf1);
	printf("Nullstelle f2: %f \n", rnf2);
	fclose(file_newton_f1);
	fclose(file_newton_f2);

	//Teilaufgabe 2: Newtonverfahren im mehrdimensionalen
	double *x0 = malloc(sizeof(double) * 2);
	double *fx = malloc(sizeof(double) * 2);

	printf("Newtonverfahren im Mehrdimensionalen. \n");
	testNewtonMultiDim(x0, fx, 0.5, 0.5);
	testNewtonMultiDim(x0, fx, -0.5, -0.5);
	testNewtonMultiDim(x0, fx, 0.5, 1.5);
}