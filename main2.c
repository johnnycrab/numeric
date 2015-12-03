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
}