/**
 * Aufgabe 1
 *
 * Student: Jonathan Pirnay
 * 
 */

#include <stdio.h>
#include <stdlib.h>

// Testdaten
// 
// LR-Zerlegung von test_matrix sollte L = {{1,0,0},{1,1,0},{3,3,1}}, R = {{1,2,3},{0,-1,-2},{0,0,-2}}
// 
int n = 3;
double test_matrix_working_nopivot[3][3] = {
	{1.0,2.0,3.0},
	{1.0,1.0,1.0},
	{3.0,3.0,1.0}
};
double test_matrix_not_working_nopivot[3][3] = {
	{1,2,3},
	{1,2,3},
	{0,0,1}
};

/**
 *
 * Aufgabe i.)
 * 
 */

// Hilfsfunktion, die eine Matrix "schön" ausgibt
void displayMatrix(double **matrix, int n) {
	printf("\n");
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			printf("%5.1f", matrix[i][j]);
		}
		printf("\n");
	}
}

// HAUTPFUNKTION VON i.) Gibt bei Erfolg "0" zurück.
int LR (int n, double **A, int *s, int flag) {

	// Zeigt an, ob die Zerlegung an einem bestimmten Punkt abbrechen musste. 
	// Wenn ja (Wert >0), wird später einfach der Wert anstatt der Matrix ausgegeben
	int brokeAtStep = 0;

	// LR-Zerlegung ohne Pivotisierung
	if (flag == 0) {
		for (int k = 0; k < n-1; k++) {
			// Diagonalelement ist 0. Breche hier ab.
			if (A[k][k] == 0) {
				brokeAtStep = k + 1;
				break;
			}

			for (int i = k+1; i < n; i++) {
				// Hier setzen wir den "Eliminationsfaktor" für den Eintrag a_i_k, aus dem später die L-Matrix bestehen wird.
				double factor = A[i][k] / A[k][k];
				A[i][k] = factor;

				// Zeile updaten und schön abziehen alter
				for (int j = k+1; j<n; j++) {
					A[i][j] = A[i][j] - (factor * A[k][j]);
				}
			}
		}

	}
	// LR-Zerlegung mit Pivotisierung
	else if (flag == 1) {

	}

	// Ausgabe
	if (brokeAtStep > 0) {
		printf("LR-Zerlegung abgebrochen an Punkt %d \n", brokeAtStep);
	}
	else {
		if (flag == 0) {
			printf("LR-Zerlegung ohne Pivotisierung liefert Matrix: \n");
			displayMatrix(A, n);
		}
	}
	
	return brokeAtStep;
}

// Hilfsfunktion, die aus den oberen hartcodierten 3x3 Testmatrizen dynamische Matrizen macht
double **makeMatrix(double m[][3]) {
	double **matrix = malloc(sizeof(double) * 9);
	for (int i=0; i<3;i++) {
		matrix[i] = malloc(sizeof(double) * 3);
	}

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			matrix[i][j] = m[i][j];
		}
	}

	return matrix;
}

//--------

int main(void) {
	//Funktionierende Matrix ohne Pivotisierung
	printf("---------- \n");
	LR(n, makeMatrix(test_matrix_working_nopivot), NULL, 0);

	//Nicht funktionierende Matrix ohne Pivotisierung. Sollte beim zweiten Schritt abbrechen.
	printf("---------- \n");
	LR(n, makeMatrix(test_matrix_not_working_nopivot), NULL, 0);
}