/**
 * Aufgabe 1
 *
 * Student: Jonathan Pirnay
 * 
 */

#include <stdio.h>
#include <stdlib.h>

// Testmatrizen
// 
// LR-Zerlegung von test_matrix sollte L = {{1,0,0},{1,1,0},{3,3,1}}, R = {{1,2,3},{0,-1,-2},{0,0,-2}}
// 
int n = 3;
double test_matrix_working_nopivot[3][3] = {
	{1.0,2.0,3.0},
	{1.0,1.0,1.0},
	{3.0,3.0,1.0}
};

// Funktioniert nicht und sollte im zweiten Schritt abbrechen
double test_matrix_not_working_nopivot[3][3] = {
	{1,2,3},
	{1,2,3},
	{0,0,1}
};

// LR-Zerlegung mit Pivotisierung sollte sein: {{2,-4,-6},{0.5,1,6},{-1,0,-1}} // Testvektor sollte sein: {1,3,2}
double test_matrix_working_pivot[3][3] = {
	{2,-4,-6},
	{-2,4,5},
	{1,-1,3}
};

// Funktioniert nicht und sollte im zweiten Schritt abbrechen
double test_matrix_notworking_pivot[3][3] = {
	{2,-4,-6},
	{-4,8,12},
	{0,0,3}
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

// Hilfsfunktion, die Betrag zurückgibt
double absolute(double a) {
	if (a < 0) {
		return -a;
	}
	else {
		return a;
	}
}

// Ballert direkt auf die Matrix, daher kein Rückgabewert
// "k" ist Nummer des Schritts
void GaussSpaltenelimination(int n, int k, double **A) {
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

			GaussSpaltenelimination(n, k, A);	
		}

	}
	// LR-Zerlegung mit Pivotisierung
	else if (flag == 1) {
		// Permutationsvektor mit normaler Reihenfolge initialisieren
		for (int k=0; k<n; k++) {
			s[k] = k;
		}

		for (int k=0; k < n-1; k++) {
			int p = k;

			for (int i=k+1; i < n; i++) {
				// Wir suchen jetzt das Pivotelement
				if (absolute(A[i][k]) > absolute(A[p][k])) {
					p = i;
				}
				if (p != k) {
					// Im Permutationsvektor vertauschen wir jetzt p und k. Merke: Im Permutationsvektor steht an der k-ten Stelle
					// welche (alte) Zeile jetzt k-te Zeile ist
					int mem_p = s[p];
					s[p] = s[k];
					s[k] = mem_p;

					// Wir vertauschen jetzt die tatsächlichen Zeilen (k und p) in der Matrix
					for (int j=0; j<n; j++) {
						double mem = A[k][j];
						A[k][j] = A[p][j];
						A[p][j] = mem;
					}
				}
			}
			
			// Hier könnte es vorkommen, dass wir unser Pivotelement 0 ist (Matrix nicht invertierbar). Wie oben abbrechen
			if (A[k][k] == 0) {
				brokeAtStep = k+1;
				break;
			}
			
			printf("Schritt: %d, Pivotzeile: %d \n", k+1, p+1);
			// Yo, alles vertauscht, jetzt Standardpivotisierung
			GaussSpaltenelimination(n, k, A);
		}
	}

	// Ausgabe
	if (brokeAtStep > 0) {
		printf("LR-Zerlegung abgebrochen an Punkt %d \n", brokeAtStep);
	}
	else {
		// Alles hat funktioniert. Ausgabe je nach Modus
		if (flag == 0) {
			printf("LR-Zerlegung ohne Pivotisierung liefert Matrix: \n");
			displayMatrix(A, n);
		}
		else if (flag == 1) {
			printf("LR-Zerlegung mit Pivotisierung liefert Matrix: \n");
			displayMatrix(A, n);
			printf("Mit Permutationsvektor: \n");
			for (int i=0; i<n; i++) {
				printf("%5d", s[i] + 1);	
			}
			printf("\n");
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

	//Funktionierende Matrix mit Pivotisierung
	printf("---------- \n");
	int *s = malloc(sizeof(int) * n);
	LR(n, makeMatrix(test_matrix_working_pivot), s, 1);

	//Nicht funktionierende Matrix mit Pivotisierung. Sollte beim zweiten Schritt abbrechen.
	printf("---------- \n");
	LR(n, makeMatrix(test_matrix_notworking_pivot), s, 1);
}