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
double test_matrix[3][3] = {
	{2.0,1.0,1.0},
	{4.0,2.0,-1.0},
	{-1,0.0,7.0}
};

/**
 *
 * Aufgabe i.)
 * 
 */

// Hilfsfunktionen, die eine Matrix/Vektor "schön" ausgibt
void displayMatrix(double **matrix, int n) {
	printf("\n");
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			printf("%5.1f", matrix[i][j]);
		}
		printf("\n");
	}
}

void displayVector(double *vector, int n) {
	printf("\n ( ");
	for (int i=0; i<n; i++) {
		printf("%5.1f", vector[i]);
	}
	printf(" ) \n");
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

/**
 *
 * Aufgabe ii.)
 * 
 */

int VwSubs(int n, double **L, double *b) {
	// erster Eintrag von b bleibt, da wir annehmen können, dass L normiert
	for (int i=1; i<n; i++) {
		// i bezeichnet die Zeile. Gehe jetzt die Spalten durch
		double mem = 0;
		for (int j=0; j<i; j++) {
			mem = mem + (b[j] * L[i][j]);
		}
		b[i] = b[i] - mem;
	}

	return 0;
}

int RwSubs(int n, double **R, double *b) {
	// ähnlich wie VwSubs, nur von hinten
	// wir müssen allerdings den ersten Eintrag berücksichtigen, weil nicht gegeben ist, dass R normiert
	for (int i=n-1; i>=0; i--) {
		double mem = 0;
		for (int j=n-1; j>i; j--) {
			mem = mem + (R[i][j] * b[j]);
		}
		if (R[i][i] == 0) {
			return 1;
		}
		b[i] = (b[i] - mem) / R[i][i];
	}

	return 0;
}

/**
 *
 * Aufgabe iii.)
 * 
 */

int Solve(int n, double **A, double *b, int flag) {
	// Leeren Permutationsvektor initialisieren
	int *s = malloc(sizeof(int) * n);

	int lr = LR(n, A, s, flag);

	if (lr > 0) {
		// abgebrochen bei LR Zerlegung
		printf("'Solve' abgebrochen bei LR-Zerlegung. Schritt: %d \n", lr);
	}
	else {
		// wenn mit Pivotisierung, dann müssen wir auch die Einträge von b permutieren, mache das mit einem Hilfsvektor
		if (flag == 1) {
			double *mem_b = malloc(sizeof(double) * n);
			for (int i=0; i<n; i++) {
				mem_b[i] = b[s[i]];
			}
			for (int i=0; i<n; i++) {
				// rüberschieben
				b[i] = mem_b[i];
			}
			free(mem_b);
		}
		
		VwSubs(n, A, b);
		int err = RwSubs(n, A, b);

		if (err) {
			printf("Fehler bei RwSubs");
		}
		else {
			return 0;
		}
	}

	return 1;
}

//-------

// Hilfsfunktion, die aus der oberen hartcodierten 3x3 Testmatrizen dynamische Matrix macht
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

	double *b = malloc(sizeof(double) * 3);
	b[0] = 1.0;
	b[1] = 2.0;
	b[2] = 3.0;

	int err = Solve(3, makeMatrix(test_matrix), b, 1);
	if (!err) {
		printf("Lösung: \n");
	}
	displayVector(b, 3);

}