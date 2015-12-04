/**
 * Aufgabe 1
 *
 * Student: Jonathan Pirnay
 * 
 */

#include <stdio.h>
#include <stdlib.h>

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

		// Zeile updaten und schön abziehen
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
		for (int k = 0; k < n; k++) {
			// Diagonalelement ist 0. Breche hier ab.
			if (A[k][k] == 0) {
				brokeAtStep = k + 1;
				break;
			}

			if (k < n-1) {
				GaussSpaltenelimination(n, k, A);	
			}
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
				if (absolute(A[k][i]) > absolute(A[k][p])) {
					p = i;
				}
				if (p != k) {
					// Im Permutationsvektor vertauschen wir jetzt p und k. Merke: Im Permutationsvektor steht an der k-ten Stelle
					// welche (alte) Spalte jetzt k-te Spalte ist
					
					int mem_p = s[p];
					s[p] = s[k];
					s[k] = mem_p;

					// Wir vertauschen jetzt die tatsächlichen Spalten (k und p) in der Matrix
					for (int j=0; j<n; j++) {
						double mem = A[j][k];
						A[j][k] = A[j][p];
						A[j][p] = mem;
					}
				}
			}
			
			// Hier könnte es vorkommen, dass wir unser Pivotelement 0 ist (Matrix nicht invertierbar). Wie oben abbrechen
			if (A[k][k] == 0) {
				brokeAtStep = k+1;
				break;
			}
			
			// Yo, alles vertauscht, jetzt Standardeliminierung
			GaussSpaltenelimination(n, k, A);
		}
		if (brokeAtStep == 0 && A[n-1][n-1] == 0) {
			brokeAtStep = n;
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
		
	}
	else {
		
		VwSubs(n, A, b);
		int err = RwSubs(n, A, b);

		if (err) {
			
		}
		else {
			// wenn mit Pivotisierung, dann müssen wir auch die Einträge der Lösung gemäß P^-1 permutieren, mache das mit einem Hilfsvektor
			if (flag == 1) {
				double *mem_b = malloc(sizeof(double) * n);
				for (int i=0; i<n; i++) {
					mem_b[s[i]] = b[i];
				}
				for (int i=0; i<n; i++) {
					// rüberschieben
					b[i] = mem_b[i];
				}
				free(mem_b);
			}
			return 0;
		}
	}

	return 1;
}