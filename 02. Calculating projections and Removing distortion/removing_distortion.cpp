#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <cmath>

// #include <Python.h>
// #include "matplotlibcpp.h"

// namespace plt = matplotlibcpp;
using namespace std;
using namespace Eigen;



void ucitajTacke(vector<pair<float, float>> &tacke, int duzina) {
	pair<float, float> tacka;
	for (int i = 0; i < duzina; i++) {
		cin >> tacka.first >> tacka.second;
		tacke.push_back(tacka);
	}
}


MatrixXf kreirajMatricuIx3(vector<pair<float, float>> &tacke) {
	MatrixXf Mat(tacke.size(), 3);

	int i = 0;
	while (i < tacke.size()) {
		Mat(i, 0) = tacke[i].first;
		Mat(i, 1) = tacke[i].second;
		Mat(i, 2) = 1;
		i++;
	}
	/*
		| t[0][p1] t[0][p2] 1 |
		| t[1][p1] t[1][p2] 1 |
		| t[2][p1] t[3][p2] 1 | 
		| t[4][p1] t[4][p2] 1 |
		| ........			  |
	*/
	return Mat;
}


// Svaka tacka i njena slika iz vektora pravi ovu matricu 2 x 9
// Za broj tacaka n, pravi matricu 2n x 9
MatrixXf kreirajMatricu2nx9(vector<pair<float, float>> &pocetne, vector<pair<float, float>> &projektovane) {
	
	MatrixXf matrix2nx9(2 * pocetne.size(), 9);
	
	auto i = pocetne.begin();
	auto j = projektovane.begin();
	int k = 0;
	
	while (i < pocetne.end()) {
		matrix2nx9(k, 0) = 0;
		matrix2nx9(k, 1) = 0;
		matrix2nx9(k, 2) = 0;
		matrix2nx9(k, 3) = -i->first;
		matrix2nx9(k, 4) = -i->second;
		matrix2nx9(k, 5) = -1;
		matrix2nx9(k, 6) = j->second * i->first;
		matrix2nx9(k, 7) = j->second * i->second;
		matrix2nx9(k, 8) = j->second;

		matrix2nx9(k + 1, 0) = i->first;
		matrix2nx9(k + 1, 1) = i->second;
		matrix2nx9(k + 1, 2) = 1;
		matrix2nx9(k + 1, 3) = 0;
		matrix2nx9(k + 1, 4) = 0;
		matrix2nx9(k + 1, 5) = 0;
		matrix2nx9(k + 1, 6) = -(j->first * i->first);
		matrix2nx9(k + 1, 7) = -(j->first * i->second);
		matrix2nx9(k + 1, 8) = -(j->first);
	
		i++; 
		j++; 
		k += 2;
		/// |  0	 0    0  -b[0] -b[1] -1   (p[1]*b[0])   (p[1]*b[1])   p[1]	|		b - par pocetnih koordinata tacaka
		///	| p[0]  p[1]  1   0     0     0  -(p[0]*b[0])  -(p[0]*b[1])  -p[0]	|		p - par projektovanih koordinata tacaka
		///	|	.	.	.	.	.	.	.	.	.	.	.	.	.	.		.	|
		/// |	.	.	.	.	.	.	.	.	.	.	.	.	.	.		.	|
		/// |	.	.	.	.	.	.	.	.	.	.	.	.	.	.		.	|
	}

	return matrix2nx9;
}

MatrixXf resiP(MatrixXf Mat, Vector3f vec, int k) {
	Vector3f stepen = ((Mat.transpose()).inverse()) * vec;

	int m = stepen.size();
	int n = k;
	MatrixXf P(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			P(i, j) = stepen(i) * Mat(i, j);
		}
	}
	return P.transpose();
}


MatrixXf naivniAlgoritam(vector<pair<float, float>> pocetne, vector<pair<float, float>> projektovane) {

	float x, y;

	x = pocetne.back().first;
	y = pocetne.back().second;
	Vector3f v(x, y, 1);

	x = projektovane.back().first;
	y = projektovane.back().second;
	Vector3f vp(x, y, 1);

	pocetne.pop_back();
	projektovane.pop_back();

	MatrixXf M = kreirajMatricuIx3(pocetne);
	MatrixXf Mp = kreirajMatricuIx3(projektovane);

	MatrixXf P1 = resiP(M, v, pocetne.size());
	MatrixXf P2 = resiP(Mp, vp, projektovane.size());


	return P2 * P1.inverse();
}


MatrixXf DLTAlgoritam(vector<pair<float, float>> pocetne, vector<pair<float, float>> projektovane) {

	int br_tacaka = pocetne.size();

	MatrixXf matrix2nx9 = kreirajMatricu2nx9(pocetne, projektovane);
	// SVD razbija matricu na proizvode 3 matrice
	JacobiSVD<MatrixXf> Svd(matrix2nx9, ComputeFullV);
	
	MatrixXf Mat(2 * br_tacaka, 9);
	MatrixXf matrica_poslednje_kolone(3, 3);
	Mat = Svd.matrixV();	// Uzimamo 3. matricu V iz A = USV*

	int br_kolona = Mat.cols();
	// Uzimamo poslednju kolonu iz dekompozicije
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			matrica_poslednje_kolone(i, j) = Mat(i * 3 + j, br_kolona - 1);
		}
	}



	return matrica_poslednje_kolone;
}

float izrSrednjeRastojanje(vector<pair<float, float>> tacke, pair<float, float> tezite) {

	float prosecno_rastojanje = 0;
	int n = tacke.size();
	for (auto i = tacke.begin(); i < tacke.end(); i++) {
		float medju_rastojanje = sqrt(pow(i->first - tezite.first, 2) + pow(i->second - tezite.second, 2));	// sqrt((x2-x1)^2 + (y2-y1)^2)
		prosecno_rastojanje += medju_rastojanje;
	}

	return prosecno_rastojanje / n;
}

Vector3f izrCentroidu(vector<pair<float, float>> tacke) {

	float tx = 0;
	float ty = 0;
	int n = tacke.size();

	// Sabiramo tezine koordinata
	for (auto i = tacke.begin(); i < tacke.end(); i++) {
		tx += i->first;
		ty += i->second;
	}

	// centroida = (tx/n, ty/n, 1)
	// homogene
	Vector3f centroida;
	centroida << tx / n, ty / n, 1;

	/*
	// Test: teziste
	cout << "Teziste je : " << centroida;
	cout << endl;
	*/

	return centroida;
}

MatrixXf izrMatricuNormalizacije(vector<pair<float, float>> tacke) {
	/*
	1. Izracunati teziste sistema tacaka
	2. Translirati teziste u koordinatni pocetak (matrica T)
	3. Skalirati tacke tako da prosecna udaljenost tacke od koordinatnog pocetka bude sqrt(2) (matrica homotetije S)
	4. Matrica normalizacije je jednaka S * T
	*/
	Vector3f vteziste = izrCentroidu(tacke);
	pair<float, float> t(vteziste(0), vteziste(1));

	// Matrica T: Translira teziste ovog skupa tacaka u O(0,0)
	MatrixXf T(3, 3);
	T << 1, 0, -vteziste(0),
		 0, 1, -vteziste(1),
		 0, 0,	    1;
	// cout << T << endl;
	float srednje_rastojanje = izrSrednjeRastojanje(tacke, t);
	float k_skaliranja = sqrt(2) / srednje_rastojanje;
	// Matrica S: Skalira, tako da prosecno rastojanje tacaka od koordinatnog pocetka bude jednak sqrt(2)
	MatrixXf S(3, 3);
	S << k_skaliranja,		0,			0,
			0,			k_skaliranja,	0,
			0,				0,			1;
	// cout << S << endl;
	// Prvo transliramo pa skaliramo
	return S * T;
}

MatrixXf nDLTAlgoritam(vector<pair<float, float>> pocetne, vector<pair<float, float>> projektovane) {
	/*
	1. Normalizuj tacke i slike:	nx[i] = T * x[i], npx[i] = Tp * px[i]
	2. Uradi DLT na:				nx[i] <-> npx[i]
	3. Denormalizuj resenje:		H = Tp_inv * nH * T
	*/
	
	MatrixXf T = izrMatricuNormalizacije(pocetne);
	MatrixXf s_tacka = (T * kreirajMatricuIx3(pocetne).transpose()).transpose();	// 3x3 * nx3.transp
	//	cout << T << endl << endl << s_tacka << endl << endl;

	MatrixXf Tp = izrMatricuNormalizacije(projektovane);
	MatrixXf p_tacka = (Tp * kreirajMatricuIx3(projektovane).transpose()).transpose();	// 3x3 * nx3.transp
	//	cout << Tp << endl << endl << p_tacka << endl << endl;


	int n = pocetne.size();
	vector<pair<float, float>> s;
	for (int i = 0; i < n; i++) {
		s.push_back(pair<float, float>(s_tacka(i, 0), s_tacka(i, 1)));
		//	cout << s_tacka(i, 0) << " " << s_tacka(i, 1) << endl;
	}
	
	vector<pair<float, float>> p;
	for (int i = 0; i < n; i++) {
		p.push_back(pair<float, float>(p_tacka(i, 0), p_tacka(i, 1)));
		//	cout << p_tacka(i, 0) << " " << p_tacka(i, 1) << endl;
	}

	/*
	////	---- Plot ----
	for (int i = 0; i < n; i++) {
		plt::plot(s_tacka(i,0), s_tacka(i, 1), 'go');
		plt::plot(p_tacka(i, 0), p_tacka(i, 1), 'g+');
	}
	*/

	MatrixXf H = DLTAlgoritam(s, p);
	return Tp.inverse() * H * T;
}

int main() {
	int opcija;
	cout << "Izaberite opciju: " << endl << "1. Rucni unos koordinata" << endl << "2. Izabir koordinata pomocu misa" << endl << ":::  ";
	cin >> opcija;
	if (opcija == 1) {				// Rucni unos koordinata
		// Unos tacaka
		int br_tacaka;
		cout << "Unesite broj korespodencija: ";
		cin >> br_tacaka;
		if (br_tacaka < 0) {
			cout << "greska: Neispravan unos! (ocekivan broj veci od 0)" << endl;
			return -1;
		
		}

		vector<pair<float, float>> pocetne_tacke;
		vector<pair<float, float>> projektovane_tacke;

		cout << "Ucitajte ulazne tacke kao parove x, y: " << endl;
		ucitajTacke(pocetne_tacke, br_tacaka);
		cout << "Ucitajte izlazne tacke kao parove x, y: " << endl;
		ucitajTacke(projektovane_tacke, br_tacaka);
		cout << endl;
		if (br_tacaka == 4) {
			cout << "--- Naivni algoritam: ---" << endl;
			cout << naivniAlgoritam(pocetne_tacke, projektovane_tacke) << endl;
			
			/*
			////	---- Plot ----
			for (int i = 0; i < br_tacaka; i++) {
				plt::plot(pocetne_tacke[i].first, pocetne_tacke[i].second, 'ro');
				plt::plot(projektovane_tacke[i].first, projektovane_tacke[i].second, 'r+');
			}
			*/

		}

		cout << "--- DLT algoritam: ---" << endl;
		cout << DLTAlgoritam(pocetne_tacke, projektovane_tacke) << endl;

		cout << "---Normalizovani DLT algoritam: ---" << endl;
		cout << nDLTAlgoritam(pocetne_tacke, projektovane_tacke);

		//	ERROR: Error	C1083	Cannot open include file : 'Python.h' : No such file or directory
		//	VS ne cita Python header iako je ubacen u Additional Include Libraries, kao i postaveljen kao header u src dir

		/*
		////	---- Plot ----
		plt::figure_size(1200, 780);
		plt::show();
		*/

	}
	else if (opcija == 2) {			// Izabir koordinata pomocu misa
		
		cout << "Nije zavrseno!";
		return 0;
		
	}
	else {							// Greska
		cout << "greska: Izabrana je nepostojeca opcija!";
		return -1;
	}


	return 0;
}