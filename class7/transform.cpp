#include "lagrangian.hpp"

void transform(double* x, double* y, int nk, int ne)
{
	double** x1 = new double*[nk+1];
	double** y1 = new double*[nk+1];
	for (int i = 0; i <= nk; i++)
	{
		x1[i] = new double[ne+1];
		y1[i] = new double[ne+1];
	}

	for (int i = 0; i <= nk; i++)
	{
		for (int j = 0; j <= ne; j++)
		{
		  //==============================================
			x1[i][j] = x[0] * (nk - i) * (ne - j) / (nk * ne) + x[1] * i * (ne - j) / (nk * ne) + x[2] * i * j / (nk * ne) + x[3] * (nk - i) * j / (nk * ne);
			y1[i][j] = y[0] * (nk - i) * (ne - j) / (nk * ne) + y[1] * i * (ne - j) / (nk * ne) + y[2] * i * j / (nk * ne) + y[3] * (nk - i) * j / (nk * ne);

            //==============================================
		}
	}

	ofstream outfile;
	outfile.open("grid.dat",ios::out);
	outfile<<"Variables=x,y,T" << endl;
	outfile << "Zone, I=" << nk + 1 << ", J=" << ne + 1 << endl;
	for(int j=0;j<=ne;j++)
		for (int i = 0; i <= nk; i++)
		{
			outfile << x1[i][j] << " " << y1[i][j] << " " << 1.0 << endl;
		}
	outfile.close();

	for (int i = 0; i <= nk; i++)
	{
		delete [] x1[i] ;
		delete [] y1[i];
	}
	delete [] x1;
	delete [] y1;
}