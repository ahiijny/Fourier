#ifndef DFT_H
#define DFT_H

#include <complex>
#include <vector>

using namespace std;
typedef complex<double> dcomp;

class Dft
{
    public:
        vector<double> fdomain(double time_span, int N);
        vector<dcomp> dft(vector<dcomp> &x, int N);
        vector<dcomp> fft(vector<dcomp> &x, int N);
        vector<dcomp> fourierCoeffs(vector<dcomp> &DFT, int N);
}
;

#endif // DFT_H
