#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include "Dft.h"

using namespace std;
typedef complex<double> dcomp;

const dcomp I(0,1);
const double PI = atan(1) * 4;

/** Returns the x-values for the domain of the DFT that 
 * is specified by the given time span and sample quantity. 
 */
vector<double> Dft::fdomain(double time_span, int N)
{
    vector<double> x(N/2 + 1, 0);
    double fmax = N / (2 * time_span);
    double df = 1 / time_span;
    
    for (int i = 0; i < (N/2 + 1); i++)
        x[i] = i * df;
        
    return x;    
}

/** Discrete Fourier Transform: 
 * Returns the DFT of the given set "x" containing N
 * periodic-discrete signals. Uses the naive DFT
 * algorithm, so complexity is O(n^2).
 */ 
vector<dcomp> Dft::dft(vector<dcomp> &x, int N)
{
    vector<dcomp> F(N, 0);
    for (int k = 0; k < N; k++)
        for (int i = 0; i < N; i++)
            F[k] = F[k] + x[i] * exp(-2*PI*k*i/N * I);
            
    return F;    
}

/** Fast Fourier Transform:
 * The return is exactly the same is dft(x,N),
 * but this algorithm has complexity O(N log_2(N)).
 */
vector<dcomp> Dft::fft(vector<dcomp> &x, int N)
{
    if (N == 1)
        return x;
    else if (N % 2 != 0) // FFT can't handle an off number of terms
        return dft(x,N);
    else
    {
        // Split even and odd indexed terms into 2 smaller sets
        int half = N/2;
        dcomp w;
        vector<dcomp> evens(half, 0);
        vector<dcomp> odds(half, 0);
        for (int k = 0; k < half; k++)
        {        
            evens[k] = x[2*k];
            odds[k] = x[2*k + 1];
        }        
        evens = fft(evens, half);           // Recursive call here
        odds = fft(odds, half);             // and here.
        for (int k = 0; k < half; k++)
        {
            w = exp(-2*PI*k/N * I);
            x[k] = evens[k] + w * odds[k];          // Recombine here
            x[k + half] = evens[k] - w * odds[k];    // and here.
        }
        return x;        
    }    
}

/** Fourier Coeffs:
 * Returns the Fourier coefficients for the Fourier series
 * for the time domain of the given discrete Fourier transform.
 * There are N/2 + 1 coefficients, where N is the number of data
 * points in the DFt. Re(Coeff) is the amplitude of the cosine 
 * for the corresponding frequency, whereas -Im(Coeff) is the
 * amplitude of the sine for the corresponding frequency.
 */
vector<dcomp> Dft::fourierCoeffs(vector<dcomp> &DFT, int N)
{
    int size = N/2 + 1;
    vector<dcomp> CF(size, 0);
    
    for (int i = 0; i < size; i++)
        CF[i] = 2.0/N * DFT[i];    
        
    CF[0] /= 2.0;
    CF[size - 1] /= 2.0;
    
    return CF;
}
