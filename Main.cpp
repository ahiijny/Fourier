#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <conio.h>

#include "Dft.h"

using namespace std;
typedef complex<double> dcomp;

string initpath = "init.cfg";
string inpath = "data.csv";
string outpath = "fourier.csv";
bool error = false;
bool batch = false;
double peak_threshold = 0.003;
int t_col = 1;
int x_col = 2;
int row_1 = 1;
int row_2 = 8;
int small_N = 8;
int skip = 1;

/** Splits a string into a string array using the given delimiter.
 * Used to aid input parsing. Delimiter occurences are not included
 * in the returned array.
 */
vector<string> split(string str, string delim)
{
    // Declare Variables
    
    vector<string> params;
    string temp;
    size_t right = str.find(delim);
    
    // Iterate through string, finding occurrences of delimiter. 
    // Substring segment to array. Erase segment from source
    // string. Continue until no more delimiter occurrences.
    
    while (right != string::npos)
    {
        params.push_back(str.substr(0, right));
        str.erase(0, right + delim.size());
        right = str.find(delim);    
    }    
    params.push_back(str);    // Append the last remaining segment
    
    return params;
}

/** Writes a template initial conditions configuration file
 * to the inpath.
 */
bool writeTemplateInit()
{
    ofstream fout(initpath.c_str(), ios_base::out | ios_base::trunc);
    if (!fout.is_open())
        return false;
        
    // Write Template
    
    fout << "# === Configuration file for Fourier.exe ===" << endl;
    fout << "inpath = \"data.csv\"" << endl;
    fout << "outpath = \"fourier.csv\"" << endl;
    fout << "t column = 1" << endl;
    fout << "x column = 2" << endl;
    fout << "rows = 1 8" << endl;
    fout << "batch = false" << endl;
    fout << "N = 8" << endl;
    fout << endl;    
    cout << "Finished creating " << initpath << endl;
    fout.close();
    
    return true;    
}

/** Reads the configuration settings from file.
 */
bool loadInit()
{
    // Read File
    
    ifstream fin (initpath.c_str(), ios::in|ios::ate);
    if (fin.is_open())
    {       
        int size = fin.tellg();
        fin.seekg (0, ios::beg);
        vector<string>lines;    
        vector<string>params;    
        string line;
        while (getline(fin, line))        
            lines.push_back(line);               
        fin.close();
        cout << "Successfully read " << initpath << ": " << size << " bytes" << endl;
        
        // Parse Configuration Data        
        // Iterate through config file line by line
        
        int i = 0;
        
        do
        {
            params = split(lines[i], " ");
            
            if (params[0] == "inpath")            
                inpath = params[2].substr(1, params[2].size() - 2);
            
            else if (params[0] == "outpath")            
                outpath = params[2].substr(1, params[2].size() - 2);
            
            else if (params[0] == "t")            
                t_col = atoi(params[3].c_str());
            
            else if (params[0] == "x")            
                x_col = atoi(params[3].c_str());
            
            else if (params[0] == "rows")
            {            
                row_1 = atoi(params[2].c_str());
                row_2 = atoi(params[3].c_str());
            }
            
            else if (params[0] == "batch")            
                batch = params[2] == "true";
                
            else if (params[0] == "N")            
                small_N = atoi(params[2].c_str());
                
            else if (params[0] == "skip")
                skip = atoi(params[2].c_str());
                
            else if (params[0] == "peak")
                peak_threshold = atof(params[3].c_str());
            
        }
        while (++i < lines.size() && params.size() > 1);
        
        return true;
    }
    else
    {
        cout << "Error: could not read " << inpath << endl;
        return false;
    }
}

/** Loads the time domain data from the specified in path
 * and returns the data as a vector of vectors. The first
 * vector contains the t data. The second vector contains
 * the x data.
 */
vector< vector<double> > loadData()
{
    // Read File
    
    bool valid = true;
    double t, x;
    int N = row_2 - row_1 + 1, maxcol = max(t_col, x_col);
    vector< vector<double> > data(2, vector<double>(0));
    ifstream fin (inpath.c_str(), ios::in|ios::ate);
    if (fin.is_open() && N > 0 && t_col > 0 && x_col > 0)
    {       
        int size = fin.tellg();
        fin.seekg (0, ios::beg);        
        vector<string>lines;   
        vector<string>params;       
        string line;
        while (getline(fin, line))        
            lines.push_back(line);               
        fin.close();
        cout << "Successfully read " << inpath << ": " << size << " bytes" << endl;
        
        // Parse Input Data From Start Column to
        // Iterate through config file line by line              
        
        for (int i = row_1 - 1; i < row_2; i++)
        {
            valid = i < lines.size();
            if (valid)
            {            
                params = split(lines[i], ",");
                valid = params.size() >= maxcol;
            }
            if (valid)
            {            
                t = atof(params[t_col - 1].c_str());
                x = atof(params[x_col - 1].c_str());
                data[0].push_back(t);
                data[1].push_back(x);
            }
            else
            {
                cout << "Error: Could not read row " << (i+1) << endl;
                error = true;
                break;
            }
        }          
    }
    else
    {
        cout << "Error: could not read " << inpath << endl;
        error = true;
    }
    return data;
}

/** Outputs the Fourier transform values to file.
 */ 
bool writeFourier(vector<double> hz, vector<dcomp> F)
{
    ofstream fout(outpath.c_str(), ios_base::out | ios_base::trunc);

    // couldn't open it (disk error?); fail
    if (!fout.is_open())
    {    
        cout << "Could not write to " << outpath << endl;
        return false;
    }
    
    // Write Headers
    
    fout << "Frequency,Amplitude,,Real,Imag" << endl;
    fout << ",,,," << endl;
    
    for (int i = 0; i < hz.size(); i++)
    {
        fout << hz[i] << ",";
        fout << abs(F[i]) << ",,";
        fout << real(F[i]) << ",";
        fout << imag(F[i]) << endl;
        
    }
    
    cout << "Wrote " << hz.size() << " rows to " << outpath << endl;
    fout.close();
}

bool writeBatchFourier(double t1, double dt, vector< vector< vector <double> > > &all_peaks)
{
    ofstream fout(outpath.c_str(), ios_base::out | ios_base::trunc);

    // couldn't open it (disk error?); fail
    if (!fout.is_open())
    {    
        cout << "Could not write to " << outpath << endl;
        return false;
    }
    
    // Write Headers
    
    fout << "t1,t2,,Peak Frequency,Amplitude" << endl;
    fout << ",,,," << endl;
    
    for (int i = 0; i < all_peaks.size(); i++)
    {
        fout << t1 + i * skip * dt << "," << t1 + (i * skip + small_N - 1) * dt << ",,";
        for (int j = 0; j < all_peaks[i][0].size(); j++)
        {
            fout << all_peaks[i][0][j] << "," << all_peaks[i][1][j] << endl;
            fout << ",,,";
        }
        fout << "," << endl;
    }
    
    cout << "Wrote " << all_peaks.size() << " steps to " << outpath << endl;
    fout.close();
    
    return true;
}

vector< vector<double> > findPeaks(vector<double> &hz, vector<dcomp> &F)
{
    vector< vector<double> > peaks(2, vector<double>(0)); // First entry is frequency data; second entry is amplitude data
    vector<double> mags;
    double normalizer = abs(F[0]) * peak_threshold;
    
    for (int i = 0; i < hz.size(); i++) // Store absolute values
        mags.push_back(abs(F[i]));
        
    for (int i = 0; i < hz.size() - 4; i++) // Find peaks
    {
        double mag1 = mags[i];
        double mag2 = mags[i+1];
        double mag3 = mags[i+2];
        double mag4 = mags[i+3];
        double mag5 = mags[i+4];
        if (mag3 > normalizer) // Make sure peak is significant and not noise
        {   
            if (mag1 < mag3 && mag2 < mag3 && mag4 < mag3 && mag5 < mag3)
            {                     
                peaks[0].push_back(hz[i+2]);
                peaks[1].push_back(mag3);
            }
        }
    }
    return peaks;
}

int main()
{
    // Variables
    
    Dft dft;
    vector< vector<double> > data;
    vector<double> hz;
    vector<dcomp> x;
    vector<dcomp> F;
    double N, dt, time_span;
    int big_N;
    
    if (!loadInit())
    {    
        writeTemplateInit();
    }
    else
    {
        // Load Data
        
        data = loadData();
        if (!error)
        {       
            big_N = row_2 - row_1 + 1;
            N = big_N;
            if (batch)
                N = small_N;
            dt = data[0][1] - data[0][0];
            time_span = N * dt;
            cout << "Loaded data: " << endl;
            cout << "t_col = " << t_col << endl;
            cout << "x_col = " << x_col << endl;
            cout << "rows = " << row_1 << " " << row_2 << endl;
            cout << "N = " << N << endl;
            cout << "dt = " << dt << endl;
            cout << "time_span = " << time_span << endl;
            cout << "batch = " << (batch ? "true" : "false") << endl;
            if (batch)
            {            
                cout << "steps = " << big_N - N + 1 << endl;
                cout << "skip = " << skip << endl;
                cout << "peak threshold = " << peak_threshold << endl;
            }
            
            if (!batch) // normal
            {           
                // Convert x values to complex
                
                for (int i = 0; i < N; i++)
                    x.push_back(data[1][i]);
                
                // Compute the Fourier transform
                
                hz = dft.fdomain(time_span, N);
                F = dft.fft(x, N);
                F = dft.fourierCoeffs(F, N); // Scale magnitudes to amplitude coefficients
                
                // Output results
                
                writeFourier(hz, F);   
            }
            else // batch
            {
                vector< vector< vector<double> > > all_peaks;
                int steps = big_N - N + 1;
                int printInterval = 10 * skip;
                
                for (int i = 0; i < steps; i+=skip)
                {
                    x.clear();
                    
                    // Convert x values to complex
                    
                    for (int j = i; j < i + N; j++)
                        x.push_back(data[1][j]);
                        
                    // Compute the Fourier transform
                
                    hz = dft.fdomain(time_span, N);
                    F = dft.fft(x, N);
                    F = dft.fourierCoeffs(F, N); // Scale magnitudes to amplitude coefficients
                    
                    // Find peaks
                    
                    vector< vector<double> > temp = findPeaks(hz, F);
                    all_peaks.push_back(temp);
                    if (i % printInterval == 0)
                        cout << "Done " << i << endl;
                }
                writeBatchFourier(data[0][1], dt, all_peaks);                
            }
        }
    }
    
    cout << "Press any key to exit.";
    getch();
}
