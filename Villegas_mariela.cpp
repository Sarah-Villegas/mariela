// Team Members: Sarah Villegas, Ali Fenton, Chanel Aquino
// Date Created: 30 October 2015
// Description: Speech Signal Analysis
/*
// 1) Table of calculated values are in comparison.txt file.
// 2) Which values for each data file match most closely?
    Mean, average power, and standard deviation values match most closely.
// 3) Which values for each data file are most different?
    Zero crossings, average magnitude, and variance values are most different.
// 4) Are there other statistical measures that you could suggest?
// 5) Can you determine if these sound recordings are from the same person? Explain.
*/

#include <iostream> // cin, cout
#include <fstream>  // file IO
#include <cmath>    // sqrt(), abs()
#include <iomanip>  // setw()
#include <cstdlib>  // exit()
#include <vector>   
using namespace std;

double mean(vector<double> vct);
double stdDev(double vct); 
double var(vector<double> vct , double& mean);  
double averagePower(vector<double> vct);
double avgMagnitude(vector<double> vct); 
int zeroCross(vector<double> vct);


//*******************************
int main()
{

    // create vectors of type double
    vector<double> sound1;
    vector<double> sound2;
    
    int 
        index = 0,
        size;
        
   double 
        meanDifference, 
        zeroCrossDifference, 
        magnitudeDifference, 
        powerDifference, 
        varDifference, 
        stdDifference,
        power1,
        power2,
        zero1, 
        zero2, 
        mag1, 
        mag2, 
        std1, 
        std2,
        value,
        value1,
        mean1,
        mean2,
        standD1,
        standD2;
    
    // create input/out files
    ifstream fin1, fin2;
    ofstream fout;
    
    // open files
    fin1.open("two_a.txt");
    fin2.open("two_b.txt");
    fout.open("comparison.txt");
    
    // exit program if fail
    if (fin1.fail() || fin2.fail() || fout.fail())
    {
        cout << "ERROR";
        exit(1);
    }    

    // read input into fin1
    while(fin1 >> value)
    {
        sound1.push_back(value);
    }
    
    // read input into fin2
    while(fin2 >> value1){
        sound2.push_back(value1);
    }
    
    mean1 = mean(sound1);
    mean2 = mean(sound2);
    standD1 = var(sound1, mean1);
    standD2 = var(sound2, mean2);
    power1 = averagePower(sound1);
    power2 = averagePower(sound2);
    zero1 = zeroCross(sound1);
    zero2 = zeroCross(sound2);
    mag1 = avgMagnitude(sound1);
    mag2 = avgMagnitude(sound2);
    std1 = stdDev(standD1);
    std2 = stdDev(standD2);
    
    powerDifference = (abs(power1 - power2)/ ((power1 + power2)/2)) * 100.0;
    meanDifference = (abs(mean1 - mean2)/ ((mean1 + mean2)/2)) * 100.0;
    zeroCrossDifference = (abs(zero1 - zero2)/ ((zero1 + zero2)/2)) * 100.0;
    magnitudeDifference = (abs(mag1 - mag2)/ ((mag1 + mag2)/2)) * 100.0;
    varDifference = (abs(standD1 - standD2)/ ((standD1 + standD2)/2)) * 100.0;
    stdDifference = (abs(std1 - std2)/ ((std1 + std2)/2)) * 100.0;
    
    fout 
        << "Team Members: Sarah Villegas, Ali Fenton, Chanel Aquino\n" << endl
        << setw(30) << "two_a.txt" << setw(20) << "two_b.txt" << setw(20) << "% difference" << endl
        << "Mean" << setw(27) << mean1 << setw(20) << mean2 << setw(15) << meanDifference << endl
        << "Zero Cross" << setw(14) << zero1 << setw(19) << zero2 << setw(23) << zeroCrossDifference <<endl
        << "Average Power" << setw(16) << power1 << setw(19) << power2 << setw(17) << powerDifference << endl
        << "Average Magnitude"  << setw(12) << mag1  << setw(20) << mag2 << setw(17) << magnitudeDifference << endl
        << "Standard Deviation" << setw(12) << std1 << setw(20) << std2 << setw(15) << stdDifference << endl
        << "Variance" << setw(24) << var(sound1, mean1) << setw(20) << var(sound2, mean2) << setw(14) << varDifference << endl;
   
   
       
    // close files
    fin1.close();
    fin2.close();
    fout.close();
    return 0;
}

//*********************************
double mean(vector<double> vct)
{
    double 
        average,
        sum = 0.0,
        size = vct.size();
    
    for (int i = 0; i < size; i++)
    {
        sum += vct[i];
    }
    
    average = sum / size;
    
    return average;
    
}

//*********************************
double stdDev(double vct)
{
    return sqrt(vct);
}

//*********************************
double var(vector<double> vct, double& mean)
{
    
    double size = vct.size();
    double results = 0; 
    double squared = 0;  
    
    
    for(int i = 0; i < size; i++ ){
    
        squared = pow(vct[1] - mean ,2);
        results += squared;
   
    }
    results = results / size;

    return results;
}

//*********************************
double averagePower(vector<double> vct)
{
    double
        avgPower,
        sum = 0.0, 
        size = vct.size();
    
    for(int k = 0; k < size; k++)
    {
        vct[k] *= vct[k];
        sum += vct[k];
    }
    
    avgPower = sum / size;
    
    return avgPower;
}

//*********************************
double avgMagnitude(vector<double> vct)
{
    double avgMag, sum = 0.0, size = vct.size();
    
    for(int k = 0; k < size; k++){
        vct[k] = abs(vct[k]);
        sum += vct[k];
    }
    
    avgMag = sum / size;
    
    return avgMag;
}

//*********************************
int zeroCross(vector<double> vct){

    int totalZero = 0; 
    int size = vct.size(); 
    for(int i = 1; i < size; i++ ){
    
        if((vct[i - 1]  > 0) && (vct[i ] < 0)){
    
            totalZero++;
        }
        else if((vct[i - 1]  < 0) && (vct[i ] > 0)){
    
            totalZero++;
        }
    }

    return totalZero; 
}






