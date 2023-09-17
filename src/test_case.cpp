#include "dl.h"
using namespace dl;
int main()
{
    double t1 = clock();
    dl::Data in_dat, out_dat; 
    //genData("test");
    ReadCSV("test_in", in_dat);  
    ReadCSV("test_out", out_dat);
    

    std::cout << "Final_mse= " << "\n";
    return 0;
}


