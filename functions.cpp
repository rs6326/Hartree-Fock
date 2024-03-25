#include "functions.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>


//count the number of lines in a file
int count_lines( const char* filename){
    
    int count=0 //number of lines in the file

    //open file
    std::ifstream overlap_ob;
    overlap_ob.open(filename);

    if(overlap_ob.is_open()){

        //count the number of lines in the file
        std::string line;
        while(overlap_ob.eof() == 0){
            getline(overlap_ob, line);
            count++;
        }
        std::cout << count << std::endl;
        overlap_ob.close();
    }
    return count;
}
