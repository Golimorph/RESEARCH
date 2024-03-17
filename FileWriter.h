//
//  FileWriter.h
//  PyrochloreG
//
//  Created by Richard Edberg on 2018-02-02.
//

#ifndef FileWriter_h
#define FileWriter_h

#include <stdio.h>
#include <string.h>
#include "Eigen/Dense"
#include <ios>
#include <fstream>
#include "Parameters.h"
#include <iomanip>


class FileWriter {
public:
    FileWriter(Parameters *parameters);

    void append(std::string delimeterLine, Eigen::MatrixXd *data);
    
    void writeDataMatrix(Eigen::MatrixXd *data, long long int elapsed_secs, int size);
    

    
private:
    Parameters *parameters;
    std::string timeString;
};
    
    
    
    
    
#endif /* FileWriter_h */


