#ifndef __DFA_H__
#define __DFA_H__

#include <exception>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <cmath>

void DFAm (const double *x, int x_length, const unsigned int *scale,
        int scale_length, const double *q, int q_length, 
        unsigned int m, void *res_ptr);       
    
class scaleException: public std::exception
{
public:
    scaleException() {};
    virtual const char* what() const throw() {
        return "Invalid scales. Check that all scales m+1 < scale < x_length";
    }
};

#endif
