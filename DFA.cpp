#include "DFA.h"

using namespace std;
using namespace Eigen;

void DFAm (const double *x, int x_length, const unsigned int *scale,
        int scale_length, const double *q, int q_length, 
        unsigned int m, void* res_ptr)
{
    Map<const ArrayXd> xvec(x, x_length);
    
    ArrayXXd DFA_result(q_length, scale_length);

    for (unsigned int i = 0; i < scale_length; ++i)
    {
        const unsigned int w = scale[i];
		if (w <= m+1 || w > x_length) {
            throw scaleException();
        }
        
		// Number of segments with this window size
        const unsigned int Ns = x_length / w;
        
        // Structure function array
        ArrayXd F2(2*Ns);
        
		// QR decomposition is same for all windows, since only dependent on window size
		// Declare indices for better polyfit and create vandermonde matrix
		const VectorXd idx = VectorXd::LinSpaced(w, -w/2, w/2);
		MatrixXd Vdm( w, m+1 );
		// Vdm = [idx^m, idx^(m-1), idx^(m-2), ..., idx^0]
		Vdm.col(m).setOnes();
		for (int j = m-1; j >=0; --j) {
			Vdm.col(j) = (idx.array() * Vdm.col(j+1).array()).matrix();
		}    

		const HouseholderQR<MatrixXd> QR (Vdm);

        // Loop through the windows (forwards and backwards)
        for (unsigned int v = 0; v < Ns; ++v) {
            // Calculate starting points of the blocks
            const unsigned int fwd_from = v*w;
            const unsigned int bw_from = x_length - (v+1)*w;
            
            // Calculate window variance after detrending
            const VectorXd x1 = xvec.segment(fwd_from, w).matrix();
            const VectorXd poly_coeff1 = QR.solve( x1 );
            F2[v] = (x1 - Vdm * poly_coeff1).array().square().sum() / w;
            
            const VectorXd x2 = xvec.segment(bw_from, w).matrix();
            const VectorXd poly_coeff2 = QR.solve( x2 );
            F2[v+Ns] = (x2 - Vdm * poly_coeff2).array().square().sum() / w;
        }
        
        // 
        for (int qc = 0; qc < q_length; ++qc) {
            if ( q[qc] != 0 ) {
                DFA_result(qc, i) = pow( F2.pow(q[qc]/2).mean(), 1/q[qc] );
            }
            else {
                DFA_result(qc, i) = exp( 0.5* F2.log().mean() );
            }
        }
    }
    // Copy result to new pointer
    memcpy (res_ptr, DFA_result.data(), sizeof(double)*q_length*scale_length);
}
    
                
        
        
                    

