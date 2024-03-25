#include "functions.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

//Functions and things used throughout
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

//count the number of lines in a file
int count_lines( const char* filename){
    
    int count=0; //number of lines in the file

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
        overlap_ob.close();
    }
    return count;
}

int comp_ind(int u, int v){
    int comp_ind;
    comp_ind = u*(u+1)/2 + v;
    return comp_ind;
}


int main()
{
    //STEP 1: import the nuclear repulsion
    std::ifstream nuc_rep_ob;
    nuc_rep_ob.open("enuc.dat");

    double nuc_rep;
    if (nuc_rep_ob.is_open())
    {
        nuc_rep_ob >> nuc_rep;
        nuc_rep_ob.close();
    }

    //STEP 2: Read in and store one-electron integrals

    //Overlap integrals
    typedef std::vector<std::vector<double>> twoD_vec;
    
    //get number of lines in the file
    int count;
    count = count_lines("s.dat");
    
    //read the file into lists
    std::ifstream overlap_ob;
    overlap_ob.open("s.dat");

    //lists of data
    std::vector<int> index_n(count);
    std::vector<int> index_m(count);
    std::vector<double> integral(count);

    //populate lists
    if(overlap_ob.is_open()){    
        for(int i = 0; i<count; i++){
            overlap_ob >> index_n[i]; overlap_ob >> index_m[i]; overlap_ob >> integral[i];
        }
        overlap_ob.close();
    }

    //put values into overlap integral matrix
    //initialize 2D vector for matrix
    int num_basis = (count-1);
    int num_aos = index_m[num_basis];
    std::vector<double> layer_1(num_aos,0);
    twoD_vec AO_s(index_n[num_basis], layer_1);

    for(int i=0; i<count; i++){
        AO_s[index_n[i]-1][index_m[i]-1] = integral[i];
    }

    //Kinetic energy
    //read the file into lists
    std::ifstream kin_ob;
    kin_ob.open("t.dat");

    //lists of data
    std::vector<int> index_w(count);
    std::vector<int> index_r(count);
    std::vector<double> integral_w(count);

    //populate lists
    if(kin_ob.is_open()){    
        for(int i = 0; i<count; i++){
            kin_ob >> index_w[i]; kin_ob >> index_r[i]; kin_ob >> integral_w[i];
        }
        kin_ob.close();
    }

    //put values into overlap integral matrix
    //initialize 2D vector for matrix
    twoD_vec kin_t(index_n[num_basis], layer_1);

    for(int i=0; i<count; i++){
        kin_t[index_w[i]-1][index_r[i]-1] = integral_w[i];
    }

    //Nuclear Attraction Integrals
    //read the file into lists
    std::ifstream nuc_ob;
    nuc_ob.open("v.dat");

    //lists of data
    std::vector<int> index_q(count);
    std::vector<int> index_v(count);
    std::vector<double> integral_v(count);

    //populate lists
    if(nuc_ob.is_open()){    
        for(int i = 0; i<count; i++){
            nuc_ob >> index_q[i]; nuc_ob >> index_v[i]; nuc_ob >> integral_v[i];
        }

        nuc_ob.close();
    }

    //put values into overlap integral matrix
    //initialize 2D vector for matrix
    twoD_vec nuc_ac(index_n[num_basis], layer_1);

    for(int i=0; i<count; i++){
        nuc_ac[index_q[i]-1][index_v[i]-1] = integral_v[i];
    }

    //Form the Core Hamiltonian
    twoD_vec H_core_vec(index_n[num_basis], layer_1);

    for(int i =0; i<index_n[num_basis]; i++){
        for(int n=0; n < index_n[num_basis]; n++){
            H_core_vec[i][n] = kin_t[i][n] + nuc_ac[i][n];
        }
    }

    //turn core hamiltonian into a matrix

    Matrix H_core(num_aos, num_aos);

    for(int i=0; i<num_aos; i++){
        for(int n=0; n<num_aos; n++){
            H_core(i,n) = H_core_vec[i][n];
            H_core(n,i) = H_core_vec[i][n];
        }
    }

    //STEP 3: Two Electron Integrals
    //count number of lines in file
    int eri_count;
    eri_count = count_lines("eri.dat");

    //import data
    std::ifstream eri_ob;
    eri_ob.open("eri.dat");

    //lists of data
    std::vector<int> eri_1(eri_count);
    std::vector<int> eri_2(eri_count);
    std::vector<int> eri_3(eri_count);
    std::vector<int> eri_4(eri_count);
    std::vector<double> eri_integral(eri_count);

    //populate lists
    if(eri_ob.is_open()){    
        for(int i = 0; i<eri_count; i++){
            eri_ob >> eri_1[i]; eri_ob >> eri_2[i]; eri_ob >> eri_3[i]; eri_ob >> eri_4[i]; eri_ob >> eri_integral[i];
        }
        eri_ob.close();
    }
    
    //initialize final array of eri integrals
    int eri_num = comp_ind(comp_ind(num_aos-1, num_aos-1), comp_ind(num_aos-1, num_aos-1));
    std::vector<double> eri(eri_num+1);

    //fill final array of eri integrals using the lists imported from the file
    for(int i = 0; i<eri_count; i++){
        //get indices
        int u = eri_1[i]-1; int v = eri_2[i]-1; int l = eri_3[i]-1; int o = eri_4[i]-1;
        int uv = (u>v) ? comp_ind(u, v) : comp_ind(v, u);
        int lo = (l>o) ? comp_ind(l, o) : comp_ind(o, l);
        int uvlo = (uv>lo) ? comp_ind(uv, lo) : comp_ind(lo, uv);

        eri[uvlo] = eri_integral[i];
    }

    //STEP 4 Build the orthogonalization matrix
    Matrix AO_s_matrix(num_aos, num_aos);

    for(int i=0; i<num_aos; i++){
        for(int n=0; n<num_aos; n++){
            AO_s_matrix(i,n) = AO_s[i][n];
            AO_s_matrix(n,i) = AO_s[i][n];
        }
    }

    Eigen::SelfAdjointEigenSolver<Matrix> solver(AO_s_matrix);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();

    //I sorta have the diagonalized S matrix, s, now I need s^(-1/2), so I need the inverse square root of the eigenvalues
    Matrix s_inv_root = Matrix::Zero(num_aos, num_aos); //initialize matrix with 0s and fill in the eigenvalues

    //fill in s^(-1/2)
    for(int i=0; i<num_aos; i++){
        s_inv_root(i,i) = 1/sqrt(evals(i));
    }

    //undiagonalize s^-(1/2) to get S^-(1/2) using the eigenvectors
    Matrix orth_mat = (evecs * s_inv_root * evecs.transpose()); 

    //STEP 5 Build the Initial Guess Density

    //initial fock matrix in orthonormal AO basis
    Matrix fock_init = (orth_mat.transpose() * H_core * orth_mat);

    //trouble shooting: making all values smaller than 1e-10 actually 0
    for (int i =0; i<num_aos; i++){
        for(int n = 0; n<num_aos; n++){
            fock_init(i,n) = (abs(fock_init(i,n))<1e-10) ? 0 : fock_init(i,n);
        }
    }

    //diagonalize fock matrix
    Eigen::SelfAdjointEigenSolver<Matrix> fock_solver(fock_init);
    Matrix coef_1_ortho = fock_solver.eigenvectors();
    Matrix energy_1 = fock_solver.eigenvalues();

    //transform eigenvectors into the original AO basis
    Matrix coef_1 = (orth_mat * coef_1_ortho);

    //trouble shooting: making all values smaller than 1e-10 actually 0
    for (int i =0; i<num_aos; i++){
        for(int n = 0; n<num_aos; n++){
            coef_1(i,n) = (abs(coef_1(i,n))<1e-10) ? 0 : coef_1(i,n);
        }
    }

    //Building the density matrix using the occupied MOs
    Matrix density_1(num_aos, num_aos);

    for(int i=0; i<num_aos; i++){
        for(int n=0; n<num_aos; n++){
            double d_in=0;
            for(int m = 0; m<5; m++){
                d_in += coef_1(i,m)*coef_1(n,m);
            }
            density_1(i,n) = d_in;
        }
    }    

    //STEP 6 Compute the Initial SCF Energy
    double ee_1=0;
    
    for (int i =0; i<num_aos; i++){
        for(int n = 0; n<num_aos; n++){
            ee_1 += density_1(i,n)*(H_core(i,n)+H_core(i,n));
        }
    }

    std::cout << "Initial Electronic energy and total energy: "<< ee_1 << "   "<< ee_1+ nuc_rep << std::endl;

    /*SCF LOOP START*/


    //STEP 7 Compute the New Fock Matrix

    //initialize loop variables
    Matrix F(num_aos, num_aos);
    Matrix C;
    Matrix D(num_aos, num_aos);
    Matrix D_prev = density_1;
    double ee_prev = ee_1;
    double ee;
    int iter = 1;

    //initialize loop termination conditions
    double delta_E = 10;
    double rms_D = 10;

    while(delta_E > 10e-12 && rms_D > 10e-11){

        for(int u = 0; u<num_aos; u++){
            for(int v = 0; v<num_aos; v++){
                double f_uv = H_core(u,v);
                double G =0;
                //loops to calculate G from the previous iteration's density matrix

                for(int l = 0; l<num_aos; l++){
                    for(int o=0; o<num_aos; o++){
                        //get (uv|lo) index
                        int uv = (u>v) ? comp_ind(u, v) : comp_ind(v, u);
                        int lo = (l>o) ? comp_ind(l, o) : comp_ind(o, l);
                        int uvlo = (uv>lo) ? comp_ind(uv, lo) : comp_ind(lo, uv); //eri[uvlo] = (uv|lo)

                        //get (ul|vo) index
                        int ul = (u>l) ? comp_ind(u, l) : comp_ind(l, u);
                        int vo = (v>o) ? comp_ind(v, o) : comp_ind(o, v);
                        int ulvo = (ul>vo) ? comp_ind(ul, vo) : comp_ind(vo, ul); //eri[ulvo] = (ul|vo)

                        G += D_prev(l,o)*(2*eri[uvlo]-eri[ulvo]);               

                    }
                }

                f_uv += G;
                F(u,v) =(f_uv);
                F(v,u) = (f_uv);
            }
        }

        //std::cout << "\nThe new fock matrix is: \n" << F << std::endl;

        //STEP 8 Build the New Density Matrix

        Matrix F_orth = (orth_mat.transpose() * F * orth_mat);

        //std::cout << "\nSecond Orthogonal Fock Matrix:\n" << F_orth << std::endl;

        //diagonalize fock matrix
        Eigen::SelfAdjointEigenSolver<Matrix> F_solver(F_orth);
        Matrix coef_2_ortho = F_solver.eigenvectors();
        Matrix energy_2 = F_solver.eigenvalues();

        //std::cout << "\nSecond Orthogonal MO Matrix:\n" << coef_2_ortho << std::endl;

        //transform eigenvectors into the original AO basis
        C = (orth_mat * coef_2_ortho);

        //std::cout << "\nSecond MO Matrix:\n" << C << std::endl;

        //Building the density matrix using the occupied MOs

        for(int i=0; i<num_aos; i++){
            for(int n=0; n<num_aos; n++){
                double d_in=0;
                for(int m = 0; m<5; m++){
                    d_in += C(i,m)*C(n,m);
                }
                D(i,n) = d_in;
            }
        }  

        //trouble shooting: manually correcting the one wrong value
        //D(4,4) = 1;

        //std::cout << "\n Second Density Matrix:\n" << D << std::endl;

        //STEP 9 Computer New SCF Energy

        //std::cout << H_core(2,3) << ' ' << fock_init(2,3) << std::endl;
        
        double ee_add = 0;
        for (int i =0; i<num_aos; i++){
            for(int n = 0; n<num_aos; n++){
                ee_add += D(i,n)*(H_core(i,n)+F(i,n));
            }
        }

        ee = ee_add;

        std::cout <<"\nIter " << iter << " energy and total energy: "<< ee << "   "<< ee + nuc_rep;

        //check termination conditions and update iteration
        iter += 1;
        delta_E = abs(ee-ee_prev);
        ee_prev = ee; //set current ee to previous ee for next iteration

        double rms_D_add=0;
        for(int u=0; u<num_aos; u++){
            for(int v=0; v<num_aos; v++){
                rms_D_add += pow((D(u,v)- D_prev(u,v)),2);
            }
        }

        rms_D = pow(rms_D_add, 0.5);
        D_prev = D; //set current D to previous D for next iteration

        std::cout << " delta ee: " << delta_E << " rms_D: "<< rms_D << std::endl;

    }


    return 0;
}