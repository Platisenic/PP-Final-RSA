#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
#include "CycleTimer.h"
#include "matrix.h"
#include <string>
#include <algorithm>
#include <mpi.h>
using namespace std;
const static uint32_t MINIMAL_BOUND = 300;

// Sieving interval length.
const static uint32_t INTERVAL_LENGTH = 65536;//65536;
std::vector<uint32_t> generateFactorBase(const mpz_class& N, uint32_t B) {
    std::vector<uint32_t> factorBase;
    /*
     * This is essentially the sieve of Eratosthenes, with the addition
     * that it also checks for (N/p) = 1.
     */
    std::vector<bool> sieve(B + 1, false);
    for (uint32_t p = 2; p <= B; ++p) {
        if (sieve[p])
            continue;

        // Add p to factor base if N is a quadratic residue modulo p.
        if (mpz_legendre(N.get_mpz_t(), mpz_class(p).get_mpz_t()) == 1)
            factorBase.push_back(p);

        // Add multiples of p to sieve.
        for (uint32_t i = p; i <= B; i += p)
            sieve[i] = true;
    }

    return factorBase;
}

void QS(string str){
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Status status;
    mpz_class N;
    int N_c;
    char* buf;
    vector<uint32_t> smooth;
    vector<vector<uint32_t>> smoothFactors;

    if(world_rank==0){
        // cin>>N;
        // str = N.get_str();
        buf = strdup(str.data());
        N_c = strlen(buf);
        // cout<<"nc: "<<N_c<<endl;
        for(int i=1;i<world_size;i++){
            MPI_Send(&N_c,1,MPI_INT,i,0,MPI_COMM_WORLD);
            MPI_Send(buf,N_c,MPI_CHAR,i,0,MPI_COMM_WORLD); //dest,tag
        }
        N = buf;
        printf("%s\n",buf);

    }else{
        MPI_Recv(&N_c,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
        buf = (char*)calloc(N_c,sizeof(char));
        MPI_Recv(buf,N_c,MPI_CHAR,0,0,MPI_COMM_WORLD,&status);
        N = buf;
        // cout<<"rank : "<<world_rank<<"N "<<N<<endl;
    }
    const float logN = mpz_sizeinbase(N.get_mpz_t(), 2) * std::log(2); // Approx.
    const float loglogN = std::log(logN);
    const mpz_class sqrtN = sqrt(N);

    // Smoothness bound B.
    const uint32_t B = MINIMAL_BOUND + std::ceil(std::exp(0.55*std::sqrt(logN * loglogN)));
    vector<uint32_t> factorBase = generateFactorBase(N,B);
    int* count;
    // printf("rank: %d size : %lu\n",world_rank,factorBase.size());
    count = (int*)calloc(world_size,sizeof(int));
    uint32_t intervalStart = world_rank*INTERVAL_LENGTH;
    int endSign = 0;
    int sumCount =0;
    int tims=0;
    uint32_t* smoothNum;
    uint32_t* factorNum;
    // printf("rank %d intervalStart %u\n",world_rank,intervalStart);
    if(world_rank==0){
        // cout<<"rank 0"<<factorBase.size()<<endl;
        int nowCount = 0;
        while(sumCount<factorBase.size()+20){
            // printf("rank %d sumCount :%d\n",world_rank,sumCount);
            uint32_t a,i;
            for(i=0,a=intervalStart;i<INTERVAL_LENGTH;i++,a++){
                mpz_class z_square = ((a+sqrtN)*(a+sqrtN))-N; 
                vector<uint32_t> factors;
                for(uint32_t j=0;j<factorBase.size();j++){
                    const int p = factorBase[j];
                    while(mpz_divisible_ui_p(z_square.get_mpz_t(), p)){
                        mpz_divexact_ui(z_square.get_mpz_t(), z_square.get_mpz_t(), p);
                        factors.push_back(j);
                    }
                }
                if(z_square ==1){
                    smooth.push_back(a); //need change
                    smoothFactors.push_back(factors);
                };
            }
            count[world_rank] = smooth.size();
            sumCount =count[world_rank];
            for(int i=1;i<world_size;i++){
                MPI_Recv(&count[i],1,MPI_INT,i,1,MPI_COMM_WORLD,&status); //tag 1
                sumCount += count[i];
            }
            if(sumCount>=factorBase.size()+20)  endSign = 1;
            for(int i=1;i<world_size;i++)  
                MPI_Send(&endSign,1,MPI_INT,i,2,MPI_COMM_WORLD);//tag 2
            // cout<<"rank "<<world_rank<<" intervalStart "<<intervalStart<<endl;            
            intervalStart += INTERVAL_LENGTH*world_size;
            // printf("rank 0 intervalStart %u\n",intervalStart);
        }
        // printf("rank %d sumCount %ld\n",world_rank,smooth.size());
        // receive smooth size , soomth member;
        int smoothSize=0; //client smooth size
        for(int i=1;i<world_size;i++){
            MPI_Recv(&smoothSize,1,MPI_INT,i,3,MPI_COMM_WORLD,&status);
            smoothNum = (unsigned int*)calloc(smoothSize,sizeof(unsigned int));
            MPI_Recv(smoothNum,smoothSize,MPI_UNSIGNED,i,4,MPI_COMM_WORLD,&status);
            for(int j=0;j<smoothSize;j++){
                // if(smooth.size()<factorBase.size()+20) 
                smooth.push_back(smoothNum[j]);
            }
            int factorMaxSize=0;
            MPI_Recv(&factorMaxSize,1,MPI_INT,i,5,MPI_COMM_WORLD,&status);//recv factormax size
            // printf("recv rank %d smoothSize %d factorMaxsize %d\n",i,smoothSize,factorMaxSize);
            factorNum = (unsigned int*)calloc(factorMaxSize,sizeof(unsigned int)); //
            for(int k=0;k<smoothSize;k++){ //receive smooth factors
                memset(factorNum,0,factorMaxSize);
                MPI_Recv(factorNum,factorMaxSize,MPI_UNSIGNED,i,6,MPI_COMM_WORLD,&status);//recv all factors
                vector<uint32_t> factors;
                bool non_zero = false;
                int stop_size= 0 ;
                for(int s=0;s<factorMaxSize;s++){
                    if(!factorNum[s] && non_zero) break;
                    factors.push_back(factorNum[s]);
                    if(factorNum[s])non_zero = true;
                }
                // if(smoothFactors.size()<factorBase.size()+20)
                smoothFactors.push_back(factors);
                // printf("recv rank %d factorNum %d th size %lu\n",i,k,factors.size());
            }
        }
        // printf("all smooth Size %ld factorSize %ld\n",smooth.size(),smoothFactors.size());
        // printf("rank %d sumCount %ld\n",world_rank,smooth.size());
    }else{
        while(!endSign){
            uint32_t a,i;
            for(i=0,a=intervalStart;i<INTERVAL_LENGTH;i++,a++){
                mpz_class z_sqaure = ((a+sqrtN)*(a+sqrtN))%N;
                vector<uint32_t> factors;
                for(uint32_t j=0;j<factorBase.size();j++){
                    const int p = factorBase[j];
                    while(mpz_divisible_ui_p(z_sqaure.get_mpz_t(),p)){
                        mpz_divexact_ui(z_sqaure.get_mpz_t(),z_sqaure.get_mpz_t(),p);
                        factors.push_back(j);
                    }
                }
                if(z_sqaure==1){
                    smooth.push_back(a);
                    smoothFactors.push_back(factors);
                }
            }
            count[world_rank] = smooth.size();
            MPI_Send(&count[world_rank],1,MPI_INT,0,1,MPI_COMM_WORLD);
            MPI_Recv(&endSign,1,MPI_INT,0,2,MPI_COMM_WORLD,&status);
            // cout<<"rank "<<world_rank<<" intervalStart "<<intervalStart<<endl;            
            intervalStart += INTERVAL_LENGTH*world_size;   
        }
        int smoothSize=0;
        for(auto i:smoothFactors){
            if(i.size()>smoothSize) smoothSize = i.size();
        }
        // printf("rank %d intervalStart %u endSign %d smoothSize %ld smoothfactorsize %d\n",
        //     world_rank,intervalStart,endSign,smooth.size(),smoothSize);
        int smoothLen = smooth.size();
        MPI_Send(&smoothLen,1,MPI_INT,0,3,MPI_COMM_WORLD);
        smoothNum = (unsigned int*)calloc(smoothLen,sizeof(unsigned int));
        for(int i=0;i<smoothLen;i++) smoothNum[i] = smooth[i];
        MPI_Send(smoothNum,smoothLen,MPI_UNSIGNED,0,4,MPI_COMM_WORLD); //smooth number
        // printf("send rank %d smoothSize %d\n",world_rank,smoothSize);
        MPI_Send(&smoothSize,1,MPI_INT,0,5,MPI_COMM_WORLD);//send factorMax size
        factorNum  = (unsigned int*)calloc(smoothSize,sizeof(unsigned int)); //factorMax size
        // printf("rank %d  smoothsize %d\n",world_rank,smoothLen,smoothFactors.size());
        for(int i=0;i<smoothLen;i++){ // all smooth
            memset(factorNum,0,smoothSize);
            for(int j=0;j<smoothSize;j++){
                if(j<smoothFactors[i].size()) factorNum[j] = smoothFactors[i][j];
                else factorNum[j] = 0;
                // cout<<factorNum[j]<<" ";
            }
            // cout<<endl;
            MPI_Send(factorNum,smoothSize,MPI_UNSIGNED,0,6,MPI_COMM_WORLD);//send factors
            // printf("Send rank %d factorNum %d th size %lu maxSize %d\n",world_rank,i,smoothFactors[i].size(),smoothSize);
        }
        
    }
    if(world_rank==0){
        // cout<<"smooth size"<<smooth.size()<<" "<<smoothFactors.size()<<endl;
        Matrix M(factorBase.size(),smoothFactors.size()+1);
        for(int i=0;i<smoothFactors.size();i++){
            // cout<<"smooth factors : "<<smoothFactors[i].size()<<endl;
            for(int j=0;j<smoothFactors[i].size();j++){
                M(smoothFactors[i][j],i).flip();//
            }
        }
        M.reduce();
        mpz_class a;
        mpz_class b;
        do {
            vector<uint32_t> x = M.solve();
            // for(auto ss:x){
            //     cout<<ss<<" ";
            // }
            // cout<<endl;
            a = 1;
            b = 1;
            vector<uint32_t> decomp(factorBase.size(),0);
            for(uint32_t i=0;i<smoothFactors.size();i++){
                if(x[i]==1){
                    for(uint32_t p=0;p<smoothFactors[i].size();p++) ++decomp[smoothFactors[i][p]];
                    b *= (smooth[i]+sqrtN);
                }
            }
            for(uint32_t p=0;p<factorBase.size();++p){
                for(uint32_t i=0;i<(decomp[p]/2);++i){
                    a *= factorBase[p];
                }
            }
        }while(a % N == b % N || a % N == (- b) % N + N);
        
        mpz_class factor = gcd(b-a,N);
        cout<<factor<<" "<<N/factor<<" "<<factor*(N/factor)<<endl;

    }


}
int main(int argc,char** argv){
    MPI_Init(NULL,NULL);
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    // double start = CycleTimer::currentSeconds();
    double start = MPI_Wtime();
    string buf="";
    if(world_rank==0){
        buf = argv[1];
    }
    QS(buf);
    // double conduct_time = CycleTimer::currentSeconds()-start;
    double conduct_time = MPI_Wtime()-start;
    printf("time : %f\n",conduct_time);
    MPI_Finalize();
    return 0;
}