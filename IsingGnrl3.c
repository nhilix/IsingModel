#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ARRAY 64
#define T_EQ 3
#define T_MEAS 8
#define T_RUN 1
#define MAXSIZE 150000
#define TCRIT 2.269
#define TDIV 0.1

double r,T,Tdiv,m_sum[MAXSIZE],m2_avg[T_RUN],m_avg[T_RUN], m2_runavg,m_runavg,prob[17];
int E,i,j,l,t,ileft,iright,jabove,jbelow,size, S[ARRAY][ARRAY],L;
Full_Run(int L,double m2_avg[],double m_avg[],int l);
Initialize(int S[ARRAY][ARRAY],int L);
Sweep(int S[ARRAY][ARRAY],double m_sum[],double prob[],int t,int L);
FILE* fout;

int main(void){
    fout=fopen("IsingPic.txt","w");
    srand(time(NULL));
    while(1){
    puts("L value(int):");
    scanf("%d",&L);
    if(L==1) break;
    Tdiv = TDIV;
    for(T=0.1;T<=2*TCRIT;T+=Tdiv){m2_runavg=0;m_runavg=0;
        for(E=-8;E<=8;E++)
            prob[E+8] = (double)1/(exp(E/T)+1);
        for(l=0;l<T_RUN;l++){m2_avg[l]=0;m_avg[l]=0;
        Full_Run(L,m2_avg,m_avg,l);
        m2_runavg += m2_avg[l]/T_RUN;
        m_runavg += m_avg[l]/T_RUN;}
        fprintf(fout,"%f \t%f \t%f\n",T,(double)m2_runavg/(T_MEAS*L*L),
                sqrt((m2_runavg/(T_MEAS*L*L) - m_runavg*m_runavg/(T_MEAS*L*L)/(T_MEAS*L*L))/(T_RUN*T_MEAS*L*L-1)));
        system("cls");
        for(j=0;j<L;j++){
        for(i=0;i<L;i++){
                if (S[i][j]<0) printf("%c%c",176,176);
                else printf("%c%c",178,178);}
        printf("\n");}
        printf("%2.0f%%",T/(2*TCRIT) * 100);
        if (T<(TCRIT-1) || T>(TCRIT+1)) Tdiv=TDIV;
        else Tdiv = (double) TDIV/10;}}
    fclose(fout);
    return 0;}
Sweep(int S[ARRAY][ARRAY],double m_sum[],double prob[],int t,int L)
{
    i=L-1;ileft=L-2;                //Use Heat Bath model to determine spin flip
    j=L-1;jabove=L-2;
    for(jbelow=0;jbelow<L;jbelow++){
        for(iright=0;iright<L;){
            E = 2*S[i][j]*(S[ileft][j]+S[iright][j]+S[i][jabove] + S[i][jbelow]);
            r = (double) rand() / RAND_MAX;
            if(r<prob[(int)E+8])
                S[i][j] *= -1;
            m_sum[t] += (double)S[i][j]/(L*L);
            ileft = i;
            i = iright++;}
        jabove = j;
        j = jbelow;}
}
Initialize(int S[ARRAY][ARRAY],int L){
    i=L-1;ileft=L-2;                  //Randomly initialize the spins of the Lattice
    j=L-1;jabove=L-2;
    for(jbelow=0;jbelow<L;jbelow++){
        for(iright=0;iright<L;){
            r = (double) rand()/RAND_MAX;
            if(r<0.5 && T>TCRIT)S[i][j] = 1;
            else S[i][j] = -1;
            ileft = i;
            i = iright++;}
        jabove = j;
        j = jbelow;}}
Full_Run(int L,double m2_avg[],double m_avg[],int l){
    Initialize(S,L);
    for(t=0;t<T_EQ*L*L;t++)Sweep(S,m_sum,prob,t,L);
    for(t=0;t<T_MEAS*L*L;t++){
        m_sum[t] = 0;
        Sweep(S,m_sum,prob,t,L);
        m2_avg[l] +=  m_sum[t]*m_sum[t];
        m_avg[l] += m_sum[t];}
}
