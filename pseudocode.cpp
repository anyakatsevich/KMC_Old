t = 0;
Tfinal = ... 
int L = ...
int *h = .... //initialize
double *rates = ...//compute from h
R = sum of rates = R_max; // (R = R_max when there's only one block)
t_stop = C/R_max;
while (t < Tfinal){
    double t_1 = 0;
    double t_2 = 0;
    while (t_1 < t_stop){
        dt = draw from Exp(R)
        
        i = draw from distribution {j w.p. r_j/R, j=0,...,L-1}
        if (i < L/2){
            jump left or right w.p. 1/2
            update values of h around i
            update rates around i
            update R
        }
        t_1 += dt;
    }
    //COMMUNICATE BETWEEN PROCESSORS
    while (t_2 < t_stop){
        dt = draw from Exp(R)
        
        i = draw from distribution {j w.p. r_j/R, j=0,...,L-1}
        if (i > L/2){
            jump left or right w.p. 1/2
            update values of h around i
            update rates around i
            update R
        }
        t_2 += dt;
    }
    t += t_stop;
    //COMMUNICATE BETWEEN PROCESSORS 
    //RECOMPUTE R_max AND t_stop


}