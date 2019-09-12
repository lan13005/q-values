#ifndef RUN_H
#define RUN_H
const int dim=1;
const int kDim=10;

double calc_distance( double phaseSpace_1[dim], double phaseSpace_2[dim] ){
	double sum = 0;
	double diff=0;
	for (int i=0; i<dim; ++i){
		diff = phaseSpace_1[i]-phaseSpace_2[i];
		sum += diff*diff;
	}
	return sum;
}


int numDOFbkg = 1;
Double_t background(Double_t *x, Double_t *par){
	return par[0];
	//return par[0]+par[1]*x[0];
	//return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
	//return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
	//return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
}

int numDOFsig = 3;
Double_t signal(Double_t *x, Double_t *par){
	return par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));
	//return (x[0]-par[0])*(x[0]-par[0]);

}

Double_t fitFunc(Double_t *x, Double_t *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}


#endif
