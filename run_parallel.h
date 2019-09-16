#ifndef RUN_PARALLEL_H
#define RUN_PARALLEL_H

const int dim=5;

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

int numDOFsig = 6;
Double_t signal(Double_t *x, Double_t *par){
	return par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])) + par[3]*exp(-0.5*((x[0]-par[4])/par[5])*((x[0]-par[4])/par[5]));
	//return (x[0]-par[0])*(x[0]-par[0]);

}

Double_t fitFunc(Double_t *x, Double_t *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}

void standardizeArray(double inputVector[], int nentries, string name){
	double max_inputVector = *max_element(inputVector, inputVector+nentries);
	double min_inputVector = *min_element(inputVector, inputVector+nentries);
	//double extrema_inputVector[2] = {max_inputVector, min_inputVector};
	for (int ientry=0; ientry<nentries; ++ientry){
		inputVector[ientry] = (inputVector[ientry]-min_inputVector)/(max_inputVector-min_inputVector);
	}
	cout << "Max,min "+name+": " << max_inputVector << "," << min_inputVector << endl;
	cout << "	Finished standardizing "+name << endl;
}

#endif
