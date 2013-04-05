package com.jove.ldpc;
public class DensityEvolution {

//constant used in Phi()
static double ALFA=-0.4527;
static double BETA= 0.0218;
static double GAMA= 0.86;
static double PI= Math.PI;

//iteration parameter
static int ITERATION= 1000;
static double SIGMA_BEGIN= 0.2;
static double SIGMA_PRECISION= 0.0001;

//check irregular: 1;
//check regular: 0;
boolean chk_irre;

// parameter used in phi_inv
// phi_inv_sample
static int PHI_SAMPLE_LENGTH = 100000;
static double PHI_SAMPLE_RANGE = 40;
static double PHI_SAMPLE_STEP = PHI_SAMPLE_RANGE/PHI_SAMPLE_LENGTH;
double[] phi_inv_in;
double[] phi_inv_out;
//phi_inv function in and out vector
//---------------------------------
//regular LDPC degree
//---------------------------------
int dv,dc;
double code_rate;
//----------------------------------
//irregular LDPC degree distribution
//----------------------------------
int vdd_cnt=0;
int cdd_cnt=0;

//edge perspective
double[][] vdd_e;
double[][] cdd_e;
//node perspective
double[][] vdd_n;
double[][] cdd_n;

double sigma;

double[][] allocDouble(int i,int j)
{
	double[][] a  = new double[i][];
	for( int k=0;k<i;k++ )
	{
		a[k]=new double[j];
		for( int l=0;l<j;l++ )
			a[k][l] = 0.0;
	}
	return a;
}

void init()
{
	int i;
	phi_inv_in = new double[PHI_SAMPLE_LENGTH];
	for( i=0;i<PHI_SAMPLE_LENGTH;i++ )
		phi_inv_in[i]=0.;
	phi_inv_out = new double[PHI_SAMPLE_LENGTH];
	for( i=0;i<PHI_SAMPLE_LENGTH;i++ )
		phi_inv_out[i]=0.;
	vdd_e = allocDouble(20,2);
	cdd_e = allocDouble(20,2);
	vdd_n = allocDouble(20,2);
	cdd_n = allocDouble(20,2);
    for (i=0;(PHI_SAMPLE_RANGE-i*PHI_SAMPLE_STEP)>0;i++) 
    { 
    	phi_inv_out[i]=PHI_SAMPLE_RANGE-(double)i*PHI_SAMPLE_STEP;
		phi_inv_in[i]=Phi(phi_inv_out[i]);
    }		

//    for( i=0;i<phi_inv_in.length;i++ )
//    {
//    	if( (i%1000)==0 )
//    		System.out.println((double)i*PHI_SAMPLE_STEP+":f("+phi_inv_out[i]+")="+phi_inv_in[i]);
//    }
	
}
//-------------------------
//function declearation
//-------------------------
//void usage();
//int check_threshold(double);
//double re_mu_update (double,double);
//double irre_mu_update (double,double);
//double Phi(double);
//double Phi_inv(double);
//int compare (const void *, const void *);
//some transfer function
//void vdd_n2vdd_e();
//void cdd_n2cdd_e();
//void vdd_e2vdd_n();
//void cdd_e2cdd_n();
//void vdd_n2cdd_n();
//void vdd_e2cdd_e();
//void print_degree_dis();
//void check_degree_dis();

void usage()
{
	System.out.println("usage: GADE_ALL_NEW");
    System.out.println("the file named argv.dat must be contain in the same fold of the program");
    System.out.println("regular code example:");
    System.out.println("---------------");
    System.out.println("-r");
    System.out.println("3");
    System.out.println("6");
    System.out.println("---------------");
    System.out.println("irregular code example:");
    System.out.println("---------------");
    System.out.println("-i-e or -i-n");
    System.out.println("3(vdd counter)");
    System.out.println("2 0.3");
    System.out.println("3 0.4");
    System.out.println("4 0.3");
    System.out.println("3(cdd counter)");
    System.out.println("6 0.8");
    System.out.println("7 0.2");
    System.out.println("---------------");
	System.out.println("-i-e-r -i-n-r");
	System.out.println("0.5(code rate)");
    System.out.println("3(vdd counter)");
    System.out.println("2 0.3");
    System.out.println("3 0.4");
    System.out.println("4 0.3");
    System.out.println("---------------");
}

//----------------
//this main is DE
//----------------

//check if the threshold is satisfied
int check_threshold(double thre) 
{
	assert(thre>=0);
	double mu0 = 2/(thre*thre);
	double mu;
	double mu_next=0;
	for (int j=0;j<ITERATION;j++) 
	{
		mu = mu_next;
		mu_next = (chk_irre)?
			re_mu_update(mu0,mu)//check regular LDPC
			:
			irre_mu_update(mu0,mu);//check irregular LDPC
//		System.out.println("mu:"+mu_next);
	if (mu_next >= (PHI_SAMPLE_RANGE-1)) 
		{
			return 1;
		}
	}
	return 0;
}

//regular mu update: with constant dv and dc
double re_mu_update (double mu0,double mu) {
	assert(mu0>=0);
	assert(mu>=0);
	double temp1=0;
	double temp2=0;
	temp1 = 1-Phi(mu0+(dv-1)*mu);
	temp2 = 1-Math.pow(temp1,(dc-1));
	return Phi_inv(temp2);
}

//irregular mu update: with different dv and different dc
//the output is wrong
double irre_mu_update (double mu0,double mu){
	assert(mu0>=0);
	assert(mu>=0);
	
	double temp1=0;
	double temp2=0;
	
	for (int v=0;v<vdd_cnt;v++) {
		temp1 += vdd_e[v][1]*Phi(mu0+(vdd_e[v][0]-1)*mu);
	}
	temp1 = 1 - temp1;
	
	for (int c=0;c<cdd_cnt;c++) {
		temp2 += cdd_e[c][1]*Phi_inv(1-Math.pow(temp1,(cdd_e[c][0]-1)));
	}
	return temp2;
}

//Phi function
double Phi(double x) 
{
	//assert x>0
	assert(x>=0);

    if(x<=10) 
    {
		return Math.exp(ALFA*Math.pow(x,GAMA)+BETA);
	}
	else 
	{
		return Math.sqrt(PI/x)*Math.exp(-x/4)*(1+1/(14*x)-3/(2*x));
	}
}

//inverse Phi function v2: using bsearch
double Phi_inv(double y) 
{
    int index = binarySearch ( phi_inv_in, 0, PHI_SAMPLE_LENGTH-1, y );
    if( index == -1 && index==PHI_SAMPLE_LENGTH-1 )
    	return  (index!=-1) ? (phi_inv_out[index]): PHI_SAMPLE_RANGE;
    else
    {
    	if( index==0 && y<phi_inv_in[0] )
    		return phi_inv_out[index];
    	double a=phi_inv_in[index+1]-phi_inv_in[index];
    	double b=phi_inv_out[index+1]-phi_inv_out[index];
    	double c=b/a*(y-phi_inv_in[index])+phi_inv_out[index];
    	return c;
    }
}

public static int binarySearch( double[] v, int from, int to, double y )
{  
	if( to==from || to==from+1 )
		return from;
	if (from > to)
		return -1;
	int mid = (from + to) / 2;
	if (y <v[mid] ) 
		return binarySearch(v, from, mid, y);
	else
		return binarySearch(v, mid, to, y);
}


//transfer node perspective degree distribution to edge perspective
//vdd_n->vdd-e
void vdd_n2vdd_e()
{
	int i;
    double vdd_n_sum=0;
    for(i=0;i<vdd_cnt;i++) 
    {
        vdd_n_sum += vdd_n[i][0]*vdd_n[i][1];
    }
    for(i=0;i<vdd_cnt;i++) 
    {
        vdd_e[i][0]=vdd_n[i][0];
        vdd_e[i][1] = (vdd_n[i][0]*vdd_n[i][1])/vdd_n_sum;
    } 
}
//cdd_n->cdd_e	
void cdd_n2cdd_e ()
{
    int i;
    double cdd_n_sum=0;
    for(i=0;i<cdd_cnt;i++) {
        cdd_n_sum += cdd_n[i][0]*cdd_n[i][1];
    }
    for(i=0;i<vdd_cnt;i++) {
        cdd_e[i][0]=cdd_n[i][0];
        cdd_e[i][1] = (cdd_n[i][0]*cdd_n[i][1])/cdd_n_sum;
    }     
}

//transfer edge perspective degree distribution to node perspective
//vdd_e->vdd_n
void vdd_e2vdd_n() 
{
	int i;
	double vdd_e_sum=0;
    for(i=0;i<vdd_cnt;i++) {
        vdd_e_sum += (vdd_e[i][1]/vdd_e[i][0]);
    }
    for(i=0;i<vdd_cnt;i++) {
        vdd_n[i][0]=vdd_e[i][0];
        vdd_n[i][1] = (vdd_e[i][1]/vdd_e[i][0])/vdd_e_sum;
    } 
}

//cdd_e->cdd_n
void cdd_e2cdd_n ()
{
    int i;
    double cdd_e_sum=0;
    for(i=0;i<cdd_cnt;i++) {
        cdd_e_sum += cdd_e[i][1]/cdd_e[i][0];
    }
    for(i=0;i<vdd_cnt;i++) {
        cdd_n[i][0]=cdd_e[i][0];
        cdd_n[i][1] = (cdd_e[i][1]/cdd_e[i][0])/cdd_e_sum;
    }     
}

//vdd_n -> cdd_n
void vdd_n2cdd_n() 
{
	double vdd_sum=0;
	double cdd_sum=0;
	int i;
	for(i=0;i<vdd_cnt;i++) {
		vdd_sum += vdd_n[i][0]*vdd_n[i][1];
	}
	cdd_sum = vdd_sum/(1-code_rate);
	cdd_n[0][0]=(int) cdd_sum;
	cdd_n[1][0]=cdd_n[0][0]+1;
	cdd_n[1][1]=cdd_sum-(double)cdd_n[0][0];
	cdd_n[0][1]=1-cdd_n[1][1];
}

//vdd_e->vdd_n->cdd_n->cdd_e
void vdd_e2cdd_e() 
{
	
	vdd_e2vdd_n();
	vdd_n2cdd_n();
	cdd_n2cdd_e();
	
}

void print_degree_dis() {
	int i;
	double vdd_sum=0;
	double cdd_sum=0;
	System.out.println("degree distribution(edge perspective):");
	System.out.println("Varible Node:");
	for(i=0;i<vdd_cnt;i++) {
		System.out.println(vdd_e[i][0]+":\t"+vdd_e[i][1]);
	}
	System.out.println("Check Node:");
	for(i=0;i<cdd_cnt;i++) 
	{
		System.out.println(cdd_e[i][0]+":\t"+cdd_e[i][1]);
	}
	System.out.println("degree distribution(node perspective):");
	System.out.println("Varible Node:");
	for(i=0;i<vdd_cnt;i++) {
		System.out.println(vdd_n[i][0]+":\t"+vdd_n[i][1]);
	}
	System.out.println("Check Node:");
	for(i=0;i<cdd_cnt;i++) {
		System.out.println(cdd_n[i][0]+":\t"+cdd_n[i][1]);
	}
	
	for(i=0;i<vdd_cnt;i++) {
		vdd_sum += vdd_n[i][0]*vdd_n[i][1];
	}
	System.out.println("Average Column Weight:"+vdd_sum);
	for(i=0;i<cdd_cnt;i++) {
		cdd_sum += cdd_n[i][0]*cdd_n[i][1];
	}
	System.out.println("Average Row Weight:"+cdd_sum);
	System.out.println(((chk_irre)?"regular":"irregular")+" code threshold " + sigma );
	double EbN0 = 10.*Math.log10(0.5/(sigma*sigma*code_rate));
	System.out.println(" EbN0: " + EbN0 );
	
}

void check_degree_dis() {
	int i;
	double vdd_n_sum=0;
	double cdd_n_sum=0;
	double vdd_e_sum=0;
	double cdd_e_sum=0;
	
	for (i=0;i<vdd_cnt;i++) 
	{
		vdd_n_sum += vdd_n[i][1];
	}
	for (i=0;i<cdd_cnt;i++) {
		cdd_n_sum += cdd_n[i][1];
	}	
	for (i=0;i<vdd_cnt;i++) {
		vdd_e_sum += vdd_e[i][1];
	}
	for (i=0;i<cdd_cnt;i++) {
		cdd_e_sum += cdd_e[i][1];
	}
	if (vdd_n_sum<0.99 || vdd_n_sum>1.01 || cdd_n_sum<0.99 || cdd_n_sum>1.01
	  ||vdd_e_sum<0.99 || vdd_e_sum>1.01 || cdd_e_sum<0.99 || cdd_e_sum>1.01)
		System.out.println("Warning: The sum of Degree Distribution not equal to 1!");
}

DensityEvolution(double rate, Degree deg)
{
	chk_irre = false;
    code_rate = rate;
	init();
    int i;
    vdd_cnt=deg.tab.length;
    
	for(i=0;i<vdd_cnt;i++) {
        vdd_n[i][0] = deg.tab[i].deg;
        vdd_n[i][1] = deg.tab[i].degFrac;
    }
    cdd_cnt=2;
	vdd_n2cdd_n();
    vdd_n2vdd_e();
	cdd_n2cdd_e();
	
	for(sigma=SIGMA_BEGIN;;sigma=sigma+SIGMA_PRECISION) 
	{
		if (check_threshold(sigma)==0)
			break;
//		System.out.println(sigma);
	}
}
public static void main(String[] args) 
{
//	double df[] = {0.521814,    0.271293,     0.0,    0.206893};
//	int dg[] = {   2,            3,            4,          5};
//	double df[] = { 0.4717,    0.3336,    0.0102,    0.0426,    0.0070,    0.0049,    0.1300 };
//	int dg[] = {2,3,4,5,8,10,12};
	int dg[] = {2,	3,	12};
	double df[] = { 0.750000, 	0.166667, 	0.083333333};


	Degree deg = new Degree( dg,df );
	
	DensityEvolution ttu = new DensityEvolution(0.125,deg);
	
	ttu.print_degree_dis();
	ttu.check_degree_dis();
//	System.out.println( sigma-SIGMA_PRECISION );
	return;
}

}
