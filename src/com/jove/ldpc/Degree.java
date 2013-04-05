package com.jove.ldpc;
import java.io.PrintWriter;

public class Degree {

	/**
	 * @param args
	 */
	class DegItem {
		int deg;
		double degFrac;
	}
	DegItem[] tab;
	
	Degree( Degree deg)
	{
		tab = new DegItem[deg.tab.length];
		for( int i=0;i<deg.tab.length;i++ )
		{
			tab[i] = new DegItem();
			tab[i].degFrac=deg.tab[i].degFrac;
			tab[i].deg=deg.tab[i].deg;
		}
	}
	public Degree( int[] deg, double[] degFrac)
	{
		assert(deg.length==degFrac.length);
		tab = new DegItem[deg.length];
		for( int i=0;i<tab.length;i++ )
		{
			tab[i] = new DegItem();
			tab[i].degFrac=degFrac[i];
			tab[i].deg=deg[i];
		}
	}
	public void log(PrintWriter out) {
		int i;		
		out.printf("# ");
		for( i=0;i<tab.length;i++ )
			out.printf("  %8d",tab[i].deg);
		out.printf("\n");
		out.printf("# ");
 		for( i=0;i<tab.length;i++ )
			out.printf("  %8f",tab[i].degFrac);
		out.printf("\n");
	}
		
	boolean check()
	{
		double sum = 0;
		for( int i=0;i<tab.length;i++ )
			sum += tab[i].degFrac;
		if( Math.abs(sum-1)>1e-6 )
			return false;
		return true;
	}
	int[] degSeq(int N)
	{
		int[] ret = new int[N];
		int i,j;
		if( !check() )
		{
			System.out.printf("Degree distrubtion Error( Sum != 1.0 )\n");
			return null;
		}
		double[] df = new double[tab.length];
		df[0] = tab[0].degFrac;
		for( i=1;i<tab.length;i++)
			df[i] = tab[i].degFrac + df[i-1];
		for(i=0;i<N;i++) 
		{
			double dtmp=(double)i/N;
		    for(j=tab.length-1;j>=0;j--) 
		    {
		      if(dtmp>df[j]) break;
		    }
		    if(dtmp<df[0]) 
		    	ret[i]=tab[0].deg;
		    else 
		    	ret[i]=tab[j+1].deg;
		}
		int[] res = new int[N];
		int[] primer={11,7,5,3,1};
		j=0;
		while(primer[j]!=1)
		{
			if((N%primer[j])!=0)
				break;
			j++;
		}
		for( i=0;i<N;i++ )
		{
			res[i]=ret[(i*primer[j])%N];
		}
		return ret;
	}
	
	Degree node2edge()
	{
		int i;
	    double vdd_n_sum=0;
	    int[] dg = new int[tab.length];
	    double[] df = new double[tab.length];
	    for(i=0;i<tab.length;i++) 
	    {
	        vdd_n_sum += tab[i].deg*tab[i].degFrac;
	    }
	    for(i=0;i<tab.length;i++) 
	    {
	    	dg[i] = tab[i].deg;
	        df[i] = (tab[i].deg*tab[i].degFrac)/vdd_n_sum;
	    }
	    Degree ret = new Degree( dg, df ); 
	    return ret;
	}

	Degree edge2node()
	{
		int i;
	    double vdd_n_sum=0;
	    int[] dg = new int[tab.length];
	    double[] df = new double[tab.length];
	    for(i=0;i<tab.length;i++) 
	    {
	        vdd_n_sum += tab[i].deg/tab[i].degFrac;
	    }
	    for(i=0;i<tab.length;i++) 
	    {
	    	dg[i] = tab[i].deg;
	        df[i] = (tab[i].deg/tab[i].degFrac)/vdd_n_sum;
	    }
	    Degree ret = new Degree( dg, df ); 
	    return ret;
	}
	Degree v2c( double rate ) 
	{
		double vdd_sum=0;
		double cdd_sum=0;
		int i;
		int[] dg = new int[2];
		double[] df = new double[2];
		for(i=0;i<tab.length;i++) {
			vdd_sum += tab[i].deg*tab[i].degFrac;
		}
		
		cdd_sum = vdd_sum/(1.-rate);
		dg[0]=(int)Math.floor(cdd_sum);
		dg[1]=dg[0]+1;
		df[1]=cdd_sum-(double)dg[0];
		df[0]=1.-df[1];
		return new Degree(dg,df);
	}
	static double rateN( Degree V, Degree C)
	{
		Degree nV = V.node2edge();
		Degree nC = C.node2edge();
		return Degree.rateE(nV,nC);
	}
	static double rateE( Degree V, Degree C)
	{
		double sumV=0.,sumC=0.;
		int i;
		for(i=0;i<V.tab.length;i++)
			sumV += V.tab[i].degFrac/V.tab[i].deg;
		for(i=0;i<C.tab.length;i++)
			sumC += C.tab[i].degFrac/C.tab[i].deg;
		return 1.-sumC/sumV;
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
