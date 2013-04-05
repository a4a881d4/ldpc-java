package com.jove.ldpc;

public class Random {
	long seed;  
	long seed_u;
	Random() 
	{
		seed=987654321;
		seed_u=123456789l;
	}
	void bubbleSort(int a[], int size)
	{
		{
			int i, temp;
			for(int pass=1; pass<size; pass++) 
			{
				for(i=0;i<size-pass;i++)
					if(a[i]>a[i+1])
					{
						temp=a[i];
						a[i]=a[i+1];
						a[i+1]=temp;
					}
			}
		}
	}
	double gauss(double sdev, double mean)
	{
		double sum=0.0;
		for (int i=1;i<=12;i++)
		{ 
			sum=sum+Math.random(); 
		}
		return (sum-6.)*sdev+mean;
	}
	double uniform(double a, double b)
	{
		double t = Math.random();
		t=a+(b-a)*t;
		return t;
	}
	int uniform(int  a, int b) // [a, b)
	{
		double t = Math.random();
		int tt;
		if(b==a+1) return(a);
		t=a+(b-a)*t;
		tt=(int)t;
		if(tt<a) tt=a;
		else if(tt>=b) tt=b-1;
		return(tt);
	}
}
