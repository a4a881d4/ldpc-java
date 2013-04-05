package com.jove.ldpc;

public class PegMask {

	/**
	 * @param args
	 */
	int M;
	int N;
	int[] degSeq;
	int[][] mask;
	int[] rowWeight;
	int Weight;
	Degree deg;
	PegMask( Degree deg )
	{
		this.deg = new Degree( deg );
	}
	void init(int M, int N)
	{
		this.M = M;
		this.N = N;
		degSeq = deg.degSeq(N);
		int i,j;
		mask = new int[M][];
		rowWeight = new int[M];
		for( i=0;i<M;i++ )
		{
			mask[i] = new int[N];
			for( j=0;j<N;j++ )
				mask[i][j]=0;
			rowWeight[i]=0;
		}
		Weight=0;
	}
	int matchTwo(int i, int j)
	{
		int ret = 0;
		for( int l = 0;l<M;l++ )
			ret += mask[l][i]*mask[l][j];
		return ret;
	}
	int random(int begin, int end, int k, int[] ret)
	{
		int size = end - begin;
		int i,j;
		for( i=1;i<k;i++ )
		{
			ret[i] = ret[0]; 
			j=0;
			while( j<i )
			{
				if( ret[i]==ret[j] )
				{
					ret[i]=getAPos(size,begin);
					j=0;
				}
				else
					j++;
			}
		}
		return k;
	}
	void show()
	{
		int i,j;
		for(i=0;i<M;i++)
		{
			for(j=0;j<N;j++)
				System.out.printf("%1d",mask[i][j]);
			System.out.printf("\n");
		}
	}
	int getAPos(int size,int begin)
	{
		int ret = (int)Math.round(Math.random()*(double)size)+begin;
		int avg =(int)Math.round((double)Weight/(double)M+0.5);
		int retry = 0;
		while( rowWeight[ret]>=avg )
		{
			ret = (int)Math.round(Math.random()*(double)size)+begin;
			retry++;
			if( retry>10 )
				break;
		}
		return ret;
	}
	int fillMask()
	{
		int i,j;
		int retry;
		int level = 2;
		mask[M-1][N-1] = 1;
		mask[M-2][N-1] = 1;
		int back = M-2;
		int[] pos=new int[M];
		
		for( i=N-2;i>=0;i-- )
		{
			retry = 0;
			
			while(true)
			{
				int num = degSeq[N-1-i];
				if( back !=0 )
				{
					pos[0]=back-1;
				}
				else
				{
					pos[0]=rlight(i+1);
				}
				int ret = random(back,M-1,num,pos);
				for(j=0;j<num;j++)
				{
					mask[pos[j]][i]=1;
				}
				for(j=i+1;j<N;j++)
					if(matchTwo(i,j)>=level)
						break;
				if(j==N)
				{
					for(j=0;j<num;j++)
					{
						rowWeight[pos[j]]++;
					}
					Weight+=num;
					if( back!=0 )
						back--;
					break;
				}
				else
				{
					retry ++;
					if( retry>100 )
					{
						level++;
						retry = 0;
					}
					for(j=0;j<num;j++)
						mask[pos[j]][i]=0;
				}
			}
		}
		return level;
	}
	int[] rweight(int print)
	{
		int[] ret = new int[M];
		int i,j;
		for(i=0;i<M;i++)
		{
			ret[i]=0;
			for(j=0;j<N;j++)
				ret[i]+=mask[i][j];
			if( print==1 )
				System.out.printf("%d : %d\n",i,ret[i]);
		}
		return ret;
	}
	int rlight(int b)
	{
		int ret=-1;
		int i,j;
		int rec=N;
		for(i=0;i<M;i++)
		{
			int cnt=0;
			for(j=b;j<N;j++)
				cnt+=mask[i][j];
			if( cnt<rec )
			{
				rec = cnt;
				ret=i;
			}
		}
		return ret;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double df[] = {0.521814,    0.271293,     0.0,    0.206893};
		int dg[] = {   2,            3,            4,          7};
		Degree deg = new Degree(dg,df);
		PegMask aM = new PegMask(deg);
		aM.init(34, 68);
		int l = aM.fillMask();
		aM.show();
		System.out.printf("Level is %d\n",l);
		aM.rweight(1);
	}

}
