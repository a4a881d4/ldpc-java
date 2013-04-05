package com.jove.ldpc;

public class BigGirth {

	  int M, N;
	  int K;
	  int EXPAND_DEPTH;
	  int[][] H;

	  int[] localGirth;
	  
	  NodesInGraph[] nodesInGraph;
	  Random myrandom;

	  BigGirth(int M, int N, int[] symbolDegSequence, int sglConcent, int tgtGirth)
	  {
		  int i, j, k, m, index;
		  int[] mid;
		  int[] localDepth= new int[1];
		  localDepth[0]=100;
		  
		  EXPAND_DEPTH=(tgtGirth-4)/2; 
		  if(EXPAND_DEPTH<0) EXPAND_DEPTH=0;

		  myrandom=new Random(); 

		  this.M=M;
		  this.N=N;
		  H=null;
		  mid=new int[M];

		  localGirth=new int[N];

		  nodesInGraph=new NodesInGraph [N];
		  for(i=0;i<N;i++)
		  {	
			  nodesInGraph[i] = new NodesInGraph();
			  nodesInGraph[i].setNumOfConnectionSymbolBit(symbolDegSequence[i]);
		  }
		  j=0;
		  for(k=0;k<N;k++) 
			  j+=symbolDegSequence[k];
		  k=j/M;
		  for(i=0;i<M;i++) 
			  mid[i]=k;
		  for(i=0;i<j-k*M;i++) 
			  mid[i]++;
		  k=0; 
		  for(i=0;i<M;i++) 
			  k+=mid[i];
		  if(k!=j) 
		  {
			  System.out.println("Wrong in computing maxDegParity!");
		  }

		  for(i=0;i<M;i++) 
		  {
			  if(sglConcent==0) 
				  nodesInGraph[i].initConnectionParityBit(mid[i]);
			  else  
				  nodesInGraph[i].initConnectionParityBit(); 
		  } 

		  for(k=0;k<N;k++)
		  {
			  m=1000000;
			  index=-1;
			  for(i=0;i<M;i++)
			  {
				  if( nodesInGraph[i].numOfConnectionParityBit<m && 
						  nodesInGraph[i].numOfConnectionParityBit<nodesInGraph[i].maxDegParity
						  ) 
				  {
					  m=nodesInGraph[i].numOfConnectionParityBit;
					  index=i;
				  }
			  }
			  nodesInGraph[k].connectionSymbolBit[0]=index;//least connections of parity bit

			  int iter=0; 
			  ITER:
				  while(true)
				  {
					  localGirth[k]=100;
					  for(m=1;m<nodesInGraph[k].numOfConnectionSymbolBit;m++)
					  {
						  nodesInGraph[k].connectionSymbolBit[m]=selectParityConnect(k, m, localDepth);
						  localGirth[k]=(localGirth[k]>localDepth[0])?localDepth[0]:localGirth[k];      
						  if(k>0 && localGirth[k]<localGirth[k-1] && iter<20) 
						  {
							  iter++; 
							  continue ITER;
						  }
					  }
					  if(localGirth[k]==0 && iter<30) 
					  {
						  iter++; 
						  continue ITER;
					  }
					  break ITER;
				  }
		  
			  System.out.print("k="+k+"  ");
			  for(m=0;m<nodesInGraph[k].numOfConnectionSymbolBit;m++)
				  System.out.print(nodesInGraph[k].connectionSymbolBit[m]+" ");
			  System.out.print("LocalGirth="+(2*localGirth[k]+4));
			  System.out.println();
			  updateConnection(k);
		  }
		  loadH();
	  }
	  BigGirth(){}

	  
	  private int selectParityConnect(int kthSymbol, int mthConnection, int[] cycle)
	  {
		  int i, j, k, index, mincycles, numCur,  cpNumCur;
		  int[] tmp, med;
		  int[] current;//take note of the covering parity bits

		  mincycles=0;
		  tmp=new int[M]; 
		  med=new int[M];

		  numCur=mthConnection;
		  current=new int[mthConnection];
		  for(i=0;i<mthConnection;i++) 
			  current[i]=nodesInGraph[kthSymbol].connectionSymbolBit[i];

		  LOOP:
			  while(true)
			  {
				  mincycles++;
				  for(i=0;i<M;i++) tmp[i]=0;
				  //maintain 
				  for(i=0;i<mthConnection;i++) 
					  tmp[nodesInGraph[kthSymbol].connectionSymbolBit[i]]=1;
				  for(i=0;i<numCur;i++)
				  {
					  for(j=0;j<nodesInGraph[current[i]].numOfConnectionParityBit;j++)
					  {
						  for(k=0;k<nodesInGraph[nodesInGraph[current[i]].connectionParityBit[j]].numOfConnectionSymbolBit;k++)
						  {
							  tmp[nodesInGraph[nodesInGraph[current[i]].connectionParityBit[j]].connectionSymbolBit[k]]=1;
						  }
					  }
				  }

				  index=0; 
				  cpNumCur=0;
				  for(i=0;i<M;i++) 
				  {
					  if(tmp[i]==1) 
						  cpNumCur++;
					  if(tmp[i]==1 || nodesInGraph[i].numOfConnectionParityBit>=nodesInGraph[i].maxDegParity) 
						  index++;   
				  }
				  if(cpNumCur==numCur) //can not expand any more
				  {  //additional handlement to select one having least connections
					  j=10000000; //dummy number
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0 && nodesInGraph[i].numOfConnectionParityBit<j && nodesInGraph[i].numOfConnectionParityBit<nodesInGraph[i].maxDegParity)
							  j=nodesInGraph[i].numOfConnectionParityBit;
					  }
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0)
						  {
							  if(nodesInGraph[i].numOfConnectionParityBit!=j || nodesInGraph[i].numOfConnectionParityBit>=nodesInGraph[i].maxDegParity)
							  {
								  tmp[i]=1;
							  }
						  }
					  }
					  index=0;
					  for(i=0;i<M;i++) 
						  if(tmp[i]==1) 
							  index++;
					  //----------------------------------------------------------------
					  j=myrandom.uniform(0, M-index)+1; //randomly selected
					  index=0;
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0) 
							  index++;
						  if(index==j) 
							  break;
					  }
					  tmp=null;
					  current=null;
					  return(i);
				  }
				  else if(index==M || mincycles>EXPAND_DEPTH)
				  {//covering all parity nodes or meet the upper bound on cycles

					  cycle[0]=mincycles-1;
					  for(i=0;i<M;i++) 
						  tmp[i]=0;
					  for(i=0;i<numCur;i++) 
						  tmp[current[i]]=1;
					  index=0;
					  for(i=0;i<M;i++) 
						  if(tmp[i]==1) 
							  index++;
					  if(index!=numCur) 
					  {
						  System.out.println("Error in the case of (index==M)");
					  }
					  //additional handlement to select one having least connections
					  j=10000000; 
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0 && nodesInGraph[i].numOfConnectionParityBit<j && nodesInGraph[i].numOfConnectionParityBit<nodesInGraph[i].maxDegParity)
							  j=nodesInGraph[i].numOfConnectionParityBit;
					  }
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0)
						  {
							  if(nodesInGraph[i].numOfConnectionParityBit!=j || nodesInGraph[i].numOfConnectionParityBit>=nodesInGraph[i].maxDegParity)
							  {
								  tmp[i]=1;
							  }
						  }
					  }

					  index=0;
					  for(i=0;i<M;i++) 
						  if(tmp[i]==1) 
							  index++;

					  j=myrandom.uniform(0, M-index)+1;
					  index=0;
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0) 
							  index++;
						  if(index==j) 
							  break;
					  }
					  tmp=null;
					  current=null;
					  return(i);
				  }
				  else if(cpNumCur>numCur && index!=M)
				  {
					  current=null;
					  numCur=cpNumCur;
					  current=new int[numCur];
					  index=0;
					  for(i=0;i<M;i++) 
					  {
						  if(tmp[i]==1) 
						  {
							  current[index]=i; 
							  index++;
						  }
					  }
					  continue LOOP;
				  }
				  else 
				  {
					  System.out.println("Should not come to this point...");
					  System.out.println("Error in BigGirth::selectParityConnect()");
					  return(-1);
				  }
			  }
	  }

	  private int selectParityConnectTr(int kthSymbol, int mthConnection, int[] cycle)
	  {
		  int i, j, k, index, mincycles, numCur,  cpNumCur;
		  int[] tmp, med;
		  int[] current;//take note of the covering parity bits

		  mincycles=0;
		  tmp=new int[M]; 
		  med=new int[M];

		  numCur=mthConnection;
		  current=new int[mthConnection];
		  for(i=0;i<mthConnection;i++) 
			  current[i]=nodesInGraph[kthSymbol].connectionSymbolBit[i];

		  LOOP:
			  while(true)
			  {
				  mincycles++;
				  for(i=0;i<M;i++) tmp[i]=0;
				  //maintain 
				  for(i=0;i<mthConnection;i++) 
					  tmp[nodesInGraph[kthSymbol].connectionSymbolBit[i]]=1;
				  for(i=0;i<numCur;i++)
				  {
					  for(j=0;j<nodesInGraph[current[i]].numOfConnectionParityBit;j++)
					  {
						  for(k=0;k<nodesInGraph[nodesInGraph[current[i]].connectionParityBit[j]].numOfConnectionSymbolBit;k++)
						  {
							  tmp[nodesInGraph[nodesInGraph[current[i]].connectionParityBit[j]].connectionSymbolBit[k]]=1;
						  }
					  }
				  }
				  if( kthSymbol<M-1 )
				  {
					  for( i=kthSymbol+2;i<M;i++ )
						  tmp[i]=1;
				  }
				  index=0; 
				  cpNumCur=0;
				  for(i=0;i<M;i++) 
				  {
					  if(tmp[i]==1) 
						  cpNumCur++;
					  if(tmp[i]==1 || nodesInGraph[i].numOfConnectionParityBit>=nodesInGraph[i].maxDegParity) 
						  index++;   
				  }
				  if(cpNumCur==numCur) //can not expand any more
				  {  //additional handlement to select one having least connections
					  j=10000000; //dummy number
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0 && nodesInGraph[i].numOfConnectionParityBit<j && nodesInGraph[i].numOfConnectionParityBit<nodesInGraph[i].maxDegParity)
							  j=nodesInGraph[i].numOfConnectionParityBit;
					  }
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0)
						  {
							  if(nodesInGraph[i].numOfConnectionParityBit!=j || nodesInGraph[i].numOfConnectionParityBit>=nodesInGraph[i].maxDegParity)
							  {
								  tmp[i]=1;
							  }
						  }
					  }
					  index=0;
					  for(i=0;i<M;i++) 
						  if(tmp[i]==1) 
							  index++;
					  //----------------------------------------------------------------
					  j=myrandom.uniform(0, M-index)+1; //randomly selected
					  index=0;
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0) 
							  index++;
						  if(index==j) 
							  break;
					  }
					  tmp=null;
					  current=null;
					  return(i);
				  }
				  else if(index==M || mincycles>EXPAND_DEPTH)
				  {//covering all parity nodes or meet the upper bound on cycles

					  cycle[0]=mincycles-1;
					  for(i=0;i<M;i++) 
						  tmp[i]=0;
					  for(i=0;i<numCur;i++) 
						  tmp[current[i]]=1;
					  index=0;
					  for(i=0;i<M;i++) 
						  if(tmp[i]==1) 
							  index++;
					  if(index!=numCur) 
					  {
						  System.out.println("Error in the case of (index==M)");
					  }
					  //additional handlement to select one having least connections
					  j=10000000; 
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0 && nodesInGraph[i].numOfConnectionParityBit<j && nodesInGraph[i].numOfConnectionParityBit<nodesInGraph[i].maxDegParity)
							  j=nodesInGraph[i].numOfConnectionParityBit;
					  }
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0)
						  {
							  if(nodesInGraph[i].numOfConnectionParityBit!=j || nodesInGraph[i].numOfConnectionParityBit>=nodesInGraph[i].maxDegParity)
							  {
								  tmp[i]=1;
							  }
						  }
					  }

					  index=0;
					  for(i=0;i<M;i++) 
						  if(tmp[i]==1) 
							  index++;

					  j=myrandom.uniform(0, M-index)+1;
					  index=0;
					  for(i=0;i<M;i++)
					  {
						  if(tmp[i]==0) 
							  index++;
						  if(index==j) 
							  break;
					  }
					  tmp=null;
					  current=null;
					  return(i);
				  }
				  else if(cpNumCur>numCur && index!=M)
				  {
					  current=null;
					  numCur=cpNumCur;
					  current=new int[numCur];
					  index=0;
					  for(i=0;i<M;i++) 
					  {
						  if(tmp[i]==1) 
						  {
							  current[index]=i; 
							  index++;
						  }
					  }
					  continue LOOP;
				  }
				  else 
				  {
					  System.out.println("Should not come to this point...");
					  System.out.println("Error in BigGirth::selectParityConnect()");
					  return(-1);
				  }
			  }
	  }

	  private void updateConnection(int kthSymbol)
	  {
		  int i, j, m;
		  int[] tmp;

		  for(i=0;i<nodesInGraph[kthSymbol].numOfConnectionSymbolBit;i++)
		  {
			  m=nodesInGraph[kthSymbol].connectionSymbolBit[i];//m [0, M) parity node
			  tmp=new int[nodesInGraph[m].numOfConnectionParityBit+1];
			  for(j=0;j<nodesInGraph[m].numOfConnectionParityBit;j++)
				  tmp[j]=nodesInGraph[m].connectionParityBit[j];
			  tmp[nodesInGraph[m].numOfConnectionParityBit]=kthSymbol;

			  nodesInGraph[m].connectionParityBit=null;
			  nodesInGraph[m].numOfConnectionParityBit++; //increase by 1
			  nodesInGraph[m].connectionParityBit=new int[nodesInGraph[m].numOfConnectionParityBit];
			  for(j=0;j<nodesInGraph[m].numOfConnectionParityBit;j++)
				  nodesInGraph[m].connectionParityBit[j]=tmp[j];
			  tmp=null;
		  }
	  }
	  void loadH()
	  {
		  int i, j;
		  if(H==null) {
			  H=new int[M][];
			  for(i=0;i<M;i++) 
				  H[i]=new int[N];
		  }

		  for(i=0;i<M;i++){
			  for(j=0;j<N;j++){
				  H[i][j]=0;
			  }
		  }
		  for(i=0;i<M;i++){
			  for(j=0;j<nodesInGraph[i].numOfConnectionParityBit;j++)
			  {
				  int k = nodesInGraph[i].connectionParityBit[j];
				  H[M-1-i][N-1-k]=1;
			  }
		  }
	  }
	  int[][] getMask()
	  {
		  return H;
	  }
	  BigGirth(int M, int N, Degree deg, int tgtGirth)
	  {
		  int i, j, k, m, index;
		  int[] mid;
		  int[] localDepth= new int[1];
		  localDepth[0]=100;
		  int[] symbolDegSequence = deg.degSeq(N);
		  EXPAND_DEPTH=(tgtGirth-4)/2; 
		  if(EXPAND_DEPTH<0) EXPAND_DEPTH=0;

		  myrandom=new Random(); 

		  this.M=M;
		  this.N=N;
		  H=null;
		  mid=new int[M];

		  localGirth=new int[N];

		  nodesInGraph=new NodesInGraph [N];
		  for(i=0;i<N;i++)
		  {	
			  nodesInGraph[i] = new NodesInGraph();
			  nodesInGraph[i].setNumOfConnectionSymbolBit(symbolDegSequence[i]);
		  }
		  j=0;
		  for(k=0;k<N;k++) 
			  j+=symbolDegSequence[k];
		  k=j/M;
		  for(i=0;i<M;i++) 
			  mid[i]=k;
		  for(i=0;i<j-k*M;i++) 
			  mid[i]++;
		  k=0; 
		  for(i=0;i<M;i++) 
			  k+=mid[i];
		  if(k!=j) 
		  {
			  System.out.println("Wrong in computing maxDegParity!");
		  }

		  for(i=0;i<M;i++) 
		  {
			  nodesInGraph[i].initConnectionParityBit();
		  } 

		  for(k=0;k<N;k++)
		  {
			  m=1000000;
			  index=-1;
			  if( k<M-1 )
			  {
				  index = k+1;
			  }
			  else
			  {
				  for(i=0;i<M;i++)
				  {
					  if( nodesInGraph[i].numOfConnectionParityBit<m && 
							  nodesInGraph[i].numOfConnectionParityBit<nodesInGraph[i].maxDegParity
						  	) 
					  {
						  m=nodesInGraph[i].numOfConnectionParityBit;
						  index=i;
					  }
				  }
			  }
			  nodesInGraph[k].connectionSymbolBit[0]=index;//least connections of parity bit

			  int iter=0; 
			  ITER:
				  while(true)
				  {
					  localGirth[k]=100;
					  for(m=1;m<nodesInGraph[k].numOfConnectionSymbolBit;m++)
					  {
						  nodesInGraph[k].connectionSymbolBit[m]=selectParityConnectTr(k, m, localDepth);
						  localGirth[k]=(localGirth[k]>localDepth[0])?localDepth[0]:localGirth[k];      
						  if(k>0 && localGirth[k]<localGirth[k-1] && iter<20) 
						  {
							  iter++; 
							  continue ITER;
						  }
					  }
					  if(localGirth[k]==0 && iter<30) 
					  {
						  iter++; 
						  continue ITER;
					  }
					  break ITER;
				  }
		  
			  System.out.print("k="+k+"  ");
			  for(m=0;m<nodesInGraph[k].numOfConnectionSymbolBit;m++)
				  System.out.print(nodesInGraph[k].connectionSymbolBit[m]+" ");
			  System.out.print("LocalGirth="+(2*localGirth[k]+4));
			  System.out.println();
			  updateConnection(k);
		  }
		  loadH();
	  }
	  
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double df[] = {0.521814,    0.271293,     0.0,    0.206893};
		int dg[] = {   2,            3,            4,          5};
		Degree deg = new Degree( dg,df );
		int N = 52;
		int M = 26;
		int[] degSeq = deg.degSeq(N);
		BigGirth g = new BigGirth( M, N, deg, 80);
		
		int i,j;
		System.out.printf("int[][] mask = {\n");
		for( i=0;i<M;i++ )
		{
			System.out.printf("{ ");
			for( j=0;j<N;j++ )
				System.out.printf("%4d, ",g.H[i][j]);
			System.out.printf("},\n");
		}
		System.out.printf("};\n");
	}

}
