package com.jove.ldpc;

public class NodesInGraph {

	  int numOfConnectionParityBit;
	  int[] connectionParityBit;
	  int numOfConnectionSymbolBit;
	  int[] connectionSymbolBit;
	  int maxDegParity;

	  NodesInGraph()
	  {
			connectionParityBit=null;
			connectionSymbolBit=null;
	  }
	  void setNumOfConnectionSymbolBit(int deg)
	  {
		  if(deg<=0) 
		  {
			  System.out.println("Wrong NodesInGraph::setNumOfConnectionSymbolBit()");
		  }
		  numOfConnectionSymbolBit=deg;
		  connectionSymbolBit=new int[deg];
	  }
	  void initConnectionParityBit()
	  {
		  maxDegParity=10000;
		  numOfConnectionParityBit=0;
		  connectionParityBit=new int[1];
	  }
	  void initConnectionParityBit(int deg)
	  {
		  maxDegParity=deg;
		  numOfConnectionParityBit=0;
		  connectionParityBit=new int[1];
	  }
	  
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
