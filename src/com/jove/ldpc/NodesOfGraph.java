package com.jove.ldpc;

public class NodesOfGraph {
	  int numOfParityConnections;
	  int[] parityConnections;
	  int numOfSymbolConnections;
	  int[] symbolConnections;
	  int numOfSymbolMapping;
	  int[] symbolMapping;
	  NodesOfGraph()
	  {
		  parityConnections=null;
		  symbolConnections=null;
	  }
	  void setParityConnections(int num, int[] value)
	  {
		  numOfParityConnections=num;
		  parityConnections=new int[num];
		  for(int i=0;i<numOfParityConnections;i++)
		  {
		    parityConnections[i]=value[i];
		  }
	  }
	  void setSymbolConnections(int num, int[] value)
	  {
		  numOfSymbolConnections=num;
		  symbolConnections=new int[num];
		  for(int i=0;i<numOfSymbolConnections;i++)
		  {
		    symbolConnections[i]=value[i];
		  }
	  }
	  void setSymbolMapping(int num, int[] value)
	  {
		  numOfSymbolMapping=num;
		  symbolMapping=new int[num];
		  for(int i=0;i<numOfSymbolMapping;i++)
		  {
		    symbolMapping[i]=value[i];
		  }
	  }

}
