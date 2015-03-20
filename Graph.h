#ifndef GRAPH_H
#define GRAPH_H

#include<iostream>
#include<vector>
#include<list>
#include<stdlib.h>
#include<unordered_map>
#include<math.h>
#include<time.h>
#include<chrono>
#include<cstring>
#include<fstream>


using namespace std;
using namespace std::chrono;

#define MIN(a,b) (((a)<(b))?(a):(b))

#define SIZE 5000
#define WTLIMIT 10000
#define DEGREE 6 
#define SPARSE_DEGREE 6

enum color{white,grey,black};

class edge {
public:
  int destination,weight;
  edge(int destination,int weight);
};

class kedge{
public:
  int u,v,weight;
  kedge(int u,int v,int weight);
};

class _Graph{
  vector<list<edge> >AdjList;
  vector<kedge> EdgeList;
  int rank[SIZE],dad[SIZE];
  color visited[SIZE]={white};
  list<edge> fringeList;

public:
  _Graph();
  void CreateSparse();
  void CreateDense();
  void AddEdge(int source,int destination,int weight);
  void PrintGraph(string filename);
  int MaxCapPath(int start,int end);
  int Kruskal(int start,int end);
  int MaxCapPath_without_Heap(int start,int end);
  void Union(int r1,int r2);
  int FindDad(int vertex);
  bool DFS(int vertex,int end,int *cap);
  void KruskalInit();

  void InsertIntoSortedList(edge E);
  edge SortedMax();
  void SortedDelete(edge E);
  void ClearGraph();
};


// Class for MaxHeap

class MaxHeap{
  vector<edge> H;
  //vector<int> positions;
  int positions[SIZE];
  vector<kedge> KH;

public:
  edge GetMaximum();
  void Insert(edge F);
  void Delete(int index);
  void AddFringe(edge F);
  int HeapSize();
  void Heapify(int index,int direction); // true is up
  void HeapSort();
  void PrintHeap();
  void PrintPositions();
  int GetPosition(int vertex);
  void UpdateWt(int,int);
  kedge KGetMaximum();
  void KInsert(kedge F);
  void KDelete(int index);
  void KAddFringe(kedge F);
  int KHeapSize();
  void KHeapify(int index,int direction); // true is up
  void KPrintHeap();
  void Swap(int index1,int index2);

  void Clear();
  MaxHeap();
  //~MaxHeap();
  
};

void swap(edge*,edge*);
void swap(int,int);
void swap(kedge*,kedge*);


#endif
