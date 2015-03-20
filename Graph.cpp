// Max Bandwidth Path in a network using Dijkstra and Kruskal
#include "Graph.h"
using namespace std;


edge::edge(int destination,int weight){
  this->destination = destination;
  this->weight = weight;
}

kedge::kedge(int u,int v,int weight):u(u),v(v),weight(weight){}

_Graph::_Graph(){
  int size = SIZE;
  //cout<<"\nGenerating Graph for "<<size<<" vertices"<<endl;
  for (int i=0; i<SIZE; i++) {
    list<edge> E;
    this->AdjList.push_back(E);
    dad[i]=i;
  }
}

void _Graph::ClearGraph(){

  fringeList.clear();
  EdgeList.clear();
  //delete visited;
  //delete rank,dad;
  for (int i=0; i<SIZE; i++) {
    AdjList[i].clear();
  }
  AdjList.clear();
}


void _Graph::AddEdge(int source,int destination,int weight){  // undirected edges
  (AdjList[source]).push_back(edge(destination,weight));
  (AdjList[destination]).push_back(edge(source,weight));
  EdgeList.push_back(kedge(source,destination,weight));
}
void _Graph::PrintGraph(string filename){
  list<edge>::iterator iter;
  ofstream myfile;
  myfile.open(filename);
  for (int i=0; i<SIZE; i++) {
    //cout <<endl;
    //cout << i<<"--->";
    myfile << i<<"--->"<<AdjList[i].size()<<endl;
   
    for (iter=AdjList[i].begin(); iter!=AdjList[i].end(); ++iter) {
      myfile <<"["<<(*iter).destination<<","<<(*iter).weight<<"]"<<"  ";
    }
    myfile<<endl;
  }
  myfile.flush();
  //cout <<endl;
  myfile.close();
}

void _Graph::CreateSparse(){

  unordered_map<int,bool> mp;
  list<edge>::iterator iter;
  int dest;
  ofstream myfile;
  myfile.open("temp.txt");
  for (int i=0; i < SIZE; i++) {

    mp.clear();
    mp.insert({i,true});
    myfile <<"["<< i<<"]-->";
    for (iter=AdjList[i].begin(); iter!=AdjList[i].end(); ++iter) {
      mp.insert({(*iter).destination,true});
      myfile << (*iter).destination<<",";
    }

    //if(AdjList[i].size() >= SPARSE_DEGREE) continue;
    int c = 0;
    while(mp.size() <= SPARSE_DEGREE){
      dest = rand()%SIZE;
      //Check if dest has enough vertices
      c++;
      if((AdjList[dest].size() >= SPARSE_DEGREE) && (c < 500)) {
	    	continue;
      }
      if(mp.find(dest) == mp.end()){
	      mp.insert({dest,true});
	      AddEdge(i,dest,(1 + rand()%WTLIMIT));
	      myfile << dest<<",";
      }
    }
    myfile <<endl;
    myfile.flush();
    //cout <<"["<< i<<"]";
  }
  myfile.close();
  int degree = SPARSE_DEGREE;
  cout << "Generated Sparse Graph with degree "<<degree<<endl;
}

void _Graph::CreateDense(){

  unordered_map<int,bool> mp;
  int dest;
  int twentyp = (int)((0.20)*SIZE);
  int fifteenp = (int)((0.14)*SIZE);

  list<edge>::iterator iter;
  cout << "Twenty percent = "<<twentyp<<endl;
  for (int i=0; i<SIZE; i++) {

    if(AdjList[i].size() >= fifteenp) continue;
    mp.clear();
    mp.insert({i,true});
    
    for (iter=AdjList[i].begin(); iter!=AdjList[i].end(); ++iter) {
      mp.insert({(*iter).destination,true});
    }
    
    while(mp.size() <= fifteenp){
      dest = rand()%SIZE;
      //Check if dest has enough vertices
      if(AdjList[dest].size() >= twentyp) continue;

      if(mp.find(dest) == mp.end()){
	      mp.insert({dest,true});
	      AddEdge(i,dest,(1 + rand()%WTLIMIT));
      }
    }
    //cout<<"Done for Vertex "<<i<<endl;
    if(i==1000 || i== 2000 || i==3000||i==4000)
      cout << i<<",";
  }
  cout <<"Created dense graph with approx degree 1000\n";
}

void _Graph::KruskalInit(){

  for (int i=0; i<SIZE; i++) {
    dad[i]=i;rank[i]=0;
    visited[i]=white;
  }

}


int _Graph::Kruskal(int start,int end){

  // cout <<"Kruskal Starts"<<endl;
  KruskalInit();
  MaxHeap Heap;
  _Graph MST;
  //clock_gettime(CLOCK_MONOTONIC, &tstart);
  //Create heap from edges
  //MaxHeap Heap;
  struct timespec tstart={0,0}, tend={0,0};
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  for (int i=0; i<EdgeList.size(); i++) {
    Heap.KInsert(EdgeList[i]);
  }
  //Heap.KPrintHeap();
  

  
  while(Heap.KHeapSize()){
    
    kedge maxedge = Heap.KGetMaximum();
    Heap.KDelete(0);

    //cout <<" Selected [" << maxedge.u << "," <<maxedge.v << "]\n";
    int udad,vdad;
    udad = FindDad(maxedge.u);
    vdad = FindDad(maxedge.v);

    // cout <<"Parents of "<< maxedge.u <<" is "<<udad<<endl;
    //cout <<"Parents of "<< maxedge.v <<" is "<<vdad<<endl;

    if(vdad != udad){
      MST.AddEdge(maxedge.u,maxedge.v,maxedge.weight);
      Union(udad,vdad);
      
    }
  }
  
  //MST.PrintGraph("MST.txt");
  
  //do DFS from start to end
  int maxcapacity = 999999;
  MST.DFS(start,end,&maxcapacity);
  //cout <<"Kruskal Ends"<<endl;
  //clock_gettime(CLOCK_MONOTONIC, &tend);

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  cout << "Kruskal with Heap: \t Time taken = "<< duration<<"\t MaxCap = "<< maxcapacity<<endl;
  return maxcapacity;
  
}


int _Graph::MaxCapPath(int start,int end){

  vector<int> dist(SIZE,0);
  vector<int> dad(SIZE,-1);
  vector<int> maxcapacity(SIZE,WTLIMIT+2);
  enum status {unseen,fringe,intree};
  vector<status> S(SIZE,unseen);
  //struct timespec tstart={0,0}, tend={0,0};
  //clock_gettime(CLOCK_MONOTONIC, &tstart);

  //cout <<"Dijkstra with Heap Starts"<<endl;
  MaxHeap Heap;
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  for (auto it=AdjList[start].begin(); it!=AdjList[start].end(); ++it) {
    dist[(*it).destination] = (*it).weight;
    maxcapacity[(*it).destination] = (*it).weight;
    S[(*it).destination] = fringe;
    Heap.Insert(*it);
    dad[(*it).destination] = start;
  }
  S[start] = intree;
  dad[start] = start;

  while(Heap.HeapSize() >0){
    //Heap.PrintHeap();
    edge maxedge = Heap.GetMaximum();
    Heap.Delete(0);

    int u = maxedge.destination;
    //cout << "\tSelected vertex "<<u<<" weight "<<maxedge.weight<<endl;
    S[u] = intree;

    if(u == end) break;

    for (auto iter=AdjList[u].begin(); iter!=AdjList[u].end(); ++iter) {
      int v = (*iter).destination;
      int wt = (*iter).weight;
      if(v == u) continue; // dont look at parent
      if(S[v] == unseen){
      	S[v] = fringe;
        dad[v] = u;
      	//dist[v] = dist[u] + wt;
      	dist[v] = wt;
      	maxcapacity[v] = MIN(maxcapacity[u],wt);
      	Heap.Insert(edge(v,maxcapacity[v]));
      }
      else if((S[v] == fringe) && (MIN(maxcapacity[u],wt) > maxcapacity[v])){// (dist[v] < wt)){//dist[u]+wt)){
        //cout <<"U";

      	//dist[v] = dist[u] + wt;
      	dist[v] = wt;
      	maxcapacity[v] = MIN(maxcapacity[u],wt);
      	int heapindex = Heap.GetPosition(v);
        {
          Heap.UpdateWt(heapindex,maxcapacity[v]);
          Heap.Heapify(heapindex,1);
        }
	      //Heap.Delete(heapindex);
	      //Heap.Insert(edge(v,maxcapacity[v]));
        dad[v] = u;
      }
    }
  }
  //cout <<"Distance to vertex "<<end<<" is " << dist[end]<<endl;
  //delete dist;
  //Heap.Clear();
  //S.clear();
  //cout <<"Completed Dijkstra "<<maxcapacity[end]<<endl;
  //clock_gettime(CLOCK_MONOTONIC, &tend);


  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  cout << "Dijkstra with Heap: \t Time taken = "<< duration<<"\t MaxCap = "<< maxcapacity[end]<<"  ";
  int path = end;
  cout <<endl;
  /*do{
    cout <<"["<<path<<","<<maxcapacity[path]<<"]<--";
    path = dad[path];
  }while(dad[path] != path);
  cout << start<<endl;
  
  for(int i =0 ;i<SIZE;i++){
    cout <<dad[i]<<",";
  }*/

  return maxcapacity[end];
}


int _Graph::MaxCapPath_without_Heap(int start,int end){

  vector<int> dist(SIZE,0);
  vector<int> dad(SIZE,-1);
  vector<int> maxcapacity(SIZE,WTLIMIT+2);
  enum status {unseen,fringe,intree};
  vector<status> S(SIZE,unseen);
  struct timespec tstart={0,0}, tend={0,0};

  //clock_gettime(CLOCK_MONOTONIC, &tstart);

  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  
  //cout <<"Dijkstra without Heap Starts"<<endl;
  for (auto it=AdjList[start].begin(); it!=AdjList[start].end(); ++it) {
    dist[(*it).destination] = (*it).weight;
    maxcapacity[(*it).destination] = (*it).weight;
    S[(*it).destination] = fringe;
    fringeList.push_back(*it);
    dad[(*it).destination] = start;
    
  }
  S[start] = intree;
  dad[start] = start;
  while(fringeList.size() > 0){
    
    // get max from the fringe list
    list<edge>::iterator maxit = fringeList.begin();
    for (auto myit=fringeList.begin(); myit!=fringeList.end(); ++myit) {
      if((*maxit).weight < (*myit).weight){
      	maxit = myit;
      }
    }
    int u = (*maxit).destination;
    fringeList.erase(maxit);
    //cout << "Selected vertex "<<u<<endl;
    S[u] = intree;

    //if(u == end) break;

    for (auto iter=AdjList[u].begin(); iter!=AdjList[u].end(); ++iter) {
      int v = (*iter).destination;
      int wt = (*iter).weight;
      if(v == u) continue; // dont look at parent
      if(S[v] == unseen){
	      S[v] = fringe;
        dad[v] = u;
	      dist[v] = wt;
	      maxcapacity[v]= MIN(maxcapacity[u],wt);
	      fringeList.push_back(edge(v,maxcapacity[v]));
      }
      else if((S[v] == fringe) && (MIN(maxcapacity[u],wt) > maxcapacity[v])){//(dist[v] < wt)){//dist[u]+wt)){
	      list<edge>::iterator del;
	      for (auto myit=fringeList.begin(); myit!=fringeList.end(); ++myit) {
	        if(v == (*myit).destination){
	          del = myit;
	        }
	      }
	      fringeList.erase(del);
	      dist[v] = wt;
	      dad[v] = u;
	      maxcapacity[v]= MIN(maxcapacity[u],wt);
	      fringeList.push_back(edge(v,maxcapacity[v]));
      }
    }
  }
  //cout <<"Dijkstra without Heap Ends"<<endl;
  /*  clock_gettime(CLOCK_MONOTONIC, &tend);
  double timetaken = 
    ((tend.tv_sec - tstart.tv_sec)*10^9) -
    (tend.tv_nsec - tstart.tv_nsec);
  */

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  cout << "Dijkstra without Heap: \t Time taken = "<<duration <<"\t MaxCap = "<< maxcapacity[end]<<"  ";
  int path = end,pathlength =0;
  cout <<endl;



  /*do{
    cout <<"["<<path<<","<<maxcapacity[path]<<"]<--";
    path = dad[path];
    if(pathlength++ > 50) break;
  }while(dad[path] != path);
  cout << start<<endl;
  cout <<endl;
  for(int i =0 ;i<SIZE;i++){
    cout <<dad[i]<<",";
  }*/
  cout <<endl;
  return maxcapacity[end];
}

void TestGraphs(bool sparse){

  int i=0;
  int dheap,kr,djwheap;
  srand(time(NULL));
  while(i++<5){
    _Graph G;
    if(sparse == true)
      G.CreateSparse();
    else
      G.CreateDense();

    G.PrintGraph("Output.txt");
    int j=0;
    while (j++<5){
      
      int rx = rand()%SIZE;
      int ry = rand()%SIZE;

      if(rx == ry) {ry = rand()%SIZE;}

      cout <<"Random start and end vertices : "<<rx<<" "<<ry<<endl;

      //G.AddEdge(rx,ry,1);
      dheap = G.MaxCapPath(rx,ry);
      kr = G.Kruskal(rx,ry);
      djwheap = G.MaxCapPath_without_Heap(rx,ry);

      if(dheap != kr || kr != djwheap || dheap != djwheap)
      {
        cout <<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Mismatch in algo>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
        exit(0);
      }
    }
    G.ClearGraph();
  }
}


int main(int argc, char*argv[]){
  
  //cout <<argc;
  if(argc == 2){
    //cout <<argv[1]; 
    if(strcmp(argv[1],"1") == 0 ){
      TestGraphs(true);
      return 0;
    }
    else if (strcmp(argv[1],"2") == 0){
      TestGraphs(false);
      return 0;
    }
    else if (strcmp(argv[1],"3") == 0){
      _Graph G;
      G.CreateSparse();
      G.PrintGraph("Output.txt");

      int x,y;
      while(1){
        cout <<"\nEnter start and end vertices for Max Cap path: ";
        cin >> x>>y;
        G.MaxCapPath(x,y);
        G.Kruskal(x,y);
        G.MaxCapPath_without_Heap(x,y);
      }

      return 0;
    }
    else if (strcmp(argv[1],"4") == 0){
      _Graph G;
      G.CreateDense();
      G.PrintGraph("Output.txt");
      
      int x,y;
      while(1){
        cout <<"\nEnter start and end vertices for Max Cap path: ";
        cin >> x>>y;
       
        G.MaxCapPath(x,y);
        G.Kruskal(x,y);
        G.MaxCapPath_without_Heap(x,y);
      }
      return 0;
    }
  }

  cout << "Usage: a.out <option> where option is {1,2,3,4}\n"
       << "1 for sparse graph random test\n"
       << "2 for dense graph random test\n"
       << "3 to manually enter start and destination vertices for sparse graph\n"
       << "4 to manually enter start and destination vertices for dense graph\n"<<endl;
  return 0;

  _Graph G;
  G.CreateSparse();
  //G.CreateDense();
  G.PrintGraph("Output.txt");
  
  int x,y;
  while(1){
    cout <<"\nEnter start and end vertices for Max Cap path:(Ctrl^C to exit) ";
    cin >> x>>y;
    //G.AddEdge(x,y,1);
    
    G.MaxCapPath(x,y);
    G.Kruskal(x,y);
    G.MaxCapPath_without_Heap(x,y);
  }

  /*
  MaxHeap Heap;
  int i=0;
  while(i++<20){
    Heap.Insert(edge(i,(1+rand()%30)));
    Heap.PrintHeap();
  }
  Heap.HeapSort();*/
}


/*****************/
//Heap functions

/*edge::edge(int index,int weight){
  this->index = index;
  this->weight = weight;''
  }*/

MaxHeap::MaxHeap(){
  //positions.assign(SIZE,-1);
  for (int i=0; i<SIZE; i++) {
    positions[i] = -1;
  }

}

/*MaxHeap::~MaxHeap(){
  H.clear();
  positions.clear();
  KH.clear();
}*/

void MaxHeap::UpdateWt(int index,int wt){
  H[index].weight = wt;
}

edge MaxHeap::GetMaximum(){ 
  if(H.size() ==0)
    cout<<"WRONG!!";
  return H[0];
}
void MaxHeap::AddFringe(edge F){
  //cout << F.destination<<endl;
  H.push_back(F);
  positions[F.destination] = H.size()-1;
  //cout <<"-->>";PrintHeap();
}
int MaxHeap::HeapSize(){return H.size();}

int MaxHeap::GetPosition(int vertex){return positions[vertex];}

void MaxHeap::PrintHeap(){
  for (int i=0; i<H.size(); i++) {
    cout<<"["<<H[i].destination<<","<< H[i].weight<<","
	<<"p["<<H[i].destination<<"]="<<
      positions[H[i].destination]<<"="<<i<<"],";
  }
  cout << endl;
}

void MaxHeap::PrintPositions(){
  cout <<"Positions "<<endl;
  for (int i=0; i<SIZE; i++) {
    cout<< positions[i]<<"\t";
  }
  cout << endl;
}


void MaxHeap::Heapify(int index,int direction){
  
  if(direction == 1){
    if(index == 0) return;

    int parent = ((index+1)>>1)-1;
    if(H[parent].weight < H[index].weight){
      
      swap(&H[parent],&H[index]);
      
      positions[H[parent].destination] = parent;
      positions[H[index].destination] = index;
      
      if(positions[H[parent].destination] == -1 ||
      	 positions[H[index].destination] == -1)
        	cout<<"XSHIT";

      Heapify(parent,1);
    }
    else return;
  }
  else{

    if(H.size() == 0) return;

    int child1 = (index<<1)+1;
    if(child1 > H.size()-1) return; // no valid children
    int child2 = child1+1;

    if(child2 > H.size()-1){ //only one valid child

      if(H[child1].weight > H[index].weight){

      	swap(&H[child1],&H[index]);
      	positions[H[child1].destination] = child1;
      	positions[H[index].destination] = index;	

    	  if(positions[H[child1].destination] == -1 ||
    	     positions[H[index].destination] == -1)
    	    cout<<"SHIT";
	
    	  Heapify(child1,0);
      }
      else return;
    }
    else { // both children are valid
      int mcindex;
      if(H[child2].weight < H[child1].weight) mcindex = child1;
      else mcindex = child2;

      if(H[mcindex].weight > H[index].weight){
      	swap(&H[mcindex],&H[index]);

      	positions[H[mcindex].destination] = mcindex;
      	positions[H[index].destination] = index;

      	if(positions[H[mcindex].destination] == -1 ||
      	   positions[H[index].destination] == -1)
      	  cout<<"SHIT";

      	Heapify(mcindex,0);
      }
      else return;
    }
  }
  for (int i=0; i<SIZE; i++) {
    if((positions[i] <=0) || positions[i] > H.size()-1)
      cout<"FXK";
  }


}
void MaxHeap::Insert(edge x){
  
  this->AddFringe(x);
  this->Heapify(H.size()-1,1);
}
void MaxHeap::Delete(int index){
  
  if(index <0 || index >= H.size()){
    cout << "Wrong index passed for deletion from heap" << index<<endl;
    return;
  }

  positions[H[index].destination] = -1;
  swap(&H[index],&H[H.size()-1]);
  positions[H[index].destination] = index;

  H.pop_back();
  this->Heapify(index,0);

}
 
void MaxHeap::Swap(int index1,int index2){
  edge temp = H[index1];
  H[index1] = H[index2];
  H[index2] = temp;
}

void MaxHeap::HeapSort(){
  if(H.size()==0) {cout<<endl;return;}
  cout<<(GetMaximum()).weight<<',';
  Delete(0);
  HeapSort();

}


/*****************/
//Kruskal Heap 

kedge MaxHeap::KGetMaximum(){ 
  if(KH.size() ==0)
    cout<<"WRONG!!";
  return KH[0];
}
void MaxHeap::KAddFringe(kedge F){KH.push_back(F);}
int MaxHeap::KHeapSize(){return KH.size();}
void MaxHeap::KPrintHeap(){
  for (int i=0; i<KH.size(); i++) {
    cout<<"["<<KH[i].u<<","<<KH[i].v<<"] -->"<< KH[i].weight<<"\n";
  }
}

void MaxHeap::KHeapify(int index,int direction){

  
  if(direction == 1){
    if(index == 0) return;
    int parent = ((index+1)>>1)-1;
    if(KH[parent].weight < KH[index].weight){
      swap(&KH[parent],&KH[index]);
      KHeapify(parent,1);
    }
    else return;
  }
  else{

    if(KH.size() == 0){ return;}
    int child1 = (index<<1)+1;
    if(child1 > KH.size()-1) return;
    int child2 = (index<<1)+2;

    if(child2 > KH.size()-1){ //only one valid child
      if(KH[child1].weight > KH[index].weight){
	      swap(&KH[child1],&KH[index]);
	      KHeapify(child1,0);
      }
      else return;
    }
    else { // both children are valid
      int mcindex;
      if(KH[child2].weight < KH[child1].weight) mcindex = child1;
      else mcindex = child2;

      if(KH[mcindex].weight > KH[index].weight){
	      swap(&KH[mcindex],&KH[index]);
	      KHeapify(mcindex,0);
      }
      else return;
    }
  }
}
void MaxHeap::KInsert(kedge x){
  this->KAddFringe(x);
  this->KHeapify(KH.size()-1,1);
  
}

void MaxHeap::KDelete(int index){
  swap(&KH[index],&KH[KH.size()-1]);
  KH.pop_back();
  this->KHeapify(index,0);
  
}
 
void _Graph::Union(int r1,int r2){
  if(rank[r1]>rank[r2]) dad[r2] = r1;
  else if(rank[r2]>rank[r1]) dad[r1] = r2;
  else { dad[r1] = r2; rank[r2]++;}
}

int _Graph::FindDad(int vertex){
  int x = vertex;
  //vector<int> updateDad;
  while(dad[x] != x){
    //updateDad.push_back(x);
    x = dad[x];
  }
  // Path Compression
  /*for (int i=0; i<updateDad.size(); i++) {
    dad[updateDad[i]] = x;
    }*/

  return x;
}
  
bool _Graph::DFS(int vertex,int end,int *maxcapacity){
  
  if(vertex == end) return true;

  visited[vertex] = grey;
  for (auto myit = AdjList[vertex].begin(); myit!=AdjList[vertex].end(); ++myit) {
    
    if(visited[(*myit).destination] == white){
      int cap = MIN(*maxcapacity,(*myit).weight);
      if(DFS((*myit).destination,end,&cap) == true){
	*maxcapacity = cap;
	return true;
      }
    }
  }
  visited[vertex] = black;
  return false;
}

void MaxHeap::Clear(){
  H.clear();
  KH.clear();
}
void swap(edge *x,edge *y){

  edge temp = *x;
  *x = *y;
  *y = temp;
}
void swap(int *x,int*y){
  int temp = *x;
  *x = *y;
  *y = temp;
}
void swap(kedge *x,kedge *y){

  kedge temp = *x;
  *x = *y;
  *y = temp;
}
